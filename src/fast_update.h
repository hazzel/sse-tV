#pragma once
#include <vector>
#include <array>
#include <iostream>
#include <Eigen/Dense>
#include "dump.h"
#include "lattice.h"
#include "parameters.h"
#include "Random.h"

class fast_update
{
	public:
		template<int n, int m>
		using matrix_t = Eigen::Matrix<double, n, m, Eigen::ColMajor>; 
		using dmatrix_t = matrix_t<Eigen::Dynamic, Eigen::Dynamic>;

		fast_update(Random& rng_, const lattice& l_, const parameters& param_)
			: rng(rng_), l(l_), param(param_)
		{}

		void initialize()
		{
			for (int i = 0; i < l.n_sites(); ++i)
				for (auto j : l.neighbors(i, "nearest neighbors"))
					if (i < j)
						lattice_bonds.push_back({i, j});
			max_order(10);
			n_non_ident = 0;
			equal_time_gf = 0.5 * dmatrix_t::Identity(l.n_sites(), l.n_sites());
			id_2 = matrix_t<2, 2>::Identity();
			id_N = dmatrix_t::Identity(l.n_sites(), l.n_sites());
			dmatrix_t expLam = vertex_matrix(param.lambda);
			dmatrix_t invExpLam = vertex_matrix(-param.lambda);
			A = vertex_block(expLam);
			invA = vertex_block(invExpLam);
			B = A - id_2;
			invB = B.inverse();
			C = invA - id_2;
			invC = C.inverse();
		}

		void max_order(int n_max_order_)
		{
			n_max_order = n_max_order_;
			if (n_max_order_ % param.n_delta != 0)
				n_max_order += param.n_delta - n_max_order_ % param.n_delta;
			current_vertex = n_max_order;
			n_intervals = n_max_order / param.n_delta;
			bond_list.resize(n_max_order, 0);
			U.resize(n_intervals + 1);
			D.resize(n_intervals + 1);
			V.resize(n_intervals + 1);
			for (int n = 0; n < n_intervals + 1; ++n)
			{
				U[n] = id_N; D[n] = id_N; V[n] = id_N;
			}
			std::cout << "max order set to " << n_max_order << std::endl
				<< "n_intervals set to " << n_intervals << std::endl;
		}
		int max_order() const
		{
			return n_max_order;
		}

		int non_ident() const
		{
			return n_non_ident;
		}

		int get_current_bond()
		{
			return bond_list[current_vertex-1];
		}

		void rebuild()
		{
			std::cout << "start rebuild" << std::endl;
			for (int n = 1; n <= n_intervals; ++n)
				store_svd_forward(n - 1);
			std::cout << "rebuild done" << std::endl;
		}

		void serialize(odump& out)
		{
		}

		void serialize(idump& in)
		{
		}

		dmatrix_t vertex_matrix(double lambda, int vertex_id)
		{
			dmatrix_t lam = dmatrix_t::Identity(l.n_sites(), l.n_sites());
			int bond_id = bond_list[vertex_id-1] - 1;
			if (bond_id < 0)
				return lam;
			std::pair<int, int> bond = lattice_bonds[bond_id];
			if (lambda > 0.)
			{
				lam(bond.first, bond.first) = A(0, 0);
				lam(bond.first, bond.second) = A(0, 1);
				lam(bond.second, bond.first) = A(1, 0);
				lam(bond.second, bond.second) = A(1, 1);
			}
			else
			{
				lam(bond.first, bond.first) = invA(0, 0);
				lam(bond.first, bond.second) = invA(0, 1);
				lam(bond.second, bond.first) = invA(1, 0);
				lam(bond.second, bond.second) = invA(1, 1);
			}
			return lam;
		}

		dmatrix_t vertex_matrix(double lambda)
		{
			dmatrix_t lam = dmatrix_t::Zero(l.n_sites(), l.n_sites());
			std::pair<int, int> bond = lattice_bonds[0];
			lam(bond.first, bond.second) = lambda;
			lam(bond.second, bond.first) = lambda;
			Eigen::SelfAdjointEigenSolver<dmatrix_t> solver(lam);
			dmatrix_t d = solver.eigenvalues().unaryExpr([](double e)
				{ return std::exp(e); }).asDiagonal();
			return solver.eigenvectors() * d * solver.eigenvectors().adjoint();
		}
		
		matrix_t<2, 2> vertex_block(const dmatrix_t& m, int vertex_id)
		{
			int bond_id = bond_list[vertex_id-1] - 1;
			if (bond_id < 0)
				return matrix_t<2, 2>::Identity();
			std::pair<int, int> bond = lattice_bonds[bond_id];
			matrix_t<2, 2> block;
			block << m(bond.first, bond.first), m(bond.first, bond.second),
				m(bond.second, bond.first), m(bond.second, bond.second);
			return block;
		}
		
		matrix_t<2, 2> vertex_block(const dmatrix_t& m,
			const std::pair<int, int>& bond)
		{
			matrix_t<2, 2> block;
			block << m(bond.first, bond.first), m(bond.first, bond.second),
				m(bond.second, bond.first), m(bond.second, bond.second);
			return block;
		}

		matrix_t<2, 2> vertex_block(const dmatrix_t& m)
		{
			std::pair<int, int> bond = lattice_bonds[0];
			matrix_t<2, 2> block(2, 2);
			block << m(bond.first, bond.first), m(bond.first, bond.second),
				m(bond.second, bond.first), m(bond.second, bond.second);
			return block;
		}

		dmatrix_t get_R()
		{
			dmatrix_t R = dmatrix_t::Identity(l.n_sites(), l.n_sites());
			for (int n = current_vertex; n >= 1; --n)
			{
				if (bond_list[n-1] == 0) continue;
				R *= vertex_matrix(param.lambda, n);
			}
			return R;
		}
		
		dmatrix_t get_L()
		{
			dmatrix_t L = dmatrix_t::Identity(l.n_sites(), l.n_sites());
			for (int n = n_max_order; n > current_vertex; --n)
			{
				if (bond_list[n-1] == 0) continue;
				L *= vertex_matrix(param.lambda, n);
			}
			return L;
		}
		
		dmatrix_t propagator(int n, int m)
		{
			dmatrix_t P = dmatrix_t::Identity(l.n_sites(), l.n_sites());
			for (int i = n; i > m; --i)
			{
				if (bond_list[i-1] == 0) continue;
				P *= vertex_matrix(param.lambda, i);
			}
			return P;
		}

		void print_gf_from_scratch()
		{
			std::cout << "current vertex: " << current_vertex
				<< ", non_ident: " << n_non_ident << std::endl;
			dmatrix_t G = (id_N + get_R()*get_L()).inverse();
			std::cout << "|G-Gp| = " << (G-equal_time_gf).norm() << std::endl;
			std::cout << "differences at" << std::endl;
			for (int i = 0; i < l.n_sites(); ++i)
				for (int j = 0; j < l.n_sites(); ++j)
					if (std::abs(G(i, j) - equal_time_gf(i, j)) > 0.000001)
						std::cout << i << ", " << j << std::endl;
			std::cout << std::endl;
			if ((G-equal_time_gf).norm() > 0.00001)
			{
				print_bonds();
				print_matrix(G);
				print_matrix(equal_time_gf);
			}
			//std::cout << "G from scratch" << std::endl;
			//print_matrix(G);
			//std::cout << "G from update" << std::endl;
			//print_matrix(equal_time_gf);
		}

		void print_bonds()
		{
			std::cout << "print_bonds()" << std::endl;
			for (auto i : bond_list)
				if (i > 0)
					std::cout << lattice_bonds[i-1].first << " , "
						<< lattice_bonds[i-1].second << std::endl;
			std::cout << std::endl;
		}
		
		double try_update_vertex(int bond_type)
		{
			// Insert bond at vertex
			if (get_current_bond() == 0)
			{
				int bond_id = rng() * l.n_bonds();
				bond_buffer = bond_id + 1;
				std::pair<int, int> bond = lattice_bonds[bond_id];
				matrix_t<2, 2> d = id_2 + B * (id_2 - vertex_block(equal_time_gf,
					bond));
				return d.determinant() * l.n_bonds() * param.beta * param.t
					/ ((n_max_order - n_non_ident) * std::sinh(param.lambda));
			}
			// Remove bond at vertex
			else
			{
				bond_buffer = 0;
				matrix_t<2, 2> d = id_2 + C * (id_2 - vertex_block(equal_time_gf,
					current_vertex));
				return d.determinant() * ((n_max_order - n_non_ident + 1.) 
					* std::sinh(param.lambda)) / (l.n_bonds() * param.beta * param.t);
			}
		}

		void finish_update_vertex(int bond_type)
		{
			// Insert bond at vertex
			if (get_current_bond() == 0)
			{
				int bond_id = bond_buffer - 1;
				bond_list[current_vertex-1] = bond_buffer;
				std::pair<int, int> bond = lattice_bonds[bond_id];
				matrix_t<2, 2> d = (invB + (id_2 - vertex_block(equal_time_gf,
					bond))).inverse();
				dmatrix_t e = dmatrix_t::Zero(l.n_sites(), l.n_sites());
				e(bond.first, bond.first) = d(0, 0);
				e(bond.first, bond.second) = d(0, 1);
				e(bond.second, bond.first) = d(1, 0);
				e(bond.second, bond.second) = d(1, 1);
				equal_time_gf = equal_time_gf - (equal_time_gf * e * (id_N
					- equal_time_gf));
				++n_non_ident;
			}
			// Remove bond at vertex
			else
			{
				int bond_id = bond_list[current_vertex-1] - 1;
				bond_list[current_vertex-1] = 0;
				std::pair<int, int> bond = lattice_bonds[bond_id];
				matrix_t<2, 2> d = (invC + (id_2 - vertex_block(equal_time_gf,
					bond))).inverse();
				dmatrix_t e = dmatrix_t::Zero(l.n_sites(), l.n_sites());
				e(bond.first, bond.first) = d(0, 0);
				e(bond.first, bond.second) = d(0, 1);
				e(bond.second, bond.first) = d(1, 0);
				e(bond.second, bond.second) = d(1, 1);
				equal_time_gf = equal_time_gf - (equal_time_gf * e * (id_N
					- equal_time_gf));
				--n_non_ident;
			}
		}

		void start_forward_sweep()
		{
			//recompute_equal_time_gf(id_N, id_N, id_N, V.front(), D.front(),
			//	U.front());
			U.front() = dmatrix_t::Identity(l.n_sites(), l.n_sites());
			D.front() = dmatrix_t::Identity(l.n_sites(), l.n_sites());
			V.front() = dmatrix_t::Identity(l.n_sites(), l.n_sites());
		}

		void start_backward_sweep()
		{
			//recompute_equal_time_gf(U.back(), D.back(), V.back(), id_N, id_N,
			//	id_N);
			U.back() = dmatrix_t::Identity(l.n_sites(), l.n_sites());
			D.back() = dmatrix_t::Identity(l.n_sites(), l.n_sites());
			V.back() = dmatrix_t::Identity(l.n_sites(), l.n_sites());
		}

		void advance_forward()
		{
			if (current_vertex == n_max_order)
				return;
			std::cout << "move forward, current vertex = " << current_vertex
				<< ", non_ident = " << n_non_ident << std::endl;
			// Wrap equal time gf forwards
			equal_time_gf = vertex_matrix(param.lambda, current_vertex + 1)
				* equal_time_gf * vertex_matrix(-param.lambda, current_vertex + 1);
			++current_vertex;
			
			dmatrix_t g = (id_N + get_R() * get_L()).inverse();
			double diff = (equal_time_gf - g).norm();
			if (diff > 0.001)
			{
				print_bonds();
				std::cout << "equal time gf:" << std::endl;
				print_matrix(equal_time_gf);
				std::cout << "correct" << std::endl;
				print_matrix(g);
				std::cout << "diff: " << (equal_time_gf - g).norm() << std::endl;
			}
			// Stabilize equal time gf
			if (current_vertex % param.n_delta == 0)
				store_svd_forward(current_vertex / param.n_delta - 1);
			std::cout << "-----" << std::endl << std::endl;
		}
		
		void advance_backward()
		{
			if (current_vertex == 1)
				return;
			std::cout << "move backward, current vertex = " << current_vertex
				<< ", non_ident = " << n_non_ident << std::endl;
			// Wrap equal time gf backwards
			equal_time_gf = vertex_matrix(-param.lambda, current_vertex)
				* equal_time_gf * vertex_matrix(param.lambda, current_vertex);
			--current_vertex;

			dmatrix_t g = (id_N + get_R() * get_L()).inverse();
			double diff = (equal_time_gf - g).norm();
			if (diff > 0.001)
			{
				std::cout << "equal time gf:" << std::endl;
				print_matrix(equal_time_gf);
				std::cout << "correct" << std::endl;
				print_matrix(g);
				std::cout << "diff: " << (equal_time_gf - g).norm() << std::endl;
			}
			// Stabilize equal time gf
			if (current_vertex % param.n_delta == 0 || current_vertex == 1)
				store_svd_backward(current_vertex / param.n_delta + 1);
			std::cout << "-----" << std::endl << std::endl;
		}

		// n = 0, ..., n_intervals - 1
		void store_svd_forward(int n)
		{
			std::cout << "n = " << n + 1 << std::endl
				<< "propagator from " << (n+1)*param.n_delta
				<< " to " << n * param.n_delta << std::endl;
			dmatrix_t U_l = U[n+1];
			dmatrix_t D_l = D[n+1];
			dmatrix_t V_l = V[n+1];
			dmatrix_t b = propagator((n+1) * param.n_delta, n * param.n_delta);
			svd_solver.compute(b * U[n] * D[n], Eigen::ComputeThinU |
				Eigen::ComputeThinV);
			U[n+1] = svd_solver.matrixU();
			D[n+1] = svd_solver.singularValues().asDiagonal();
			V[n+1] = svd_solver.matrixV().adjoint() * V[n];

			
			dmatrix_t b1 = propagator((n+1) * param.n_delta, 0);
			dmatrix_t b2 = propagator(n_max_order, (n+1) * param.n_delta);
			svd_solver.compute(b2, Eigen::ComputeThinU | Eigen::ComputeThinV);
			V_l = svd_solver.matrixU();
			D_l = svd_solver.singularValues().asDiagonal();
			U_l = svd_solver.matrixV().adjoint();
			svd_solver.compute(b1, Eigen::ComputeThinU | Eigen::ComputeThinV);
			U[n+1] = svd_solver.matrixU();
			D[n+1] = svd_solver.singularValues().asDiagonal();
			V[n+1] = svd_solver.matrixV().adjoint();

			equal_time_gf = (id_N + b1 * b2).inverse();

			/*
			std::cout << "B(beta, (n+1)delta)" << std::endl;
			dmatrix_t t = propagator(n_max_order, (n+1) * param.n_delta);
			print_matrix(t);
			std::cout << "V_l * D_l * U_l" << std::endl;
			print_matrix(V_l * D_l * U_l);
			std::cout << "B((n+1)delta, 0)" << std::endl;
			t = propagator((n+1) * param.n_delta, 0);
			print_matrix(t);
			std::cout << "U * D * V" << std::endl;
			print_matrix(U[n+1] * D[n+1] * V[n+1]);
			*/
			//recompute_equal_time_gf(U_l, D_l, V_l, U[n+1], D[n+1], V[n+1]);
		}
	
		//n = n_intervals, ..., 1 
		void store_svd_backward(int n)
		{
			std::cout << "n = " << n - 1 << std::endl
				<< "propagator from " << n*param.n_delta
				<< " to " << (n-1) * param.n_delta << std::endl;
			dmatrix_t b = propagator(n * param.n_delta, (n-1) * param.n_delta);
			svd_solver.compute(D[n] * U[n] * b, Eigen::ComputeThinU |
				Eigen::ComputeThinV);
			dmatrix_t U_r = U[n-1];
			dmatrix_t D_r = D[n-1];
			dmatrix_t V_r = V[n-1];
			V[n-1] = V[n] * svd_solver.matrixU();
			D[n-1] = svd_solver.singularValues().asDiagonal();
			U[n-1] = svd_solver.matrixV().adjoint();
			
			
			dmatrix_t b1 = propagator(current_vertex, 0);
			dmatrix_t b2 = propagator(n_max_order, current_vertex);
			svd_solver.compute(b2, Eigen::ComputeThinU | Eigen::ComputeThinV);
			U_r = svd_solver.matrixU();
			D_r = svd_solver.singularValues().asDiagonal();
			V_r = svd_solver.matrixV().adjoint();
			svd_solver.compute(b1, Eigen::ComputeThinU | Eigen::ComputeThinV);
			V[n-1] = svd_solver.matrixU();
			D[n-1] = svd_solver.singularValues().asDiagonal();
			U[n-1] = svd_solver.matrixV().adjoint();

			equal_time_gf = (id_N + b1 * b2).inverse();
			
			/*
			std::cout << "B((n-1)delta, 0)" << std::endl;
			dmatrix_t t = propagator((n-1) * param.n_delta, 0);
			print_matrix(t);
			std::cout << "U_r * D_r * V_r" << std::endl;
			print_matrix(U_r * D_r * V_r);
			t = propagator(n_max_order, (n-1) * param.n_delta);
			std::cout << "B(beta, (n-1)delta)" << std::endl;
			print_matrix(t);
			std::cout << "V * D * U" << std::endl;
			print_matrix(V[n-1] * D[n-1] * U[n-1]);
			*/
			//recompute_equal_time_gf(U[n-1], D[n-1], V[n-1], U_r, D_r, V_r);
		}

		void recompute_equal_time_gf(const dmatrix_t& U_l, const dmatrix_t& D_l,
			const dmatrix_t& V_l, const dmatrix_t& U_r, const dmatrix_t& D_r,
			const dmatrix_t& V_r)
		{
			//return;
			svd_solver.compute(U_r.adjoint() * U_l.adjoint() + D_r * (V_r * V_l)
				* D_l, Eigen::ComputeThinU | Eigen::ComputeThinV);
			dmatrix_t D = svd_solver.singularValues().unaryExpr([](double s)
				{ return 1. / s; }).asDiagonal();
			dmatrix_t old_gf = equal_time_gf;
			equal_time_gf = (U_l.adjoint() * svd_solver.matrixV()) * D
				* (svd_solver.matrixU().adjoint() * U_r.adjoint());

//			std::cout << std::endl;
//			std::cout << "current vertex: " << current_vertex << std::endl;
			//std::cout << "recompute_equal_time_gf():" << std::endl;
			//print_matrix(equal_time_gf);
			//std::cout << "correct:" << std::endl;
			dmatrix_t g = (id_N + get_R() * get_L()).inverse();
			//print_matrix(g);
			double diff = (equal_time_gf - g).norm();
			if (diff > 0.001)
			{
				std::cout << "recomputed equal time gf:" << std::endl;
				print_matrix(equal_time_gf);
				std::cout << "id - getR * getL" << std::endl;
				print_matrix(g);
			}
			std::cout << "diff: " << (equal_time_gf - g).norm() << std::endl;
			std::cout << std::endl;
		}

		double measure_M2()
		{
			double M2 = 0.;
			// Use translational symmetry and half-filling here
			//for (int i = 0; i < l.n_sites(); ++i)
			int i = rng() * l.n_sites();
			for (int j = 0; j < l.n_sites(); ++j)
				M2 += equal_time_gf(i, j) * equal_time_gf(i, j)
					/ std::pow(l.n_sites(), 1.0);
			return M2;
		}
		
		void measure_density_correlations(std::vector<double>& corr)
		{
			// Use translational symmetry and half-filling here
			int i = rng() * l.n_sites();
			for (int j = 0; j < l.n_sites(); ++j)
				corr[l.distance(i, j)] += l.parity(i) * l.parity(j)
					* equal_time_gf(i, j) * equal_time_gf(i, j);
		}
	private:
		void print_matrix(const dmatrix_t& m)
		{
			Eigen::IOFormat clean(4, 0, ", ", "\n", "[", "]");
			std::cout << m.format(clean) << std::endl << std::endl;
		}
	private:
		Random& rng;
		const lattice& l;
		const parameters& param;
		std::vector<std::pair<int, int>> lattice_bonds;
		std::vector<int> bond_list;
		int current_vertex;
		int bond_buffer;
		int n_max_order;
		int n_non_ident;
		int n_intervals;
		dmatrix_t equal_time_gf;
		dmatrix_t A; //exp(Lambda_b)
		dmatrix_t invA; //exp(-Lambda_b)
		dmatrix_t B; //exp(Lambda_b) - I
		dmatrix_t invB; //(exp(Lambda_b) - I)^-1
		dmatrix_t C; //exp(Lambda_b) - I
		dmatrix_t invC; //(exp(-Lambda_b) - I)^-1
		dmatrix_t id_2;
		dmatrix_t id_N;
		std::vector<dmatrix_t> U;
		std::vector<dmatrix_t> D;
		std::vector<dmatrix_t> V;
		Eigen::JacobiSVD<dmatrix_t> svd_solver;
};
