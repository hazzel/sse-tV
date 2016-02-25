#pragma once
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "measurements.h"
#include "dump.h"
#include "lattice.h"
#include "parameters.h"
#include "Random.h"

template<typename stabilizer_t>
class fast_update
{
	public:
		template<int n, int m>
		using matrix_t = Eigen::Matrix<double, n, m, Eigen::ColMajor>; 
		using dmatrix_t = matrix_t<Eigen::Dynamic, Eigen::Dynamic>;
		using sparse_t = Eigen::SparseMatrix<double, Eigen::ColMajor>;

		fast_update(Random& rng_, const lattice& l_, const parameters& param_,
			measurements& measure_)
			: rng(rng_), l(l_), param(param_), measure(measure_),
				update_time_displaced_gf(false),
				equal_time_gf{}, time_displaced_gf{},
				stabilizer{measure, equal_time_gf, time_displaced_gf}
		{}

		void initialize()
		{
			for (int i = 0; i < l.n_sites(); ++i)
				for (auto j : l.neighbors(i, "nearest neighbors"))
					if (i < j)
						lattice_bonds.push_back({i, j});
			n_non_ident = 0;
			id_2 = matrix_t<2, 2>::Identity();
			id_N = dmatrix_t::Identity(l.n_sites(), l.n_sites());
			equal_time_gf = 0.5 * id_N; 
			time_displaced_gf = 0.5 * id_N; 
			dmatrix_t expLam = vertex_matrix(param.lambda);
			dmatrix_t invExpLam = vertex_matrix(-param.lambda);
			A = vertex_block(expLam);
			invA = vertex_block(invExpLam);
			B = A - id_2;
			invB = B.inverse();
			C = invA - id_2;
			invC = C.inverse();
			max_order(500);
		}

		void max_order(int n_max_order_)
		{
			int old_max_order = n_max_order;
			n_max_order = n_max_order_;
			if (n_max_order_ % param.n_delta != 0)
				n_max_order += param.n_delta - n_max_order_ % param.n_delta;
			current_vertex = n_max_order;
			n_intervals = n_max_order / param.n_delta;
			if (old_max_order == n_max_order)
				return;
			else if (old_max_order < n_max_order)
				bond_list.resize(n_max_order, 0);
			else
			{
				bond_list.erase(std::remove(bond_list.begin(), bond_list.end(), 0),
					bond_list.end());
				bond_list.resize(n_max_order, 0);
			}
			stabilizer.resize(n_intervals, l.n_sites());
			rebuild();
			std::cout << "Max order set from " << old_max_order << " to "
				<< n_max_order << ", current order is " << n_non_ident << std::endl;
		}

		int max_order() const
		{
			return n_max_order;
		}

		int non_ident() const
		{
			return n_non_ident;
		}

		int get_current_bond() const
		{
			return bond_list[current_vertex-1];
		}

		int get_current_vertex() const
		{
			return current_vertex;
		}

		void enable_time_displaced_gf()
		{
			update_time_displaced_gf = true;
		}

		void disable_time_displaced_gf()
		{
			update_time_displaced_gf = false;
		}

		void rebuild()
		{
			for (int n = 1; n <= n_intervals; ++n)
			{
				dmatrix_t b = propagator(n * param.n_delta, (n-1) * param.n_delta);
				stabilizer.set(n, b);
			}
		}

		void serialize(odump& out)
		{
		}

		void serialize(idump& in)
		{
		}

		dmatrix_t vertex_matrix_dense(double lambda, int vertex_id)
		{
			dmatrix_t lam = id_N; 
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

		sparse_t vertex_matrix(double lambda, int vertex_id)
		{
			using triplet_t = Eigen::Triplet<double>;
			std::vector<triplet_t> triplet_list;
			triplet_list.reserve(l.n_sites() + 2);
			for (int i = 0; i < l.n_sites(); ++i)
				triplet_list.push_back({i, i, 1.});
			sparse_t lam(l.n_sites(), l.n_sites());
			int bond_id = bond_list[vertex_id-1] - 1;
			if (bond_id < 0)
			{
				lam.setFromTriplets(triplet_list.begin(), triplet_list.end());
				return lam;
			}
			std::pair<int, int> bond = lattice_bonds[bond_id];
			if (lambda > 0.)
			{
				triplet_list.push_back({bond.first, bond.second, A(0, 1)});
				triplet_list.push_back({bond.second, bond.first, A(1, 0)});
				for (auto& t : triplet_list)
					if (t.row() == bond.first && t.col() == bond.first)
						t = {bond.first, bond.first, A(0, 0)};
					else if(t.row() == bond.second && t.col() == bond.second)
						t = {bond.second, bond.second, A(1, 1)};
			}
			else
			{
				triplet_list.push_back({bond.first, bond.second, invA(0, 1)});
				triplet_list.push_back({bond.second, bond.first, invA(1, 0)});
				for (auto& t : triplet_list)
					if (t.row() == bond.first && t.col() == bond.first)
						t = {bond.first, bond.first, invA(0, 0)};
					else if(t.row() == bond.second && t.col() == bond.second)
						t = {bond.second, bond.second, invA(1, 1)};
			}
			lam.setFromTriplets(triplet_list.begin(), triplet_list.end());
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
				return id_2; 
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
			dmatrix_t R = id_N; 
			for (int n = current_vertex; n >= 1; --n)
			{
				if (bond_list[n-1] == 0) continue;
				R *= vertex_matrix(param.lambda, n);
			}
			return R;
		}
		
		dmatrix_t get_L()
		{
			dmatrix_t L = id_N; 
			for (int n = n_max_order; n > current_vertex; --n)
			{
				if (bond_list[n-1] == 0) continue;
				L *= vertex_matrix(param.lambda, n);
			}
			return L;
		}
		
		dmatrix_t propagator(int n, int m)
		{
			dmatrix_t P = id_N; 
			for (int i = n; i > m; --i)
			{
				if (bond_list[i-1] == 0) continue;
				P *= vertex_matrix(param.lambda, i);
			}
			return P;
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
					* std::sinh(param.lambda)) / (l.n_bonds() * param.beta
					* param.t);
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

		void advance_forward()
		{
			if (current_vertex == n_max_order)
				return;
			if (update_time_displaced_gf)
			{
				// Wrap time displaced gf forwards
				if (bond_list[current_vertex] > 0)
					time_displaced_gf = time_displaced_gf * vertex_matrix(
						param.lambda, current_vertex + 1);
			}
			else
			{
				// Wrap equal time gf forwards
				if (bond_list[current_vertex] > 0)
					equal_time_gf = vertex_matrix(param.lambda, current_vertex + 1)
						* equal_time_gf * vertex_matrix(-param.lambda,
						current_vertex + 1);
			}
			++current_vertex;
		}
		
		void advance_backward()
		{
			if (current_vertex == 0)
				return;
			if (update_time_displaced_gf)
			{
				// Wrap time displaced gf forwards
				if (bond_list[current_vertex - 1] > 0)
					time_displaced_gf = vertex_matrix(-param.lambda, current_vertex)
						* time_displaced_gf;
			}
			else
			{
				// Wrap equal time gf backwards
				if (bond_list[current_vertex - 1] > 0)
					equal_time_gf = vertex_matrix(-param.lambda, current_vertex)
						* equal_time_gf * vertex_matrix(param.lambda, current_vertex);
			}
			--current_vertex;
		}
		
		void stabilize_forward()
		{
			if (current_vertex % param.n_delta != 0)
				return;
			// n = 0, ..., n_intervals - 1
			int n = current_vertex / param.n_delta - 1;
			dmatrix_t b = propagator((n+1) * param.n_delta, n * param.n_delta);
			stabilizer.stabilize_forward(n, b);
		}
	
		void stabilize_backward()
		{
			if (current_vertex % param.n_delta != 0)
				return;
			//n = n_intervals, ..., 1 
			int n = current_vertex / param.n_delta + 1;
			dmatrix_t b = propagator(n * param.n_delta, (n-1) * param.n_delta);
			stabilizer.stabilize_backward(n, b);
		}

		double measure_M2()
		{
			double M2 = 0.;
			// Use translational symmetry and half-filling here
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

		const dmatrix_t& measure_time_displaced_gf()
		{
			return time_displaced_gf;
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
		measurements& measure;
		std::vector<std::pair<int, int>> lattice_bonds;
		std::vector<int> bond_list;
		int current_vertex;
		int bond_buffer;
		int n_max_order;
		int n_non_ident;
		int n_intervals;
		bool update_time_displaced_gf;
		dmatrix_t equal_time_gf;
		dmatrix_t time_displaced_gf;
		dmatrix_t A; //exp(Lambda_b)
		dmatrix_t invA; //exp(-Lambda_b)
		dmatrix_t B; //exp(Lambda_b) - I
		dmatrix_t invB; //(exp(Lambda_b) - I)^-1
		dmatrix_t C; //exp(Lambda_b) - I
		dmatrix_t invC; //(exp(-Lambda_b) - I)^-1
		dmatrix_t id_2;
		dmatrix_t id_N;
		stabilizer_t stabilizer;
};
