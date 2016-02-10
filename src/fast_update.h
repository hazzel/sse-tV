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
			n_max_order = 10;
			n_non_ident = 0;
			current_vertex = 0;
			bond_list.resize(n_max_order, 0);
			greens_function = 0.5 * dmatrix_t::Identity(l.n_sites(), l.n_sites());
			id_2 = matrix_t<2, 2>::Identity();
			id_N = dmatrix_t::Identity(l.n_sites(), l.n_sites());
			dmatrix_t V = vertex_matrix(param.lambda);
			dmatrix_t invV = vertex_matrix(-param.lambda);
			A = vertex_block(V);
			invA = vertex_block(invV);
			B = A - id_2;
			//invB = vertex_block((V - id_N).inverse());
			invB = B.inverse();
			std::cout << "B" << std::endl;
			print_matrix(B);
			std::cout << "invB" << std::endl;
			print_matrix(invB);
			C = invA - id_2;
			//invC = vertex_block((invV - id_N).inverse());
			invC = C.inverse();
			std::cout << "C" << std::endl;
			print_matrix(C);
			std::cout << "invC" << std::endl;
			print_matrix(invC);
			
			dmatrix_t lam(l.n_sites(), l.n_sites());
			auto bond = lattice_bonds[0];
			lam(bond.first, bond.second) = param.lambda;
			lam(bond.second, bond.first) = param.lambda;
			Eigen::SelfAdjointEigenSolver<dmatrix_t> solver(lam);
			dmatrix_t d = solver.eigenvalues().unaryExpr([](double e)
				{ return std::exp(e); }).asDiagonal();
			lam = solver.eigenvectors() * d * solver.eigenvectors().adjoint()-id_N;
			lam = lam.inverse().eval();
			std::cout << "(exp(lambda) - 1)^-1" << std::endl;
			print_matrix(lam);
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
			return bond_list[current_vertex];
		}

		void rebuild()
		{
		}

		void serialize(odump& out)
		{
		}

		void serialize(idump& in)
		{
		}

		dmatrix_t vertex_matrix(double lambda, int vertex_id)
		{
			dmatrix_t lam = dmatrix_t::Zero(l.n_sites(), l.n_sites());
			int bond_id = bond_list[vertex_id] - 1;
			std::pair<int, int> bond = lattice_bonds[bond_id];
			if (bond_id < 0)
				return dmatrix_t::Identity(l.n_sites(), l.n_sites());
			lam(bond.first, bond.second) = lambda;
			lam(bond.second, bond.first) = lambda;
			Eigen::SelfAdjointEigenSolver<dmatrix_t> solver(lam);
			dmatrix_t d = solver.eigenvalues().unaryExpr([](double e)
				{ return std::exp(e); }).asDiagonal();
			return solver.eigenvectors() * d * solver.eigenvectors().adjoint();
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

		matrix_t<2, 2> vertex_block(const dmatrix_t& m)
		{
			std::pair<int, int> bond = lattice_bonds[0];
			matrix_t<2, 2> block(2, 2);
			block << m(bond.first, bond.first),
				m(bond.first, bond.second),
				m(bond.second, bond.first),
				m(bond.second, bond.second);
			return block;
		}

		void print_gf_from_scratch()
		{
			std::cout << "current vertex: " << current_vertex
				<< ", non_ident: " << n_non_ident << std::endl;
			dmatrix_t R = dmatrix_t::Identity(l.n_sites(), l.n_sites());
			for (int n = 0; n <= current_vertex; ++n)
			{
				if (bond_list[n] == 0) continue;
				R *= vertex_matrix(param.lambda, n);
			}
			dmatrix_t L = dmatrix_t::Identity(l.n_sites(), l.n_sites());
			for (int n = current_vertex+1; n < max_order(); ++n)
			{
				if (bond_list[n] == 0) continue;
				L *= vertex_matrix(param.lambda, n);
			}

			dmatrix_t G = (id_N + R * L).inverse();
			std::cout << "|G-Gp| = " << (G-greens_function).norm() << std::endl;
			//std::cout << "G from scratch" << std::endl;
			//print_matrix(G);
			//std::cout << "G from update" << std::endl;
			//print_matrix(greens_function);
			std::cout << "=================" << std::endl;
		}
		
		double try_update_vertex(int bond_type)
		{
			// Insert bond at vertex
			if (get_current_bond() == 0)
			{
				matrix_t<2, 2> g = id_2;
				int bond_id = rng() * l.n_bonds();
				bond_buffer = bond_id + 1;
				std::pair<int, int> bond = lattice_bonds[bond_id];
				g(0, 0) -= greens_function(bond.first, bond.first);
				g(0, 1) -= greens_function(bond.first, bond.second);
				g(1, 0) -= greens_function(bond.second, bond.first);
				g(1, 1) -= greens_function(bond.second, bond.second);
				matrix_t<2, 2> d = id_2 + B * g;
				return d.determinant() * l.n_bonds() * param.beta * param.t
					/ ((n_max_order - n_non_ident) * std::sinh(param.lambda));
			}
			// Remove bond at vertex
			else
			{
				matrix_t<2, 2> g = id_2;
				int bond_id = bond_list[current_vertex] - 1;
				bond_buffer = 0;
				std::pair<int, int> bond = lattice_bonds[bond_id];
				g(0, 0) -= greens_function(bond.first, bond.first);
				g(0, 1) -= greens_function(bond.first, bond.second);
				g(1, 0) -= greens_function(bond.second, bond.first);
				g(1, 1) -= greens_function(bond.second, bond.second);
				matrix_t<2, 2> d = id_2 + C * g;
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
				bond_list[current_vertex] = bond_buffer;
				std::pair<int, int> bond = lattice_bonds[bond_id];
				matrix_t<2, 2> g = id_2;
				g(0, 0) -= greens_function(bond.first, bond.first);
				g(0, 1) -= greens_function(bond.first, bond.second);
				g(1, 0) -= greens_function(bond.second, bond.first);
				g(1, 1) -= greens_function(bond.second, bond.second);
				matrix_t<2, 2> d = (invB + g).inverse();
				dmatrix_t e = dmatrix_t::Zero(l.n_sites(), l.n_sites());
				e(bond.first, bond.first) = d(0, 0);
				e(bond.first, bond.second) = d(0, 1);
				e(bond.second, bond.first) = d(1, 0);
				e(bond.second, bond.second) = d(1, 1);
				greens_function = greens_function - (greens_function * e * (id_N - greens_function));
				++n_non_ident;
			}
			// Remove bond at vertex
			else
			{
				int bond_id = bond_list[current_vertex] - 1;
				bond_list[current_vertex] = 0;
				std::pair<int, int> bond = lattice_bonds[bond_id];
				matrix_t<2, 2> g = id_2;
				g(0, 0) -= greens_function(bond.first, bond.first);
				g(0, 1) -= greens_function(bond.first, bond.second);
				g(1, 0) -= greens_function(bond.second, bond.first);
				g(1, 1) -= greens_function(bond.second, bond.second);
				matrix_t<2, 2> d = (invC + g).inverse();
				dmatrix_t e = dmatrix_t::Zero(l.n_sites(), l.n_sites());
				e(bond.first, bond.first) = d(0, 0);
				e(bond.first, bond.second) = d(0, 1);
				e(bond.second, bond.first) = d(1, 0);
				e(bond.second, bond.second) = d(1, 1);
				std::cout << "finish: matrix e" << std::endl;
				print_matrix(e);
				greens_function = greens_function - (greens_function * e * (id_N - greens_function));
				--n_non_ident;
			}
		}

		void advance_forward()
		{
			if (current_vertex == n_max_order - 1)
				return;
			dmatrix_t e = dmatrix_t::Identity(l.n_sites(), l.n_sites());
			dmatrix_t f = dmatrix_t::Identity(l.n_sites(), l.n_sites());
			int bond_id = bond_list[current_vertex+1] - 1;
			std::pair<int, int> bond = lattice_bonds[bond_id];
			if (bond_id >= 0)
			{
				e(bond.first, bond.first) = A(0, 0);
				e(bond.first, bond.second) = A(0, 1);
				e(bond.second, bond.first) = A(1, 0);
				e(bond.second, bond.second) = A(1, 1);
				f(bond.first, bond.first) = invA(0, 0);
				f(bond.first, bond.second) = invA(0, 1);
				f(bond.second, bond.first) = invA(1, 0);
				f(bond.second, bond.second) = invA(1, 1);
			}
			greens_function = e * greens_function * f;
			++current_vertex;
		}
		
		void advance_backward()
		{
			if (current_vertex == 0)
				return;
			dmatrix_t e = dmatrix_t::Identity(l.n_sites(), l.n_sites());
			dmatrix_t f = dmatrix_t::Identity(l.n_sites(), l.n_sites());
			int bond_id = bond_list[current_vertex-1] - 1;
			std::pair<int, int> bond = lattice_bonds[bond_id];
			if (bond_id >= 0)
			{
				e(bond.first, bond.first) = A(0, 0);
				e(bond.first, bond.second) = A(0, 1);
				e(bond.second, bond.first) = A(1, 0);
				e(bond.second, bond.second) = A(1, 1);
				f(bond.first, bond.first) = invA(0, 0);
				f(bond.first, bond.second) = invA(0, 1);
				f(bond.second, bond.first) = invA(1, 0);
				f(bond.second, bond.second) = invA(1, 1);
			}
			greens_function = f * greens_function * e;
			--current_vertex;
		}

		double measure_M2()
		{
			double M2 = 0.;
			int i = 0;
			for (int j = 0; j < l.n_sites(); ++j)
			{
				double delta_ij = (i == j ? 1. : 0.);
				M2 += l.parity(i) * l.parity(j) * ((greens_function(i, i) - 1.)
					* (greens_function(j, j) - 1.) - greens_function(i, j)
					* (greens_function(j, i) - delta_ij)) / l.n_sites();
			}
			return M2 - 0.25;
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
		dmatrix_t greens_function;
		dmatrix_t A; //exp(Lambda_b)
		dmatrix_t invA; //exp(-Lambda_b)
		dmatrix_t B; //exp(Lambda_b) - I
		dmatrix_t invB; //(exp(Lambda_b) - I)^-1
		dmatrix_t C; //exp(Lambda_b) - I
		dmatrix_t invC; //(exp(-Lambda_b) - I)^-1
		dmatrix_t id_2;
		dmatrix_t id_N;
};
