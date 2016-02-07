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
			n_max_order = 1000;
			n_non_ident = 0;
			current_vertex = 0;
			bond_list.resize(n_max_order, 0);
			greens_function = 0.5 * dmatrix_t::Identity(l.n_sites(), l.n_sites());
			A.resize(2, 2);
			A << -1., std::exp(param.lambda), std::exp(param.lambda), -1.;
			invA = A.inverse();
			B.resize(2, 2);
			B << -1., std::exp(-param.lambda), std::exp(-param.lambda), -1.;
			invB = B.inverse();
			C.resize(2, 2);
			C << 0., std::exp(param.lambda), std::exp(param.lambda), 0.;
			invC = C.inverse();
			id_2 = matrix_t<2, 2>::Identity();
			id_N = dmatrix_t::Identity(l.n_sites(), l.n_sites());
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
		
		double try_update_vertex(int bond_type)
		{
			// Insert bond at vertex
			if (get_current_bond() == 0)
			{
				matrix_t<2, 2> g = id_2;
				int bond_id = rng() * l.n_bonds();
				bond_buffer = bond_id + 1;
				g(0, 0) -= greens_function(lattice_bonds[bond_id].first,
					lattice_bonds[bond_id].first);
				g(0, 1) -= greens_function(lattice_bonds[bond_id].first,
					lattice_bonds[bond_id].second);
				g(1, 0) -= greens_function(lattice_bonds[bond_id].second,
					lattice_bonds[bond_id].first);
				g(1, 1) -= greens_function(lattice_bonds[bond_id].second,
					lattice_bonds[bond_id].second);
				matrix_t<2, 2> d = id_2 + A * g;
				return d.determinant() * l.n_bonds() * param.beta * param.t
					/ ((n_max_order - n_non_ident) * std::sinh(param.lambda));
			}
			// Remove bond at vertex
			else
			{
				matrix_t<2, 2> g = id_2;
				int bond_id = bond_list[current_vertex]-1;
				bond_buffer = 0;
				g(0, 0) -= greens_function(lattice_bonds[bond_id].first,
					lattice_bonds[bond_id].first);
				g(0, 1) -= greens_function(lattice_bonds[bond_id].first,
					lattice_bonds[bond_id].second);
				g(1, 0) -= greens_function(lattice_bonds[bond_id].second,
					lattice_bonds[bond_id].first);
				g(1, 1) -= greens_function(lattice_bonds[bond_id].second,
					lattice_bonds[bond_id].second);
				matrix_t<2, 2> d = id_2 + B * g;
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
				matrix_t<2, 2> g = id_2;
				g(0, 0) -= greens_function(lattice_bonds[bond_id].first,
					lattice_bonds[bond_id].first);
				g(0, 1) -= greens_function(lattice_bonds[bond_id].first,
					lattice_bonds[bond_id].second);
				g(1, 0) -= greens_function(lattice_bonds[bond_id].second,
					lattice_bonds[bond_id].first);
				g(1, 1) -= greens_function(lattice_bonds[bond_id].second,
					lattice_bonds[bond_id].second);
				matrix_t<2, 2> d = (invA + g).inverse();
				dmatrix_t e = dmatrix_t::Zero(l.n_sites(), l.n_sites());
				e(lattice_bonds[bond_id].first, lattice_bonds[bond_id].first)
					= d(0, 0);
				e(lattice_bonds[bond_id].first, lattice_bonds[bond_id].second)
					= d(0, 1);
				e(lattice_bonds[bond_id].second, lattice_bonds[bond_id].first)
					= d(1, 0);
				e(lattice_bonds[bond_id].second, lattice_bonds[bond_id].second)
					= d(1, 1);
				greens_function -= (greens_function * e * (id_N - greens_function))
					.eval();
				++n_non_ident;
			}
			// Remove bond at vertex
			else
			{
				int bond_id = bond_list[current_vertex] - 1;
				bond_list[current_vertex] = 0;
				matrix_t<2, 2> g = id_2;
				g(0, 0) -= greens_function(lattice_bonds[bond_id].first,
					lattice_bonds[bond_id].first);
				g(0, 1) -= greens_function(lattice_bonds[bond_id].first,
					lattice_bonds[bond_id].second);
				g(1, 0) -= greens_function(lattice_bonds[bond_id].second,
					lattice_bonds[bond_id].first);
				g(1, 1) -= greens_function(lattice_bonds[bond_id].second,
					lattice_bonds[bond_id].second);
				matrix_t<2, 2> d = (invB + g).inverse();
				dmatrix_t e = dmatrix_t::Zero(l.n_sites(), l.n_sites());
				e(lattice_bonds[bond_id].first, lattice_bonds[bond_id].first)
					= d(0, 0);
				e(lattice_bonds[bond_id].first, lattice_bonds[bond_id].second)
					= d(0, 1);
				e(lattice_bonds[bond_id].second, lattice_bonds[bond_id].first)
					= d(1, 0);
				e(lattice_bonds[bond_id].second, lattice_bonds[bond_id].second)
					= d(1, 1);
				greens_function -= (greens_function * e * (id_N - greens_function))
					.eval();
				--n_non_ident;
			}
			print_matrix(greens_function);
			std::cout << "-----" << std::endl;
		}

		void advance_forward()
		{
			if (current_vertex == n_max_order - 1)
				return;
			dmatrix_t e = dmatrix_t::Zero(l.n_sites(), l.n_sites());
			dmatrix_t f = dmatrix_t::Zero(l.n_sites(), l.n_sites());
			int bond_id = bond_list[current_vertex+1];
			if (bond_id > 0)
			{
				e(lattice_bonds[bond_id-1].first, lattice_bonds[bond_id-1].first)
					= C(0, 0);
				e(lattice_bonds[bond_id-1].first, lattice_bonds[bond_id-1].second)
					= C(0, 1);
				e(lattice_bonds[bond_id-1].second, lattice_bonds[bond_id-1].first)
					= C(1, 0);
				e(lattice_bonds[bond_id-1].second, lattice_bonds[bond_id-1].second)
					= C(1, 1);
				f(lattice_bonds[bond_id-1].first, lattice_bonds[bond_id-1].first)
					= invC(0, 0);
				f(lattice_bonds[bond_id-1].first, lattice_bonds[bond_id-1].second)
					= invC(0, 1);
				f(lattice_bonds[bond_id-1].second, lattice_bonds[bond_id-1].first)
					= invC(1, 0);
				f(lattice_bonds[bond_id-1].second, lattice_bonds[bond_id-1].second)
					= invC(1, 1);
			}
			greens_function = e * greens_function * f;
			++current_vertex;
		}
		
		void advance_backward()
		{
			if (current_vertex == 0)
				return;
			dmatrix_t e = dmatrix_t::Zero(l.n_sites(), l.n_sites());
			dmatrix_t f = dmatrix_t::Zero(l.n_sites(), l.n_sites());
			int bond_id = bond_list[current_vertex-1];
			if (bond_id > 0)
			{
				e(lattice_bonds[bond_id-1].first, lattice_bonds[bond_id-1].first)
					= C(0, 0);
				e(lattice_bonds[bond_id-1].first, lattice_bonds[bond_id-1].second)
					= C(0, 1);
				e(lattice_bonds[bond_id-1].second, lattice_bonds[bond_id-1].first)
					= C(1, 0);
				e(lattice_bonds[bond_id-1].second, lattice_bonds[bond_id-1].second)
					= C(1, 1);
				f(lattice_bonds[bond_id-1].first, lattice_bonds[bond_id-1].first)
					= invC(0, 0);
				f(lattice_bonds[bond_id-1].first, lattice_bonds[bond_id-1].second)
					= invC(0, 1);
				f(lattice_bonds[bond_id-1].second, lattice_bonds[bond_id-1].first)
					= invC(1, 0);
				f(lattice_bonds[bond_id-1].second, lattice_bonds[bond_id-1].second)
					= invC(1, 1);
			}
			greens_function = f * greens_function * e;
			--current_vertex;
		}
	private:
		void print_matrix(const dmatrix_t& m)
		{
			Eigen::IOFormat clean(4, 0, ", ", "\n", "[", "]");
			std::cout << m.format(clean) << std::endl;
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
		dmatrix_t A; //exp(Lambda_b) - I
		dmatrix_t invA; //(exp(Lambda_b) - I)^-1
		dmatrix_t B; //exp(Lambda_b) - I
		dmatrix_t invB; //(exp(-Lambda_b) - I)^-1
		dmatrix_t C; //exp(Lambda_b)
		dmatrix_t invC; //exp(-Lambda_b)
		dmatrix_t id_2;
		dmatrix_t id_N;
};
