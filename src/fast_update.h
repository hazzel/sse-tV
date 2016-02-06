#pragma once
#include <vector>
#include <array>
#include <iostream>
#include <Eigen/Dense>
#include "dump.h"
#include "lattice.h"
#include "parameters.h"

class fast_update
{
	public:
		template<int n, int m>
		using matrix_t = Eigen::Matrix<double, n, m, Eigen::ColMajor>; 
		using dmatrix_t = matrix_t<Eigen::Dynamic, Eigen::Dynamic>;

		fast_update(const lattice& l_, const parameters& param_)
			: l(l_), param(param_)
		{
			for (int i = 0; i < l.n_sites(); ++i)
				for (auto j : l.neighbors(i, "nearest neighbors"))
					if (i < j)
						lattice_bonds.push_back({i, j});
			n_max_order = 100000;
			n_non_ident = 0;
			bond_list.resize(n_max_order, 0);
			greens_function = dmatrix_t::Identity(l.n_sites(), l.n_sites());
			bond_exp << 0., std::exp(param.lambda), std::exp(param.lambda), 0.;
			inv_bond_exp << 0., std::exp(-param.lambda), std::exp(-param.lambda), 0.;
			id_2 = matrix_t<2, 2>::Identity();
			id_N = dmatrix_t::Identity(l.n_sites(), l.n_sites());
		}

		int perturbation_order() const
		{
			return 0;
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
		
		double try_insert_bond()
		{
			matrix_t<2, 2> a = id_2;
			int bond_id = bond_list[vertex_id]-1;
			a(0, 0) -= greens_function(lattice_bonds[bond_id].first,
				lattice_bonds[bond_id].first);
			a(0, 1) -= greens_function(lattice_bonds[bond_id].first,
				lattice_bonds[bond_id].second);
			a(1, 0) -= greens_function(lattice_bonds[bond_id].second,
				lattice_bonds[bond_id].first);
			a(1, 1) -= greens_function(lattice_bonds[bond_id].second,
				lattice_bonds[bond_id].second);
			matrix_t<2, 2> b = id_2 + bond_exp * a;
			return b.determinant() * l.n_bonds() * param.beta * param.t
				/ ((n_max_order - n_non_ident) * std::sinh(param.lambda));
		}

		void finish_insert_bond()
		{

		}

		double try_remove_bond()
		{
			matrix_t<2, 2> a = id_2;
			int bond_id = bond_list[vertex_id]-1;
			a(0, 0) -= greens_function(lattice_bonds[bond_id].first,
				lattice_bonds[bond_id].first);
			a(0, 1) -= greens_function(lattice_bonds[bond_id].first,
				lattice_bonds[bond_id].second);
			a(1, 0) -= greens_function(lattice_bonds[bond_id].second,
				lattice_bonds[bond_id].first);
			a(1, 1) -= greens_function(lattice_bonds[bond_id].second,
				lattice_bonds[bond_id].second);
			matrix_t<2, 2> b = id_2 + inv_bond_exp * a;
			return b.determinant() * ((n_max_order - n_non_ident + 1.) 
				* std::sinh(param.lambda)) / (l.n_bonds() * param.beta * param.t);
		}

		void finish_remove_bond()
		{
		}

	private:
		void print_matrix(const dmatrix_t& m)
		{
			Eigen::IOFormat clean(4, 0, ", ", "\n", "[", "]");
			std::cout << m.format(clean) << std::endl;
		}
	private:
		const lattice& l;
		const parameters& param;
		std::vector<std::pair<int, int>> lattice_bonds;
		std::vector<int> bond_list;
		int vertex_id;
		int n_max_order;
		int n_non_ident;
		dmatrix_t greens_function;
		dmatrix_t bond_exp;
		dmatrix_t inv_bond_exp;
		dmatrix_t id_2;
		dmatrix_t id_N;
};
