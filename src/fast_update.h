#pragma once
#include <vector>
#include <array>
#include <iostream>
#include <Eigen/Dense>
#include "dump.h"
#include "lattice.h"
#include "parameters.h"

template<typename function_t, typename arg_t>
class fast_update
{
	public:
		template<int n, int m>
		using matrix_t = Eigen::Matrix<double, n, m, Eigen::ColMajor>; 
		using dmatrix_t = matrix_t<Eigen::Dynamic, Eigen::Dynamic>;

		fast_update(const function_t& entry_function_, const lattice& l_,
			const parameters& param_)
			: entry_function(entry_function_), l(l_), param(param_)
		{
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
		
		template<int N>
		double try_insert_bond(arg_t& arg)
		{
			matrix_t<2, 2> a = id_2;
			a(0, 0) -= greens_function(arg.bond.first, arg.bond.first);
			a(0, 1) -= greens_function(arg.bond.first, arg.bond.second);
			a(1, 0) -= greens_function(arg.bond.second, arg.bond.first);
			a(1, 1) -= greens_function(arg.bond.second, arg.bond.second);
			matrix_t<2, 2> b = id_2 + bond_exp * a;
			return b.determinant() * l.n_bonds() * param.beta * param.t
				/ ((n_max_order - n_non_ident) * std::sinh(param.lambda));
		}

		template<int N>
		void finish_insert_bond()
		{
		}

		template<int N>
		double try_remove_bond(arg_t& arg)
		{
			matrix_t<2, 2> a = id_2;
			a(0, 0) -= greens_function(arg.bond.first, arg.bond.first);
			a(0, 1) -= greens_function(arg.bond.first, arg.bond.second);
			a(1, 0) -= greens_function(arg.bond.second, arg.bond.first);
			a(1, 1) -= greens_function(arg.bond.second, arg.bond.second);
			matrix_t<2, 2> b = id_2 + inv_bond_exp * a;
			return b.determinant() * ((n_max_order - n_non_ident + 1.) 
				* std::sinh(param.lambda)) / (l.n_bonds() * param.beta * param.t);
		}

		template<int N>
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
		function_t entry_function;
		const lattice& l;
		const parameters& param;
		int n_max_order;
		int n_non_ident;
		dmatrix_t greens_function;
		matrix_t<2, 2> bond_exp;
		matrix_t<2, 2> inv_bond_exp;
		matrix_t<2, 2> id_2;
		dmatrix_t id_N;
};
