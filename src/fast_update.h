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
		using matrix_t = Eigen::Matrix<double, n, m>; 
		using dmatrix_t = matrix_t<Eigen::Dynamic, Eigen::Dynamic>;
		using sparse_t = Eigen::SparseMatrix<double>;

		fast_update(Random& rng_, const lattice& l_, const parameters& param_,
			measurements& measure_)
			: rng(rng_), l(l_), param(param_), measure(measure_),
				n_non_ident{0, 0},
				update_time_displaced_gf(false),
				equal_time_gf{}, time_displaced_gf{},
				stabilizer{measure, equal_time_gf, time_displaced_gf}
		{}

		void initialize()
		{
			for (int bond_type = 0; bond_type < 2; ++bond_type)
				for (int i = 0; i < l.n_sites(); ++i)
					for (auto j : l.neighbors(i, "nearest neighbors"))
						if (i < j)
							lattice_bonds.push_back({i, j});
			id_2 = matrix_t<2, 2>::Identity();
			id_N = dmatrix_t::Identity(l.n_sites(), l.n_sites());
			equal_time_gf = 0.5 * id_N; 
			time_displaced_gf = 0.5 * id_N; 
			
			A.resize(4, dmatrix_t(2, 2));
			A[0] << std::cosh(param.lambda), std::sinh(param.lambda),
				std::sinh(param.lambda), std::cosh(param.lambda);
			A[1] << std::cosh(param.lambda), -std::sinh(param.lambda),
				-std::sinh(param.lambda), std::cosh(param.lambda);
			A[2] << -1., 0., 0., -1;
			A[3] = A[2];

			B.resize(4, dmatrix_t(2, 2));
			for (int i = 0; i < B.size()/2; i+=2)
			{
				B[2*i] = A[2*i] - id_2;
				B[2*i+1] = B[2*i].inverse();
			}
			
			C.resize(4, dmatrix_t(2, 2));
			for (int i = 0; i < C.size()/2; i+=2)
			{
				C[2*i] = A[2*i+1] - id_2;
				C[2*i+1] = C[2*i].inverse();
			}
			
			build_vertex_matrices();
			max_order(5000);
		}

		void build_vertex_matrices()
		{
			vertex_matrices.resize(3. * l.n_bonds() + 1);
			vertex_matrices[0] = sparse_t(l.n_sites(), l.n_sites());
			vertex_matrices[0].setIdentity();
			for (int type = 0; type < 3; ++type)
				for (int i = 1; i <= l.n_bonds(); ++i)
					vertex_matrices[type*l.n_bonds() + i] =
						create_sparse_vertex_matrix(type, i);
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
				<< n_max_order << ", current order is: n_1 = " << n_non_ident[0]
				<< ", n_2 = " << n_non_ident[1] << "." << std::endl;
			std::cout << "bond list size " << bond_list.size() << std::endl;
		}

		int max_order() const
		{
			return n_max_order;
		}

		int non_ident(int bond_type) const
		{
			return n_non_ident[bond_type];
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

		sparse_t create_sparse_vertex_matrix(int type, int bond_id)
		{
			using triplet_t = Eigen::Triplet<double>;
			std::vector<triplet_t> triplet_list;
			triplet_list.reserve(l.n_sites() + 2);
			for (int i = 0; i < l.n_sites(); ++i)
				triplet_list.push_back({i, i, 1.});
			sparse_t lam(l.n_sites(), l.n_sites());
			--bond_id;
			if (bond_id < 0)
			{
				lam.setFromTriplets(triplet_list.begin(), triplet_list.end());
				return lam;
			}
			auto& bond = lattice_bonds[bond_id];
			triplet_list.push_back({bond.first, bond.second, A[type](0, 1)});
			triplet_list.push_back({bond.second, bond.first, A[type](1, 0)});
			for (auto& t : triplet_list)
				if (t.row() == bond.first && t.col() == bond.first)
					t = {bond.first, bond.first, A[type](0, 0)};
				else if(t.row() == bond.second && t.col() == bond.second)
					t = {bond.second, bond.second, A[type](1, 1)};
			lam.setFromTriplets(triplet_list.begin(), triplet_list.end());
			return lam;
		}
		
		sparse_t create_sparse_from_block(const std::pair<int, int>& bond,
			const dmatrix_t& m)
		{
			using triplet_t = Eigen::Triplet<double>;
			std::vector<triplet_t> triplet_list = {{bond.first, bond.first,
				m(0, 0)}, {bond.first, bond.second, m(0, 1)}, {bond.second,
				bond.first, m(1, 0)}, {bond.second, bond.second, m(1, 1)}};
			sparse_t s(l.n_sites(), l.n_sites());
			s.setFromTriplets(triplet_list.begin(), triplet_list.end());
			return s;
		}

		const sparse_t& vertex_matrix(int type, int vertex_id)
		{
			if (type == 3)
				--type;
			if (bond_list[vertex_id-1] == 0)
				return vertex_matrices[0];
			else
				return vertex_matrices[type * l.n_bonds() + bond_list[vertex_id-1]];
		}

		matrix_t<2, 2> vertex_block(const dmatrix_t& m, int vertex_id)
		{
			int bond_id = bond_list[vertex_id-1] - 1;
			if (bond_id < 0)
				return id_2; 
			auto& bond = lattice_bonds[bond_id];
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

		dmatrix_t propagator(int n, int m)
		{
			dmatrix_t P = id_N; 
			for (int i = n; i > m; --i)
			{
				if (bond_list[i-1] == 0) continue;
				P *= vertex_matrix(0, i);
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
				int bond_id = (rng() + bond_type) * l.n_bonds();
				bond_buffer = bond_id + 1;
				auto& bond = lattice_bonds[bond_id];
				matrix_t<2, 2> d = id_2 + B[2*bond_type] * (id_2 - vertex_block(
					equal_time_gf, bond));
				// Insert V1 vertex
				if (bond_type == 0)
					return d.determinant() * l.n_bonds() * param.beta * param.t
						/ ((n_max_order - n_non_ident[0] - n_non_ident[1])
						* std::sinh(param.lambda));
				// Insert V2 vertex
				else
					return d.determinant() * l.n_bonds() * param.beta * param.V2/4.
						/ ((n_max_order - n_non_ident[0] - n_non_ident[1]));
			}
			// Remove bond at vertex
			else if ((bond_type == 0 && get_current_bond() <= l.n_bonds()) ||
				(bond_type == 1 && get_current_bond() > l.n_bonds()))
			{
				bond_buffer = 0;
				matrix_t<2, 2> d = id_2 + C[2*bond_type] * (id_2 - vertex_block(
					equal_time_gf, current_vertex));
				// Insert V1 vertex
				if (bond_type == 0)
					return d.determinant() * ((n_max_order - n_non_ident[0]
						- n_non_ident[1] + 1.) * std::sinh(param.lambda))
						/ (l.n_bonds() * param.beta * param.t);
				// Insert V2 vertex
				else
					return d.determinant() * ((n_max_order - n_non_ident[0]
						- n_non_ident[1] + 1.)) / (l.n_bonds() * param.beta
						* param.V2 / 4.);
			}
			else
				return 0.;
		}

		void finish_update_vertex(int bond_type)
		{
			// Insert bond at vertex
			if (get_current_bond() == 0)
			{
				int bond_id = bond_buffer - 1;
				bond_list[current_vertex-1] = bond_buffer;
				auto& bond = lattice_bonds[bond_id];
				matrix_t<2, 2> d = (B[2*bond_type+1] + (id_2 - vertex_block(
					equal_time_gf, bond))).inverse();
				equal_time_gf = equal_time_gf - (equal_time_gf
					* create_sparse_from_block(bond, d)	* (id_N - equal_time_gf));
				++n_non_ident[bond_type];
			}
			// Remove bond at vertex
			else
			{
				int bond_id = bond_list[current_vertex-1] - 1;
				bond_list[current_vertex-1] = 0;
				auto& bond = lattice_bonds[bond_id];
				matrix_t<2, 2> d = (C[2*bond_type+1] + (id_2 - vertex_block(
					equal_time_gf, bond))).inverse();
				equal_time_gf = equal_time_gf - (equal_time_gf
					* create_sparse_from_block(bond, d)  * (id_N - equal_time_gf));
				--n_non_ident[bond_type];
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
					time_displaced_gf = time_displaced_gf * vertex_matrix(0,
						current_vertex + 1);
			}
			else
			{
				// Wrap equal time gf forwards
				if (bond_list[current_vertex] > 0)
					equal_time_gf = vertex_matrix(0, current_vertex + 1)
						* equal_time_gf * vertex_matrix(1, current_vertex + 1);
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
					time_displaced_gf = vertex_matrix(1, current_vertex)
						* time_displaced_gf;
			}
			else
			{
				// Wrap equal time gf backwards
				if (bond_list[current_vertex - 1] > 0)
					equal_time_gf = vertex_matrix(1, current_vertex)
						* equal_time_gf * vertex_matrix(0, current_vertex);
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
		std::vector<int> n_non_ident;
		int n_intervals;
		bool update_time_displaced_gf;
		dmatrix_t equal_time_gf;
		dmatrix_t time_displaced_gf;
		std::vector<sparse_t> vertex_matrices;
		std::vector<sparse_t> inv_vertex_matrices;
		std::vector<dmatrix_t> A; //exp(Lambda_b), exp(-Lambda_b)
		std::vector<dmatrix_t> B; //exp(Lambda_b) - I, (exp(Lambda_b) - I)^-1
		std::vector<dmatrix_t> C; //exp(-Lambda_b) - I, (exp(-Lambda_b) - I)^-1
		dmatrix_t id_2;
		dmatrix_t id_N;
		stabilizer_t stabilizer;
};
