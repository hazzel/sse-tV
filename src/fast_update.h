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
#include "wick_base.h"

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
		
		void serialize(odump& out)
		{
			int size = bond_list.size();
			out.write(size);
			for (auto b : bond_list)
				out.write(b);
			out.write(current_vertex);
			out.write(n_max_order);
			size = n_non_ident.size();
			out.write(size);
			for (auto b : n_non_ident)
				out.write(b);
		}

		void serialize(idump& in)
		{
			int size = 0; in.read(size);
			bond_list.resize(size);
			for (int i = 0; i < size; ++i)
			{
				int b = 0; in.read(b);
				bond_list[i] = b;
			}
			int n = 0; in.read(n); current_vertex = n;
			in.read(n); int max_order_ = n;
			in.read(size);
			n_non_ident.resize(size);
			for (int i = 0; i < size; ++i)
			{
				int b = 0; in.read(b);
				n_non_ident[i] = b;
			}
			rebuild();
		}

		void initialize(int max_order_)
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
			for (int i = 0; i < B.size()/2; ++i)
			{
				B[2*i] = A[2*i] - id_2;
				B[2*i+1] = B[2*i].inverse();
			}
			
			C.resize(4, dmatrix_t(2, 2));
			for (int i = 0; i < C.size()/2; ++i)
			{
				C[2*i] = A[2*i+1] - id_2;
				C[2*i+1] = C[2*i].inverse();
			}

			build_vertex_matrices();
			max_order(max_order_);
		}

		void build_vertex_matrices()
		{
			vertex_matrices.resize(2. * l.n_bonds() + 1);
			vertex_matrices[0] = sparse_t(l.n_sites(), l.n_sites());
			vertex_matrices[0].setIdentity();
			inv_vertex_matrices.resize(2. * l.n_bonds() + 1);
			inv_vertex_matrices[0] = sparse_t(l.n_sites(), l.n_sites());
			inv_vertex_matrices[0].setIdentity();
			for (int type = 0; type < 2; ++type)
				for (int i = 1; i <= l.n_bonds(); ++i)
				{
					vertex_matrices[type*l.n_bonds() + i] =
						create_sparse_vertex_matrix(2*type, type*l.n_bonds() + i);
					inv_vertex_matrices[type*l.n_bonds() + i] =
						create_sparse_vertex_matrix(2*type+1, type*l.n_bonds() + i);
				}
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
			for (auto& i : n_non_ident)
				i = 0;
			for (auto i : bond_list)
				if (i > 0 && i <= l.n_bonds())
					++n_non_ident[0];
				else if(i > l.n_bonds())
					++n_non_ident[1];
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

		void enable_time_displaced_gf(int direction)
		{
			update_time_displaced_gf = true;
			stabilizer.enable_time_displaced_gf(direction);
		}

		void disable_time_displaced_gf()
		{
			update_time_displaced_gf = false;
			stabilizer.disable_time_displaced_gf();
		}

		void rebuild()
		{
			for (int n = 1; n <= n_intervals; ++n)
			{
				dmatrix_t b = propagator(n * param.n_delta, (n-1) * param.n_delta);
				stabilizer.set(n, b);
			}
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

		const sparse_t& vertex_matrix(int vertex_id)
		{
			if (bond_list[vertex_id-1] == 0)
				return vertex_matrices[0];
			else
				return vertex_matrices[bond_list[vertex_id-1]];
		}
		
		const sparse_t& inv_vertex_matrix(int vertex_id)
		{
			if (bond_list[vertex_id-1] == 0)
				return inv_vertex_matrices[0];
			else
				return inv_vertex_matrices[bond_list[vertex_id-1]];
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
				P *= vertex_matrix(i);
			}
			return P;
		}

		dmatrix_t gf_from_scratch()
		{
			return (id_N + propagator(current_vertex, 0) * propagator(n_max_order,
				current_vertex)).inverse();
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
				bond_buffer = (rng() + bond_type) * l.n_bonds() + 1;
				auto& bond = lattice_bonds[bond_buffer - 1];
				matrix_t<2, 2> d = id_2 + B[2*bond_type] * (id_2 - vertex_block(
					equal_time_gf, bond));
				// Insert V1 vertex
				if (bond_type == 0)
					return d.determinant() * l.n_bonds() * param.beta * param.t
						/ ((n_max_order - n_non_ident[0] - n_non_ident[1])
						* std::sinh(param.lambda));
				// Insert V2 vertex
				else
					return -d.determinant() * l.n_bonds() * param.beta * param.V2/4.
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
					return -d.determinant() * ((n_max_order - n_non_ident[0]
						- n_non_ident[1] + 1.)) / (l.n_bonds() * param.beta
						* param.V2 / 4.);
			}
			else
				return 0.;
		}

		void finish_update_vertex(int bond_type)
		{
			std::pair<int, int> bond;
			matrix_t<2, 2> denom;
			// Insert bond at vertex
			if (get_current_bond() == 0)
			{
				bond = lattice_bonds[bond_buffer - 1];
				bond_list[current_vertex-1] = bond_buffer;
				denom = (B[2*bond_type+1] + (id_2 - vertex_block(equal_time_gf,
					bond))).inverse();
				++n_non_ident[bond_type];
			}
			// Remove bond at vertex
			else
			{
				bond = lattice_bonds[bond_list[current_vertex-1] - 1];
				bond_list[current_vertex-1] = 0;
				denom = (C[2*bond_type+1] + (id_2 - vertex_block(equal_time_gf,
					bond))).inverse();
				--n_non_ident[bond_type];
			}
			dmatrix_t GP(l.n_sites(), 2);
			GP.col(0) = equal_time_gf.col(bond.first);
			GP.col(1) = equal_time_gf.col(bond.second);
			dmatrix_t PIG(2, l.n_sites());
			PIG.row(0) = -equal_time_gf.row(bond.first);
			PIG(0, bond.first) += 1.;
			PIG.row(1) = -equal_time_gf.row(bond.second);
			PIG(1, bond.second) += 1.;
			equal_time_gf.noalias() -= GP * denom * PIG;
		}

		void advance_forward()
		{
			if (current_vertex == n_max_order)
				return;
			if (update_time_displaced_gf)
			{
				// Wrap time displaced gf forwards
				if (bond_list[current_vertex] > 0)
					time_displaced_gf = vertex_matrix(current_vertex + 1)
						* time_displaced_gf;
			}
			// Wrap equal time gf forwards
			if (bond_list[current_vertex] > 0)
				equal_time_gf = vertex_matrix(current_vertex + 1)
					* equal_time_gf * inv_vertex_matrix(current_vertex + 1);
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
					time_displaced_gf = time_displaced_gf
						* vertex_matrix(current_vertex);
			}
			// Wrap equal time gf backwards
			if (bond_list[current_vertex - 1] > 0)
				equal_time_gf = inv_vertex_matrix(current_vertex)
					* equal_time_gf * vertex_matrix(current_vertex);
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
			return measure_M2(equal_time_gf);
		}
		
		double measure_M2(const dmatrix_t& gf)
		{
			double M2 = 0.;
			// Use translational symmetry and half-filling here
			int i = rng() * l.n_sites();
			for (int j = 0; j < l.n_sites(); ++j)
				M2 += gf(i, j) * gf(i, j) / l.n_sites();
			return M2;
		}
		
		double measure_epsilon()
		{
			measure_epsilon(equal_time_gf);
		}
		
		double measure_epsilon(const dmatrix_t& gf)
		{
			double ep = 0.;
			int i = rng() * l.n_sites();
			for (int j : l.neighbors(i, "nearest neighbors"))	
				for (int m = 0; m < l.n_sites(); ++m)
					for (int n : l.neighbors(m, "nearest neighbors"))
					{
						ep += (equal_time_gf(j, i) * equal_time_gf(m, n)
							- time_displaced_gf(m, i) * time_displaced_gf(j, n))
							/ l.n_bonds() * 2./3.;
					}
			return ep;
		}
	
		void measure_density_correlations(std::vector<double>& corr)
		{
			std::fill(corr.begin(), corr.end(), 0.0);
			// Use translational symmetry and half-filling here
			int i = rng() * l.n_sites();
			for (int j = 0; j < l.n_sites(); ++j)
				corr[l.distance(i, j)] += l.parity(i) * l.parity(j)
					* equal_time_gf(i, j) * equal_time_gf(i, j);
		}

		void measure_dynamical_observable(int omega_n_max,
			const std::vector<double>& time_grid,
			std::vector<std::vector<double>>& dyn_mat,
			std::vector<std::vector<double>>& dyn_tau,
			const std::vector<wick_base<dmatrix_t>>& obs)
		{
			// 1 = forward, -1 = backward
			int direction = current_vertex == 0 ? 1 : -1;
			// Time grid for Matsubara frequency measurement
			std::vector<std::vector<double>> random_times(omega_n_max);
			for (auto& vec : random_times)
			{
				vec.resize(n_non_ident[0] + n_non_ident[1] + 1);
				std::for_each(vec.begin(), vec.end(), [&](double& t)
					{ t = rng() * param.beta; } );
				std::sort(vec.begin(), vec.end());
				vec.front() = 0.;
				vec.back() = param.beta;
			}
			// Time grid for imaginary time measurement
			std::vector<int> time_pos(time_grid.size(), 0);
			assign_random_times(time_grid, time_pos);

			std::vector<dmatrix_t> td_gf(time_grid.size());
			std::vector<dmatrix_t> et_gf(time_grid.size());

			enable_time_displaced_gf(direction);
			time_displaced_gf = equal_time_gf;
			int tau_pt = 0, pos_pt = 0, t = 0;
			for (int n = 0; n <= n_max_order; ++n)
			{
				if (current_vertex > 0 && bond_list[current_vertex - 1] > 0)
				{
					// Matsubara frequency measurement
					for (int i = 0; i < dyn_mat.size(); ++i)
					{
						// omega_n = 0 is a special case
						if (omega_n_max > 0)
							dyn_mat[i][0] += obs[i].get_obs(equal_time_gf,
								equal_time_gf, time_displaced_gf, time_displaced_gf)
								* (random_times[0][t+1] - random_times[0][t]);
						for (int omega_n = 1; omega_n < omega_n_max; ++omega_n)
						{
							double omega = 2.*4.*std::atan(1.) * omega_n / param.beta;
							dyn_mat[i][omega_n] += obs[i].get_obs(equal_time_gf,
								equal_time_gf, time_displaced_gf, time_displaced_gf)
								* (std::sin(omega * random_times[omega_n][t+1])
								- std::sin(omega * random_times[omega_n][t])) / omega;
						}
					}
					++t;
				}
				if (current_vertex == 0 || (current_vertex > 0 &&
					bond_list[current_vertex - 1] > 0))
				{
					// Imaginary time measurement
					while (pos_pt < time_pos.size() && tau_pt == time_pos[pos_pt])
					{
//						for (int i = 0; i < dyn_tau.size(); ++i)
//							dyn_tau[i][pos_pt] = obs[i].get_obs(equal_time_gf,
//								time_displaced_gf);
						et_gf[pos_pt] = equal_time_gf;
						td_gf[pos_pt] = time_displaced_gf;
						++pos_pt;
					}
					++tau_pt;
				}
				if (direction == 1 && current_vertex < n_max_order)
				{
					advance_forward();
					stabilize_forward();
				}
				else if (direction == -1 && current_vertex > 0)
				{
					advance_backward();
					stabilize_backward();
				}
			}
			disable_time_displaced_gf();
			if (direction == 1)
				current_vertex = 0;
			else if (direction == -1)
				current_vertex = n_max_order;

			for (int t = 0; t < time_grid.size(); ++t)
				for (int i = 0; i < dyn_tau.size(); ++i)
					dyn_tau[i][t] = obs[i].get_obs(et_gf[0], et_gf[t], td_gf[t],
						td_gf[time_grid.size()-1-t]);
		}
		
		const dmatrix_t& measure_time_displaced_gf()
		{
			return time_displaced_gf;
		}
	private:
		void assign_random_times(const std::vector<double>& time_grid,
			std::vector<int>& time_pos)
		{
			std::vector<double> random_times(n_non_ident[0] + n_non_ident[1]);
			std::for_each(random_times.begin(), random_times.end(), [&](double& t)
				{ t = rng() * param.beta; } );
			std::sort(random_times.begin(), random_times.end());
			for (int i = 0; i < time_pos.size(); ++i)
				time_pos[i] = std::lower_bound(random_times.begin(),
					random_times.end(), time_grid[i]) - random_times.begin();
		}

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
