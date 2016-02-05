#pragma once
#include <vector>
#include <array>
#include <iostream>
#include <Eigen/Dense>
#include "dump.h"
#include "lattice.h"

template<int N, typename T>
struct helper_matrices
{
	template<int n, int m>
	using matrix_t = Eigen::Matrix<T, n, m, Eigen::ColMajor>;
	
	//TODO: fix allignment issues when fixed size matrix_t<N, N>
	matrix_t<Eigen::Dynamic, Eigen::Dynamic> u;
	matrix_t<Eigen::Dynamic, Eigen::Dynamic> v;
	matrix_t<Eigen::Dynamic, Eigen::Dynamic> Mu;
	matrix_t<Eigen::Dynamic, Eigen::Dynamic> a;
	matrix_t<Eigen::Dynamic, Eigen::Dynamic> S;
	matrix_t<Eigen::Dynamic, Eigen::Dynamic> m;

};

template<typename function_t, typename arg_t>
class fast_update
{
	public:
		template<int n, int m>
		using matrix_t = Eigen::Matrix<double, n, m, Eigen::ColMajor>; 
		using dmatrix_t = matrix_t<Eigen::Dynamic, Eigen::Dynamic>;

		fast_update(const function_t& entry_function_, const lattice& l_,
			int n_flavors_)
			: entry_function(entry_function_), l(l_), flavor_cnt(n_flavors_, 0)
		{}

		int perturbation_order(int flavor) const
		{
			return flavor_cnt[flavor] / 2;
		}

		const arg_t& vertex(int index, int flavor)
		{
			int offset = 0;
			for (int f = 0; f < flavor; ++f)
				offset += flavor_cnt[f];
			return vertices[offset + index]; 
		}

		void print_vertices() const
		{
			std::cout << "print vertices:" << std::endl;
			int offset = 0;
			for (int f = 0; f < flavor_cnt.size(); f++)
			{
				std::cout << "flavor " << f << std::endl;
				for (int i = 0; i < flavor_cnt[f]/2; ++i)
				{
					std::cout << 2*i << ": " << vertices[offset+2*i].tau << ", "
						<< vertices[offset+2*i].site << " , w: "
						<< vertices[offset+2*i].worm << " ; "
						<< vertices[offset+2*i+1].tau
						<< ", " << vertices[offset+2*i+1].site << " , w: "
						<< vertices[offset+2*i+1].worm << std::endl;
				}
				offset += flavor_cnt[f];
			}
		}

		void build(std::vector<arg_t>& args, std::vector<int>&& flavors)
		{
			vertices = std::move(args);
			flavor_cnt = std::move(flavors);
			M.resize(vertices.size(), vertices.size());
			rebuild();
		}

		void rebuild()
		{
			if (M.rows() == 0) return;
			for (int i = 0; i < M.rows(); ++i)
			{
				M(i, i) = 0.0;
				for (int j = i+1; j < M.cols(); ++j)
				{
					M(i, j) = entry_function(vertices[i], vertices[j]);
					M(j, i) = -M(i, j) * l.parity(vertices[i].site)
						* l.parity(vertices[j].site);
				}
			}
			M = M.inverse().eval();
		}

		void serialize(odump& out)
		{
			int size = vertices.size();
			out.write(size);
			for (arg_t& v : vertices)
				v.serialize(out);
			size = flavor_cnt.size();
			out.write(size);
			for (int f : flavor_cnt)
				out.write(f);
		}

		void serialize(idump& in)
		{
			vertices.clear();
			int size; in.read(size);
			for (int i = 0; i < size; ++i)
			{
				arg_t v;
				v.serialize(in);
				vertices.push_back(v);
			}
			in.read(size);
			assert(size == flavor_cnt.size() && "serialization error");
			for (int i = 0; i < size; ++i)
			{
				int f; in.read(f);
				flavor_cnt[i] = f;
			}
			M.resize(vertices.size(), vertices.size());
			rebuild();
		}
		
		template<int N>
		std::complex<double> matsubara_gf(int omega_n, double beta,
			std::vector<arg_t>& args, int flavor=0)
		{
			typedef std::complex<double> C;
			int k = M.rows();
			const int n = 2*N;
			last_flavor = flavor;
			
			arg_buffer = std::move(args);
			helper<n, C>().u.resize(k, n);
			helper<n, C>().v.resize(n, k);
			helper<n, C>().a.resize(n, n);
			fill_matsubara_matrices<N>(omega_n);
			helper<n, C>().Mu.noalias() = M * helper<n, C>().u / beta;
			helper<n, C>().S = helper<n, C>().a;
			helper<n, C>().S.noalias() -= helper<n, C>().v * helper<n, C>().Mu;

			return helper<n, C>().S.determinant();
		}

		template<int N>
		double try_add(std::vector<arg_t>& args, int flavor=0)
		{
			int k = M.rows();
			const int n = 2*N;
			last_flavor = flavor;
			
			arg_buffer = std::move(args);
			helper<n>().u.resize(k, n);
			helper<n>().v.resize(n, k);
			helper<n>().a.resize(n, n);
			fill_helper_matrices<N>();
			helper<n>().Mu.noalias() = M * helper<n>().u;
			helper<n>().S = helper<n>().a;
			helper<n>().S.noalias() -= helper<n>().v * helper<n>().Mu;

			return helper<n>().S.determinant();
		}

		template<int N>
		void finish_add()
		{
			int k = M.rows();
			const int n = 2*N;

			helper<n>().S = helper<n>().S.inverse().eval();
			dmatrix_t vM = M.transpose() * helper<n>().v.transpose();
			vM.transposeInPlace();
			M.conservativeResize(k + n, k + n);

			M.block(k, 0, n, k).noalias() = -helper<n>().S * vM;
			M.block(0, 0, k, k).noalias() -= helper<n>().Mu
				* M.block(k, 0, n, k);
			M.block(0, k, k, n).noalias() = -helper<n>().Mu * helper<n>().S;
			M.template block<n, n>(k, k) = helper<n>().S;
			
			vertices.insert(vertices.end(), arg_buffer.begin(), arg_buffer.end());
			pos_buffer.resize(n/2);
			for (int i = 0; i < n/2; ++i)
				pos_buffer[i] = vertices.size() - n + 2*i;
			if(last_flavor == 0)
				permute_backward();
			flavor_cnt[last_flavor] += n;
		}

		template<int N>
		double try_remove(std::vector<int>& pos, int flavor=0)
		{
			const int n = 2*N;
			if (flavor_cnt[flavor] < n)
				return 0.0;
			last_flavor = flavor;
			pos_buffer = std::move(pos);
			int pos_offset = 0;
			for (int f = 0; f < flavor; ++f)
				pos_offset += flavor_cnt[f];
			for (int& p : pos_buffer)
				p = 2*p + pos_offset;
			fill_S_matrix<N>();
			return helper<n>().S.determinant();
		}

		template<int N>
		void finish_remove()
		{
			permute_forward();
			int k = M.rows();
			const int n = 2*N;
			
			helper<n>().S.transposeInPlace();
			dmatrix_t t = M.block(k - n, 0, n, k - n).transpose()
				* helper<n>().S.inverse();
			t.transposeInPlace();
			M.block(0, 0, k - n, k - n).noalias()
				-= M.block(0, k - n, k - n, n) * t;
			M.conservativeResize(k - n, k - n);

			for (int i = 0; i < n; ++i)
				vertices.erase(vertices.end() - 1);
			flavor_cnt[last_flavor] -= n;
		}

		//Assume that vertices to shift are already located at the end
		//TODO: include permute operation
		template<int N>
		double try_shift(std::vector<arg_t>& args)
		{
			int k = M.rows() - args.size();
			const int n = 2*N;
			last_flavor = 1; //shift worm vertices
			
			arg_buffer = std::move(args);
			helper<n>().u.resize(k, n);
			helper<n>().v.resize(n, k);
			helper<n>().a.resize(n, n);
			helper<n>().m.resize(k, k);
			fill_helper_matrices<N>();
			
			helper<n>().m = M.block(0, 0, k, k);
			helper<n>().m.noalias() -= M.block(0, k, k, n)
				* M.template block<n, n>(k, k).inverse() * M.block(k, 0, n, k);
			helper<n>().S = helper<n>().a;
			helper<n>().S.noalias() -= helper<n>().v * helper<n>().m
				* helper<n>().u;
			return helper<n>().S.determinant()
				* M.template block<n, n>(k, k).determinant();
		}
		
		template<int N>
		void finish_shift()
		{
			int k = M.rows() - arg_buffer.size();
			const int n = 2*N;
			
			helper<n>().S = helper<n>().S.inverse().eval();
			dmatrix_t vM = helper<n>().v * helper<n>().m;
			helper<n>().Mu.noalias() = helper<n>().m * helper<n>().u;
			M.block(k, 0, n, k).noalias() = -helper<n>().S * vM;
			M.block(0, 0, k, k) = helper<n>().m;
			M.block(0, 0, k, k).noalias() -= helper<n>().Mu
				* M.block(k, 0, n, k);
			M.block(0, k, k, n).noalias() = -helper<n>().Mu * helper<n>().S;
			M.template block<n, n>(k, k) = helper<n>().S;
			
			for (int i = 0; i < arg_buffer.size(); ++i)
				vertices[vertices.size() - arg_buffer.size() + i] = arg_buffer[i];
		}
	private:
		template<int N, typename T> struct type {};
		template<int N, typename T=double>
		helper_matrices<N, T>& helper()
		{
			return helper(type<N, T>());
		}
		helper_matrices<2, double>& helper(type<2, double>) { return helper_2; }
		helper_matrices<2, std::complex<double>>& helper(type<2,
			std::complex<double>>) { return helper_2_c; }
		helper_matrices<4, double>& helper(type<4, double>) { return helper_4; }

		template<int N>
		void fill_helper_matrices()
		{
			const int n = 2*N;
			for (int i = 0; i < helper<n>().u.rows(); ++i)
			{
				for (int j = 0; j < n; ++j)
				{
					helper<n>().u(i, j) = entry_function(vertices[i],
						arg_buffer[j]);
					helper<n>().v(j, i) = -helper<n>().u(i, j)
						* l.parity(vertices[i].site) * l.parity(arg_buffer[j].site);
				}
			}
			for (int i = 0; i < n; ++i)
			{
				helper<n>().a(i, i) = 0.0;
				for (int j = i+1; j < n; ++j)
				{
					helper<n>().a(i, j) = entry_function(arg_buffer[i],
						arg_buffer[j]);
					helper<n>().a(j, i) = -helper<n>().a(i, j)
						* l.parity(arg_buffer[i].site) * l.parity(arg_buffer[j].site);
				}
			}
		}
	
		template<int N>	
		void fill_matsubara_matrices(int omega_n)
		{
			typedef std::complex<double> C;
			const int n = 2*N;
			for (int i = 0; i < helper<n, C>().u.rows(); ++i)
			{
				for (int j = 0; j < n; ++j)
				{
					helper<n, C>().u(i, j) = entry_function.matsubara_frequency(
						omega_n, vertices[i], arg_buffer[j]);
					helper<n, C>().v(j, i) = -helper<n, C>().u(i, j)
						* l.parity(vertices[i].site) * l.parity(arg_buffer[j].site);
				}
			}
			for (int i = 0; i < n; ++i)
			{
				helper<n, C>().a(i, i) = 0.0;
				for (int j = i+1; j < n; ++j)
				{
					helper<n, C>().a(i, j) = entry_function.matsubara_frequency(
						omega_n, arg_buffer[i], arg_buffer[j]);
					helper<n, C>().a(j, i) = -helper<n, C>().a(i, j)
						* l.parity(arg_buffer[i].site) * l.parity(arg_buffer[j].site);
				}
			}
		}

		template<int N>
		void fill_S_matrix()
		{
			helper<2*N>().S.resize(2*N, 2*N);
			for (int i = 0; i < N; ++i)
			{
				for (int j = 0; j < N; ++j)
				{
					helper<2*N>().S.template block<2, 2>(2*i, 2*j) =
						M.template block<2, 2>(pos_buffer[i], pos_buffer[j]);
				}
			}
		}

		void swap_rows_cols(int i, int j)
		{
			int k = M.rows();
			dmatrix_t cols = M.block(0, i, k, 2);
			M.block(0, i, k, 2) = M.block(0, j, k, 2);
			M.block(0, j, k, 2) = cols;
			dmatrix_t rows = M.block(i, 0, 2, k);
			M.block(i, 0, 2, k) = M.block(j, 0, 2, k);
			M.block(j, 0, 2, k) = rows;
		}
		
		void permute_forward()
		{
			int n = 2*pos_buffer.size();
			int block_end = 0;
			for (int f = 0; f < last_flavor; ++f)
				block_end += flavor_cnt[f];
			for (int f = last_flavor; f < flavor_cnt.size(); ++f)
			{
				block_end += flavor_cnt[f];
				for (int i = 0; i < n/2; ++i)
				{
					swap_rows_cols(pos_buffer[n/2 - i - 1], block_end - 2*i - 2);
					std::swap(vertices[pos_buffer[n/2 - i - 1]],
						vertices[block_end - 2*i - 2]);
					std::swap(vertices[pos_buffer[n/2 - i - 1] + 1],
						vertices[block_end - 2*i - 1]);
					pos_buffer[n/2 - i - 1] = block_end - 2*i - 2;
				}
			}
		}

		void permute_backward()
		{
			int n = 2*pos_buffer.size();
			int block_end = 0;
			for (int f = 0; f < flavor_cnt.size(); ++f)
				block_end += flavor_cnt[f];
			for (int f = flavor_cnt.size()-1; f > last_flavor; --f)
			{
				block_end -= flavor_cnt[f];
				for (int i = 0; i < n/2; ++i)
				{
					swap_rows_cols(pos_buffer[i], block_end + 2*i);
					std::swap(vertices[pos_buffer[i]], vertices[block_end + 2*i]);
					std::swap(vertices[pos_buffer[i] + 1],
						vertices[block_end + 2*i + 1]);
					pos_buffer[i] = block_end + 2*i;
				}
			}
		}
		
		void print_matrix(const dmatrix_t& m)
		{
			Eigen::IOFormat clean(4, 0, ", ", "\n", "[", "]");
			std::cout << m.format(clean) << std::endl;
		}
	private:
		function_t entry_function;
		const lattice& l;
		std::vector<arg_t> vertices;
		std::vector<int> flavor_cnt;
		std::vector<arg_t> arg_buffer;
		std::vector<int> pos_buffer;
		int last_flavor;
		dmatrix_t M;
		helper_matrices<2, double> helper_2;
		helper_matrices<4, double> helper_4;
		helper_matrices<2, std::complex<double>> helper_2_c;
};
