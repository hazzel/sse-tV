#pragma once
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include "measurements.h"
#include "dump.h"
#include "lattice.h"
#include "parameters.h"
#include "Random.h"

class svd_stabilizer
{
	public:
		template<int n, int m>
		using matrix_t = Eigen::Matrix<double, n, m, Eigen::ColMajor>; 
		using dmatrix_t = matrix_t<Eigen::Dynamic, Eigen::Dynamic>;
		
		svd_stabilizer(measurements& measure_, dmatrix_t& equal_time_gf_, dmatrix_t&
			time_displaced_gf_)
			: measure(measure_), update_time_displaced_gf(false), equal_time_gf(equal_time_gf_),
				time_displaced_gf(time_displaced_gf_)
		{}
		
		void resize(int n_intervals_, int dimension)
		{
			id_N = dmatrix_t::Identity(dimension, dimension);
			n_intervals = n_intervals_;
			U.resize(n_intervals + 1);
			D.resize(n_intervals + 1);
			V.resize(n_intervals + 1);
			for (int n = 0; n < n_intervals + 1; ++n)
			{
				U[n] = id_N; D[n] = id_N; V[n] = id_N;
			}
		}
		
		void set(int n, const dmatrix_t& b)
		{
			svd_solver.compute(b * U[n-1] * D[n-1], Eigen::ComputeThinU
				| Eigen::ComputeThinV);
			U[n] = svd_solver.matrixU();
			D[n] = svd_solver.singularValues().asDiagonal();
			V[n] = svd_solver.matrixV().adjoint() * V[n-1];
			if (n == n_intervals)
			{
				recompute_time_displaced_gf(id_N, id_N, id_N, U.back(),
					D.back(), V.back());
				U.back() = id_N; D.back() = id_N; V.back() = id_N;
			}
		}
		
		// n = 0, ..., n_intervals - 1
		void stabilize_forward(int n, const dmatrix_t& b)
		{
			if (n == 0)
			{
				U.front() = id_N; D.front() = id_N; V.front() = id_N;
			}
			
			svd_solver.compute(b * U[n] * D[n], Eigen::ComputeThinU |
				Eigen::ComputeThinV);
			U_l = U[n+1]; D_l = D[n+1]; V_l = V[n+1];
			U[n+1] = svd_solver.matrixU();
			D[n+1] = svd_solver.singularValues().asDiagonal();
			V[n+1] = svd_solver.matrixV().adjoint() * V[n];
			
			if (update_time_displaced_gf)
				recompute_time_displaced_gf(U_l, D_l, V_l, U[n+1], D[n+1],
					V[n+1]);
			else
				recompute_equal_time_gf(U_l, D_l, V_l, U[n+1], D[n+1], V[n+1]);
		}
	
		//n = n_intervals, ..., 1 
		void stabilize_backward(int n, const dmatrix_t& b)
		{
			if (n == n_intervals)
			{
				U.back() = id_N; D.back() = id_N; V.back() = id_N;
			}
			
			svd_solver.compute(D[n] * U[n] * b, Eigen::ComputeThinU |
				Eigen::ComputeThinV);
			U_r = U[n-1]; D_r = D[n-1]; V_r = V[n-1];
			V[n-1] = V[n] * svd_solver.matrixU();
			D[n-1] = svd_solver.singularValues().asDiagonal();
			U[n-1] = svd_solver.matrixV().adjoint();
		
			if (update_time_displaced_gf)
				recompute_time_displaced_gf(U[n-1], D[n-1], V[n-1], U_r, D_r,
					V_r);
			else
				recompute_equal_time_gf(U[n-1], D[n-1], V[n-1], U_r, D_r, V_r);
		}
		
		void recompute_equal_time_gf(const dmatrix_t& U_l_, const dmatrix_t& D_l_,
			const dmatrix_t& V_l_, const dmatrix_t& U_r_, const dmatrix_t& D_r_,
			const dmatrix_t& V_r_)
		{
			dmatrix_t old_gf = equal_time_gf;
			svd_solver.compute(U_r_.adjoint() * U_l_.adjoint() + D_r_ * (V_r_ * V_l_)
				* D_l_, Eigen::ComputeThinU | Eigen::ComputeThinV);
			dmatrix_t D = svd_solver.singularValues().unaryExpr([](double s)
				{ return 1. / s; }).asDiagonal();
			equal_time_gf = (U_l_.adjoint() * svd_solver.matrixV()) * D
				* (svd_solver.matrixU().adjoint() * U_r_.adjoint());
			measure.add("norm error", (old_gf - equal_time_gf).norm());
			measure.add("max error", (old_gf - equal_time_gf).lpNorm<Eigen::
				Infinity>());
			measure.add("avg error", (old_gf - equal_time_gf).lpNorm<1>()
				/ old_gf.rows() / old_gf.cols());
		}
		
		void recompute_time_displaced_gf(const dmatrix_t& U_l_,
			const dmatrix_t& D_l_, const dmatrix_t& V_l_, const dmatrix_t& U_r_,
			const dmatrix_t& D_r_, const dmatrix_t& V_r_)
		{
			dmatrix_t M(2 * id_N.rows(), 2 * id_N.cols());
			M.topLeftCorner(id_N.rows(), id_N.cols()) = V_l_.adjoint()
				* V_r_.adjoint();
			M.topRightCorner(id_N.rows(), id_N.cols()) = D_l_;
			M.bottomLeftCorner(id_N.rows(), id_N.cols()) = -D_r_;
			M.bottomRightCorner(id_N.rows(), id_N.cols()) = U_r_.adjoint()
				* U_l_.adjoint();
			dmatrix_t L = dmatrix_t::Zero(2 * id_N.rows(), 2 * id_N.cols());
			L.topLeftCorner(id_N.rows(), id_N.cols()) = V_r_.adjoint();
			L.bottomRightCorner(id_N.rows(), id_N.cols()) = U_l_.adjoint();
			dmatrix_t R = dmatrix_t::Zero(2 * id_N.rows(), 2 * id_N.cols());
			R.topLeftCorner(id_N.rows(), id_N.cols()) = V_l_.adjoint();
			R.bottomRightCorner(id_N.rows(), id_N.cols()) = U_r_.adjoint();

			svd_solver.compute(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
			dmatrix_t D = svd_solver.singularValues().unaryExpr([](double s)
				{ return 1. / s; }).asDiagonal();
			dmatrix_t GF = (L * svd_solver.matrixV()) * D * (svd_solver.matrixU().
				adjoint() * R);
			time_displaced_gf = GF.bottomLeftCorner(id_N.rows(), id_N.cols());
		}
	private:
		measurements& measure;
		bool update_time_displaced_gf;
		int n_intervals;
		dmatrix_t& equal_time_gf;
		dmatrix_t& time_displaced_gf;
		dmatrix_t id_N;
		std::vector<dmatrix_t> U;
		std::vector<dmatrix_t> D;
		std::vector<dmatrix_t> V;
		dmatrix_t U_l;
		dmatrix_t D_l;
		dmatrix_t V_l;
		dmatrix_t U_r;
		dmatrix_t D_r;
		dmatrix_t V_r;
		Eigen::JacobiSVD<dmatrix_t> svd_solver;
};
