#pragma once
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/QR>
#include "measurements.h"
#include "dump.h"
#include "lattice.h"
#include "parameters.h"
#include "Random.h"

class qr_stabilizer
{
	public:
		template<int n, int m>
		using matrix_t = Eigen::Matrix<double, n, m>; 
		using dmatrix_t = matrix_t<Eigen::Dynamic, Eigen::Dynamic>;
		
		qr_stabilizer(measurements& measure_, dmatrix_t& equal_time_gf_,
			dmatrix_t& time_displaced_gf_)
			: measure(measure_), update_time_displaced_gf(false),
			equal_time_gf(equal_time_gf_), time_displaced_gf(time_displaced_gf_)
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
			qr_solver.compute((b * U[n-1]) * D[n-1]);
			U[n] = qr_solver.matrixQ();
			D[n] = qr_solver.matrixQR().diagonal().asDiagonal();
			V[n] = (D[n].inverse() * qr_solver.matrixQR().triangularView<Eigen
				::Upper>()) * (qr_solver.colsPermutation().transpose() * V[n-1]);
			if (n == n_intervals)
			{
				recompute_equal_time_gf(id_N, id_N, id_N, U.back(),
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
			
			qr_solver.compute((b * U[n]) * D[n]);
			U_l = U[n+1]; D_l = D[n+1]; V_l = V[n+1];
			U[n+1] = qr_solver.matrixQ();
			D[n+1] = qr_solver.matrixQR().diagonal().asDiagonal();
			V[n+1] = (D[n+1].inverse() * qr_solver.matrixQR().triangularView<Eigen
				::Upper>()) * (qr_solver.colsPermutation().transpose() * V[n]);
			
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
			
			qr_solver.compute(D[n] * (U[n] * b));
			U_r = U[n-1]; D_r = D[n-1]; V_r = V[n-1];
			V[n-1] = V[n] * qr_solver.matrixQ();
			D[n-1] = qr_solver.matrixQR().diagonal().asDiagonal();
			U[n-1] = (D[n-1].inverse() * qr_solver.matrixQR().triangularView<Eigen
				::Upper>()) * qr_solver.colsPermutation().transpose();
		
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
			qr_solver.compute(U_r_.transpose() * U_l_.inverse() + D_r_ * (V_r_
				* V_l_) * D_l_);
			dmatrix_t r = qr_solver.matrixQR().triangularView<Eigen::Upper>();
			dmatrix_t D = qr_solver.matrixQR().diagonal().asDiagonal();
			equal_time_gf = (U_l_.inverse() * (qr_solver.colsPermutation()
				* r.inverse())) * (qr_solver.matrixQ().transpose()
				* U_r_.transpose());
//			if ((old_gf - equal_time_gf).norm() > 0.00001)
//				std::cout << (old_gf - equal_time_gf).norm() << std::endl;
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
		}
	private:
		void print_matrix(const dmatrix_t& m)
		{
			Eigen::IOFormat clean(4, 0, ", ", "\n", "[", "]");
			std::cout << m.format(clean) << std::endl << std::endl;
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
		Eigen::ColPivHouseholderQR<dmatrix_t> qr_solver;
};
