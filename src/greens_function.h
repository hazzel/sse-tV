#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include "boost/multi_array.hpp"
#include "interpolation.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 
#include "lattice.h"

class greens_function
{
	public:
		typedef Eigen::MatrixXd matrix_t;
		typedef Eigen::MatrixXcd cmatrix_t;
		typedef Eigen::VectorXd vector_t;
		typedef std::complex<double> complex_t;
		greens_function() {}

		void generate_mesh(lattice* l_, double beta_, int n_time_slices_,
			int n_matsubara_)
		{
			l = l_; beta = beta_; n_time_slices = n_time_slices_;
			n_matsubara = n_matsubara_;
			dtau = beta / (2.0 * static_cast<double>(n_time_slices));
			
			matrix_t K(l->n_sites(), l->n_sites());
			for (int i = 0; i < l->n_sites(); ++i)
				for (int j = 0; j < l->n_sites(); ++j)
					K(i, j) = (l->distance(i, j) == 1) ? -1.0 : 0.0;
			Eigen::SelfAdjointEigenSolver<matrix_t> solver(K);
			vector_t ev = solver.eigenvalues();
			matrix_t V = solver.eigenvectors();
		
			generate_time_index_map(ev, V);
			generate_matsubara_index_map(ev, V);
			fill_time_mesh(ev, V);
			fill_matsubara_mesh(ev, V);
		}
		
		double imaginary_time(double tau, int i, int j) const
		{
			double tau_p;
			if (std::abs(tau) > beta/2.0)
				tau_p = beta - std::abs(tau);
			else
				tau_p = std::abs(tau);
			int x = time_index_map[i][j];
			double g = alglib::spline1dcalc(time_mesh_spline[x], tau_p);
			double sign = 1.0;
			bool same_sl = l->sublattice(i) == l->sublattice(j);
			if (std::abs(tau) > beta/2.0 && (!same_sl))
				sign *= -1.0;
			if (tau < 0.0 && same_sl)
				sign *= -1.0;
			return sign * g;
		}
		
		complex_t matsubara_frequency(int n, int i, int j) const
		{
			complex_t omega = matsubara_decimal(n);
			double re = alglib::spline1dcalc(matsubara_re_mesh_spline[n],
				std::imag(omega));
			double im = alglib::spline1dcalc(matsubara_im_mesh_spline[n],
				std::imag(omega));
			return {re, im};
		}

		complex_t matsubara_decimal(int n) const
		{
			return {0., (2.*n + 1.) * 4. * std::atan(1.) / beta};
		}
	private:
		matrix_t bare_time_gf(double tau, const vector_t& ev, const matrix_t& V)
		{
			matrix_t D = ev.unaryExpr([&](double e) { return std::exp(-tau * e)
				/ (1.0 + std::exp(-beta * e)); }).asDiagonal();
			return V * D * V.adjoint();
		}
		
		cmatrix_t bare_matsubara_gf(int n, const vector_t& ev, const matrix_t& V)
		{
			complex_t omega = {0., (2.*n + 1.) * 4. * std::atan(1.) / beta};
			cmatrix_t D = ev.cast<complex_t>().unaryExpr([&](complex_t e)
				{ return 1.0 / (omega - e); }).asDiagonal();
			return V.cast<complex_t>() * D * V.adjoint().cast<complex_t>();
		}
		
		void generate_time_index_map(const vector_t& ev, const matrix_t& V)
		{
			time_index_map.resize(boost::extents[l->n_sites()][l->n_sites()]);
			double threshold = std::pow(10.0, -13.0);
			std::vector<double> values;
			
			matrix_t g0 = bare_time_gf(0.1 * beta, ev, V);
			for (int i = 0; i < l->n_sites(); ++i)
			{
				for (int j = i; j < l->n_sites(); ++j)
				{
					bool is_stored = false;
					for (int k = 0; k < values.size(); ++k)
					{
						if (std::abs(values[k] - g0(i, j)) < threshold)
						{
							is_stored = true;
							time_index_map[i][j] = k;
							time_index_map[j][i] = k;
						}
					}
					if (!is_stored)
					{
						values.push_back(g0(i, j));
						time_index_map[i][j] = values.size() - 1;
						time_index_map[j][i] = values.size() - 1;
					}
				}
			}
			//mesh.resize(boost::extents[values.size()][n_time_slices + 1]);
			time_mesh_y.resize(values.size());
			for (auto& y : time_mesh_y)
				y.setlength(n_time_slices + 1);
			time_mesh_spline.resize(values.size());
		}
		
		void generate_matsubara_index_map(const vector_t& ev, const matrix_t& V)
		{
			matsubara_index_map.resize(boost::extents[l->n_sites()][l->n_sites()]);
			double threshold = std::pow(10.0, -11.0);
			std::vector<complex_t> values;
			
			cmatrix_t g0 = bare_matsubara_gf(1, ev, V);
			for (int i = 0; i < l->n_sites(); ++i)
			{
				for (int j = i; j < l->n_sites(); ++j)
				{
					bool is_stored = false;
					for (int k = 0; k < values.size(); ++k)
					{
						if (std::abs(values[k] - g0(i, j)) < threshold)
						{
							is_stored = true;
							matsubara_index_map[i][j] = k;
							matsubara_index_map[j][i] = k;
						}
					}
					if (!is_stored)
					{
						values.push_back(g0(i, j));
						matsubara_index_map[i][j] = values.size() - 1;
						matsubara_index_map[j][i] = values.size() - 1;
					}
				}
			}
			matsubara_re_mesh_y.resize(values.size());
			matsubara_im_mesh_y.resize(values.size());
			for (auto& y : matsubara_re_mesh_y)
				y.setlength(n_matsubara);
			for (auto& y : matsubara_im_mesh_y)
				y.setlength(n_matsubara);
			matsubara_re_mesh_spline.resize(values.size());
			matsubara_im_mesh_spline.resize(values.size());
		}
		
		void fill_time_mesh(const vector_t& ev, const matrix_t& V)
		{
			// Imaginary time mesh
			alglib::real_1d_array time_mesh_x;
			time_mesh_x.setlength(n_time_slices + 1);
			for (int t = 0; t <= n_time_slices; ++t)
			{
				time_mesh_x[t] = dtau * t;
				matrix_t g0 = bare_time_gf(dtau * t, ev, V);
				for (int i = 0; i < l->n_sites(); ++i)
					for (int j = i; j < l->n_sites(); ++j)
						time_mesh_y[time_index_map[i][j]][t] = g0(i, j);
			}
			for (int i = 0; i < time_mesh_y.size(); ++i)
				alglib::spline1dbuildakima(time_mesh_x, time_mesh_y[i],
					time_mesh_spline[i]);
		}

		void fill_matsubara_mesh(const vector_t& ev, const matrix_t& V)
		{
			// Matsubara frequency mesh
			alglib::real_1d_array matsubara_mesh_x;
			matsubara_mesh_x.setlength(n_matsubara);
			for (int n = 0; n < n_matsubara; ++n)
			{
				matsubara_mesh_x[n] = (2.*n + 1.) * 4. * std::atan(1.) / beta;
				cmatrix_t g0 = bare_matsubara_gf(n, ev, V);
				for (int i = 0; i < l->n_sites(); ++i)
					for (int j = i; j < l->n_sites(); ++j)
					{
						matsubara_re_mesh_y[matsubara_index_map[i][j]][n] = 
							std::real(g0(i, j));
						matsubara_im_mesh_y[matsubara_index_map[i][j]][n] = 
							std::imag(g0(i, j));
					}
			}
			for (int i = 0; i < matsubara_re_mesh_y.size(); ++i)
			{
				alglib::spline1dbuildakima(matsubara_mesh_x, matsubara_re_mesh_y[i],
					matsubara_re_mesh_spline[i]);
				alglib::spline1dbuildakima(matsubara_mesh_x, matsubara_im_mesh_y[i],
					matsubara_im_mesh_spline[i]);
			}
		}
	private:
		lattice* l;
		double beta;
		double dtau;
		int n_time_slices;
		int n_matsubara;
		boost::multi_array<int, 2> time_index_map;
		std::vector<alglib::real_1d_array> time_mesh_y;
		std::vector<alglib::spline1dinterpolant> time_mesh_spline;
		boost::multi_array<int, 2> matsubara_index_map;
		std::vector<alglib::real_1d_array> matsubara_re_mesh_y;
		std::vector<alglib::real_1d_array> matsubara_im_mesh_y;
		std::vector<alglib::spline1dinterpolant> matsubara_re_mesh_spline;
		std::vector<alglib::spline1dinterpolant> matsubara_im_mesh_spline;
};
