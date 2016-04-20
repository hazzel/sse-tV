#pragma once
#include <vector>
#include <functional>
#include <utility>
#include <memory>
#include <ostream>
#include <iostream>
#include "measurements.h"
#include "configuration.h"

typedef fast_update<qr_stabilizer>::dmatrix_t matrix_t;

// M2(tau) = sum_ij <(n_i(tau) - 1/2)(n_j - 1/2)>
struct wick_M2
{
	configuration& config;
	Random& rng;

	wick_M2(configuration& config_, Random& rng_)
		: config(config_), rng(rng_)
	{}
	
	double get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf_t, const matrix_t& td_gf_mt)
	{
		double M2 = 0.;
		for (int i = 0; i < config.l.n_sites(); ++i)
			for (int j = 0; j < config.l.n_sites(); ++j)
				M2 += td_gf_t(i, j) * td_gf_t(i, j)
					/ std::pow(config.l.n_sites(), 2.);
		return M2;
	}
};

// kekule(tau) = sum_{kekule} <c_i^dag(tau) c_j(tau) c_n^dag c_m>
struct wick_kekule
{
	configuration& config;
	Random& rng;

	wick_kekule(configuration& config_, Random& rng_)
		: config(config_), rng(rng_)
	{}
	
	double get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf_t, const matrix_t& td_gf_mt)
	{
		double kek = 0.;
		for (auto& a : config.l.bonds("kekule"))
			for (auto& b : config.l.bonds("kekule"))
			{
				kek += (et_gf_0(a.second, a.first) * et_gf_0(b.first,
					b.second) - td_gf_t(b.first, a.first)
					* td_gf_t(a.second, b.second));
			}
		return kek;
	}
};

// ep(tau) = sum_{<ij>,<mn>} <c_i^dag(tau) c_j(tau) c_n^dag c_m>
struct wick_epsilon
{
	configuration& config;
	Random& rng;

	wick_epsilon(configuration& config_, Random& rng_)
		: config(config_), rng(rng_)
	{}
	
	double get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf_t, const matrix_t& td_gf_mt)
	{
		double ep = 0.;
		for (int i = 0; i < config.l.n_sites(); ++i)
			for (int j : config.l.neighbors(i, "nearest neighbors"))	
				for (int m = 0; m < config.l.n_sites(); ++m)
					for (int n : config.l.neighbors(m, "nearest neighbors"))
					{
						ep += (et_gf_0(j, i) * et_gf_0(m, n)
							+ config.l.parity(i) * config.l.parity(m) * td_gf_t(i, m)
								* td_gf_t(j, n)) / std::pow(config.l.n_bonds(), 2.);
					}
		return ep;
	}
};

// sp(tau) = sum_ij e^{-i K (r_i - r_j)} <c_i(tau) c_j^dag>
struct wick_sp
{
	configuration& config;
	Random& rng;

	wick_sp(configuration& config_, Random& rng_)
		: config(config_), rng(rng_)
	{}
	
	double get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf_t, const matrix_t& td_gf_mt)
	{
		double sp = 0.;
		double pi = 4.*std::atan(1.);
		Eigen::Vector2d K(2.*pi/9., 2.*pi/9.*(2.-1./std::sqrt(3.)));
		for (int i = 0; i < config.l.n_sites(); ++i)
			for (int j = 0; j < config.l.n_sites(); ++j)
			{
				auto& r_i = config.l.real_space_coord(i);
				auto& r_j = config.l.real_space_coord(j);
				sp +=	std::cos(K.dot(r_j - r_i)) * td_gf_t(i, j);
			}
		return sp;
	}
};

// tp(tau) = sum_ijmn e^{-i K (r_i - r_j + r_m - r_n)}
			//		<c_i(tau) c_j(tau) c_n^dag c_m^dag>
struct wick_tp
{
	configuration& config;
	Random& rng;

	wick_tp(configuration& config_, Random& rng_)
		: config(config_), rng(rng_)
	{}
	
	double get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf_t, const matrix_t& td_gf_mt)
	{
		double tp = 0.;
		double pi = 4.*std::atan(1.);
		Eigen::Vector2d K(2.*pi/9., 2.*pi/9.*(2.-1./std::sqrt(3.)));
		for (int i = 0; i < config.l.n_sites(); ++i)
			for (int j = 0; j < config.l.n_sites(); ++j)
				for (int m = 0; m < config.l.n_sites(); ++m)
					for (int n = 0; n < config.l.n_sites(); ++n)
					{
						auto& r_i = config.l.real_space_coord(i);
						auto& r_j = config.l.real_space_coord(j);
						auto& r_m = config.l.real_space_coord(m);
						auto& r_n = config.l.real_space_coord(n);
						tp += std::cos(K.dot(r_j - r_i + r_m - r_n))
							* (td_gf_t(i, m) * td_gf_t(j, n)
							- td_gf_t(i, n) * td_gf_t(j, m));
					}
		return tp;
	}
};
