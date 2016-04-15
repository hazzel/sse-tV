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
	
	double get_obs(const matrix_t& equal_time_gf,
		const matrix_t& time_displaced_gf)
	{
		double M2 = 0.; int i = rng() * config.l.n_sites();
		for (int j = 0; j < config.l.n_sites(); ++j)
			M2 += time_displaced_gf(i, j) * time_displaced_gf(i, j)
				/ config.l.n_sites();
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
	
	double get_obs(const matrix_t& equal_time_gf,
		const matrix_t& time_displaced_gf)
	{
		double kek = 0.;
		for (auto& a : config.l.bonds("kekule"))
			for (auto& b : config.l.bonds("kekule"))
			{
				kek += (equal_time_gf(a.second, a.first) * equal_time_gf(b.first,
					b.second) - time_displaced_gf(b.first, a.first)
					* time_displaced_gf(a.second, b.second));
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
	
	double get_obs(const matrix_t& equal_time_gf,
		const matrix_t& time_displaced_gf)
	{
		double ep = 0.;
		int i = rng() * config.l.n_sites();
		for (int j : config.l.neighbors(i, "nearest neighbors"))	
			for (int m = 0; m < config.l.n_sites(); ++m)
				for (int n : config.l.neighbors(m, "nearest neighbors"))
				{
					ep += (equal_time_gf(j, i) * equal_time_gf(m, n)
						- time_displaced_gf(m, i) * time_displaced_gf(j, n))
						/ config.l.n_bonds() * 2./3.;
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
	
	double get_obs(const matrix_t& equal_time_gf,
		const matrix_t& time_displaced_gf)
	{
		double sp = 0.;
		double pi = 4.*std::atan(1.);
		Eigen::Vector2d K(2.*pi/9., 2.*pi/9.*(2.-1./std::sqrt(3.)));
		int i = rng() * config.l.n_sites();
		for (int j = 0; j < config.l.n_sites(); ++j)
		{
			auto& r_i = config.l.real_space_coord(i);
			auto& r_j = config.l.real_space_coord(j);
//						sp +=	std::cos(K.dot(r_j - r_i)) * config.l.parity(i) * config.l.parity(j)
//							* time_displaced_gf(j, i) * config.l.n_sites();
			sp +=	std::cos(K.dot(r_j - r_i))
				* time_displaced_gf(i, j) * config.l.n_sites();
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
	
	double get_obs(const matrix_t& equal_time_gf,
		const matrix_t& time_displaced_gf)
	{
		double tp = 0.;
		double pi = 4.*std::atan(1.);
		Eigen::Vector2d K(2.*pi/9., 2.*pi/9.*(2.-1./std::sqrt(3.)));
		auto wick = [&] (int i, int j, int m, int n)->double
		{
			auto& r_i = config.l.real_space_coord(i);
			auto& r_j = config.l.real_space_coord(j);
			auto& r_m = config.l.real_space_coord(m);
			auto& r_n = config.l.real_space_coord(n);
			return std::cos(K.dot(r_j - r_i + r_m - r_n))
				* (time_displaced_gf(i, m) * time_displaced_gf(j, n)
				- time_displaced_gf(i, n) * time_displaced_gf(j, m));
		};
		int i = rng() * config.l.n_sites();
		for (int j = 0; j < config.l.n_sites(); ++j)
			for (int m = j+1; m < config.l.n_sites(); ++m)
				for (int n = m+1; n < config.l.n_sites(); ++n)
					tp += 6. * wick(i, j, m, n);
		for (int n = 0; n < config.l.n_sites(); ++n)
		{
			int j = 0, m = 0;
			tp += (3.*config.l.n_sites()-2) * wick(i, j, m, n);
		}
		/*
		int i = rng() * config.l.n_sites();
		for (int j = 0; j < config.l.n_sites(); ++j)
			for (int m = 0; m < config.l.n_sites(); ++m)
				for (int n = 0; n < config.l.n_sites(); ++n)
					tp += wick(i, j, m, n) * l.n_sites();
		*/
		return tp;
	}
};
