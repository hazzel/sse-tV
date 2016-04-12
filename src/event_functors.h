#pragma once
#include <functional>
#include <boost/algorithm/string.hpp>
#include "measurements.h"
#include "configuration.h"

struct event_rebuild
{
	configuration& config;
	measurements& measure;

	void trigger()
	{
		config.M.rebuild();
	}
};

struct event_build
{
	configuration& config;
	Random& rng;

	void trigger()
	{}
};

struct event_max_order
{
	configuration& config;
	Random& rng;

	void trigger()
	{
//		if (config.M.non_ident(0)+config.M.non_ident(1) >= 0.75
//			* config.M.max_order())
//			config.M.max_order(std::max(10., 4./3. * config.M.max_order()));
//		else if (config.M.non_ident(0) + config.M.non_ident(1) <= 0.5
//			* config.M.max_order())
//			config.M.max_order(std::max(10., 0.8 * config.M.max_order()));
	}
};

struct event_dynamic_measurement
{
	typedef fast_update<qr_stabilizer>::dmatrix_t matrix_t;
	typedef std::function<double(const matrix_t&, const matrix_t&, Random&,
		const lattice&, const parameters&)> function_t;

	configuration& config;
	Random& rng;
	std::vector<double> time_grid;
	std::vector<std::vector<double>> dyn_mat;
	std::vector<std::vector<double>> dyn_tau;
	std::vector<double> dyn_tau_avg;
	std::vector<function_t> obs;
	std::vector<std::string> names;

	event_dynamic_measurement(configuration& config_, Random& rng_, int n_prebin,
		std::initializer_list<std::string> observables)
		: config(config_), rng(rng_)
	{
		time_grid.resize(2 * config.param.n_discrete_tau + 1);
		for (int t = 0; t < time_grid.size(); ++t)
			time_grid[t] = static_cast<double>(t) / static_cast<double>(2*config.
				param.n_discrete_tau) * config.param.beta;
		for (int i = 0; i < observables.size(); ++i)
		{
			dyn_mat.push_back(std::vector<double>(config.param.n_matsubara, 0.));
			dyn_tau.push_back(std::vector<double>(2 * config.param.n_discrete_tau
				+ 1, 0.));
		}
		dyn_tau_avg.resize(config.param.n_discrete_tau + 1);
		typedef std::initializer_list<std::string> list_t;
		if (boost::algorithm::contains(observables, list_t{"M2"}))
		{
			// M2(tau) = sum_ij <(n_i(tau) - 1/2)(n_j - 1/2)>
			obs.emplace_back([] (const matrix_t& equal_time_gf, const matrix_t&
				time_displaced_gf, Random& rng, const lattice& l,
				const parameters& param)
				{
					double M2 = 0.; int i = rng() * l.n_sites();
					for (int j = 0; j < l.n_sites(); ++j)
						M2 += time_displaced_gf(i, j) * time_displaced_gf(i, j)
							/ l.n_sites();
					return M2;
				});
			names.push_back("dyn_M2");
			config.measure.add_vectorobservable("dyn_M2_mat",
				config.param.n_matsubara, n_prebin);
			config.measure.add_vectorobservable("dyn_M2_tau",
				config.param.n_discrete_tau + 1, n_prebin);
		}
		if (boost::algorithm::contains(observables, list_t{"epsilon"}))
		{
			// ep(tau) = sum_{<ij>,<mn>} <c_i^dag(tau) c_j(tau) c_n^dag c_m>
			obs.emplace_back([] (const matrix_t& equal_time_gf, const matrix_t&
				time_displaced_gf, Random& rng, const lattice& l,
				const parameters& param)
				{
					double ep = 0.;
					int i = rng() * l.n_sites();
					for (int j : l.neighbors(i, "nearest neighbors"))	
						for (int m = 0; m < l.n_sites(); ++m)
							for (int n : l.neighbors(m, "nearest neighbors"))
							{
								ep += time_displaced_gf(m, i)
									* time_displaced_gf(j, n) / l.n_bonds() * 2./3.;
							}
					return ep;
				});
			names.push_back("dyn_epsilon");
			config.measure.add_vectorobservable("dyn_epsilon_mat",
				config.param.n_matsubara, n_prebin);
			config.measure.add_vectorobservable("dyn_epsilon_tau",
				config.param.n_discrete_tau + 1, n_prebin);
		}
		if (boost::algorithm::contains(observables, list_t{"sp"}))
		{
			// sp(tau) = sum_ij e^{-i K (r_i - r_j)} <c_i(tau) c_j^dag>
			obs.emplace_back([] (const matrix_t& equal_time_gf, const matrix_t&
				time_displaced_gf, Random& rng, const lattice& l,
				const parameters& param)
				{
					double sp = 0.;
					double pi = 4.*std::atan(1.);
					Eigen::Vector2d K(2.*pi/9., 2.*pi/9.*(2.-1./std::sqrt(3.)));
					int i = rng() * l.n_sites();
					for (int j = 0; j < l.n_sites(); ++j)
					{
						auto& r_i = l.real_space_coord(i);
						auto& r_j = l.real_space_coord(j);
//						sp +=	std::cos(K.dot(r_j - r_i)) * l.parity(i) * l.parity(j)
//							* time_displaced_gf(j, i) * l.n_sites();
						sp +=	std::cos(K.dot(r_j - r_i))
							* time_displaced_gf(i, j) * l.n_sites();
					}
					return sp;
				});
			names.push_back("dyn_sp");
			config.measure.add_vectorobservable("dyn_sp_mat",
				config.param.n_matsubara, n_prebin);
			config.measure.add_vectorobservable("dyn_sp_tau",
				2*config.param.n_discrete_tau + 1, n_prebin);
		}
		if (boost::algorithm::contains(observables, list_t{"tp"}))
		{
			// sp(tau) = sum_ij e^{-i K (r_i - r_j)} <c_i(tau) c_j^dag>
			obs.emplace_back([] (const matrix_t& equal_time_gf, const matrix_t&
				time_displaced_gf, Random& rng, const lattice& l,
				const parameters& param)
				{
					double tp = 0.;
					double pi = 4.*std::atan(1.);
					Eigen::Vector2d K(2.*pi/9., 2.*pi/9.*(2.-1./std::sqrt(3.)));
					auto wick = [&] (int i, int j, int m, int n)->double
					{
						auto& r_i = l.real_space_coord(i);
						auto& r_j = l.real_space_coord(j);
						auto& r_m = l.real_space_coord(m);
						auto& r_n = l.real_space_coord(n);
						return std::cos(K.dot(r_j - r_i + r_m - r_n))
							* (time_displaced_gf(i, m) * time_displaced_gf(j, n)
							- time_displaced_gf(i, n) * time_displaced_gf(j, m));
					};
					/*
					int i = rng() * l.n_sites();
					for (int j = 0; j < l.n_sites(); ++j)
						for (int m = j+1; m < l.n_sites(); ++m)
							for (int n = m+1; n < l.n_sites(); ++n)
								tp += 6. * wick(i, j, m, n);
					for (int n = 0; n < l.n_sites(); ++n)
					{
						int j = 0, m = 0;
						tp += (3.*l.n_sites()-2) * wick(i, j, m, n);
					}
					*/
					int i = rng() * l.n_sites();
					for (int j = 0; j < l.n_sites(); ++j)
						for (int m = 0; m < l.n_sites(); ++m)
							for (int n = 0; n < l.n_sites(); ++n)
								tp += wick(i, j, m, n) * l.n_sites();
					return tp;
				});
			names.push_back("dyn_tp");
			config.measure.add_vectorobservable("dyn_tp_mat",
				config.param.n_matsubara, n_prebin);
			config.measure.add_vectorobservable("dyn_tp_tau",
				config.param.n_discrete_tau + 1, n_prebin);
		}
	}

	void trigger()
	{
		for (int i = 0; i < dyn_mat.size(); ++i)
		{
			std::fill(dyn_mat[i].begin(), dyn_mat[i].end(), 0.);
			std::fill(dyn_tau[i].begin(), dyn_tau[i].end(), 0.);
		}
		config.M.measure_dynamical_observable(config.param.n_matsubara, time_grid,
			dyn_mat, dyn_tau, obs);

		for (int i = 0; i < dyn_mat.size(); ++i)
		{
			config.measure.add(names[i]+"_mat", dyn_mat[i]);
			// Average imaginary time measurements from 0..beta/2 and beta/2..beta
			for (int j = 0; j < dyn_tau_avg.size(); ++j)
//				dyn_tau_avg[j] = (dyn_tau[i][j] + dyn_tau[i][dyn_tau[i].size() - 1
//					- j]) / 2.;
//				dyn_tau_avg[j] = dyn_tau[i][dyn_tau[i].size() - 1 - j];
//				dyn_tau_avg[j] = dyn_tau[i][j];
			config.measure.add(names[i]+"_tau", dyn_tau[i]);
		}
	}
};
