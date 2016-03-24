#pragma once
#include <functional>
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
	typedef std::function<double(const fast_update<qr_stabilizer>::dmatrix_t&,
		Random&, const lattice&, const parameters&)> function_t;

	configuration& config;
	Random& rng;
	std::vector<double> time_grid;
	std::vector<std::vector<double>> dyn_mat;
	std::vector<std::vector<double>> dyn_tau;
	std::vector<double> dyn_tau_avg;
	std::vector<function_t> obs;
	std::vector<std::string> names;

	event_dynamic_measurement(configuration& config_, Random& rng_)
		: config(config_), rng(rng_)
	{
		time_grid.resize(2 * config.param.n_discrete_tau + 1);
		for (int t = 0; t < time_grid.size(); ++t)
			time_grid[t] = static_cast<double>(t) / static_cast<double>(2*config.
				param.n_discrete_tau) * config.param.beta;
		for (int i = 0; i < 2; ++i)
		{
			dyn_mat.push_back(std::vector<double>(config.param.n_matsubara, 0.));
			dyn_tau.push_back(std::vector<double>(2 * config.param.n_discrete_tau
				+ 1, 0.));
		}
		dyn_tau_avg.resize(config.param.n_discrete_tau + 1);
		obs.emplace_back([] (const fast_update<qr_stabilizer>::dmatrix_t& gf,
			Random& rng, const lattice& l, const parameters& param)
			{
				double M2 = 0.; int i = rng() * l.n_sites();
				for (int j = 0; j < l.n_sites(); ++j)
					M2 += gf(i, j) * gf(i, j) / l.n_sites();
				return M2;
			});
		obs.emplace_back([] (const fast_update<qr_stabilizer>::dmatrix_t& gf,
			Random& rng, const lattice& l, const parameters& param)
			{
				double ep = 0.;
				for (int i = 0; i < l.n_sites(); ++i)
					for (int j : l.neighbors(i, "nearest neighbors"))
						for (int m = 0; m < l.n_sites(); ++m)
							for (int n : l.neighbors(m, "nearest neighbors"))
							{
								double d_in = (i == n) ? 1. : 0.;
								ep += gf(j, i)*gf(n, m) + (d_in - gf(n, i))*gf(j, m);
							}
				return ep;
			});
		names.push_back("dyn_M2");
		names.push_back("dyn_epsilon");
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
				dyn_tau_avg[j] = (dyn_tau[i][j] + dyn_tau[i][dyn_tau[i].size() - 1
					- j]) / 2.;
			config.measure.add(names[i]+"_tau", dyn_tau_avg);
		}
	}
};
