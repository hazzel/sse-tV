#pragma once
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

struct event_dyn_M2
{
	configuration& config;
	Random& rng;

	void trigger()
	{
		std::vector<double> dyn_M2_mat(config.param.n_matsubara, 0.);
		std::vector<double> time_grid(2 * config.param.n_discrete_tau + 1);
		for (int t = 0; t < time_grid.size(); ++t)
			time_grid[t] = static_cast<double>(t) / static_cast<double>(2*config.
				param.n_discrete_tau) * config.param.beta;
		std::vector<double> dyn_M2_tau(time_grid.size(), 0.);
		config.M.measure_dynamical_observable(config.param.n_matsubara,
			dyn_M2_mat, time_grid, dyn_M2_tau,
			[&] (const fast_update<qr_stabilizer>::dmatrix_t& gf)
			{
				double M2 = 0.; int i = rng() * config.l.n_sites();
				for (int j = 0; j < config.l.n_sites(); ++j)
					M2 += gf(i, j) * gf(i, j) / config.l.n_sites();
				return M2;
			});
		std::vector<double> dyn_M2_tau_avg(config.param.n_discrete_tau + 1);
		for (int i = 0; i < dyn_M2_tau_avg.size(); ++i)
			dyn_M2_tau_avg[i] = (dyn_M2_tau[i] + dyn_M2_tau[dyn_M2_tau.size()
				- 1 - i]) / 2.;
		config.measure.add("dynamical_M2_mat", dyn_M2_mat);
		config.measure.add("dynamical_M2_tau", dyn_M2_tau_avg);
	}
};
