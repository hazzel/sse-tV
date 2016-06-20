#pragma once
#include <functional>
#include <boost/algorithm/string.hpp>
#include "measurements.h"
#include "configuration.h"
#include "wick_base.h"
#include "wick_functors.h"

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

	configuration& config;
	Random& rng;
	std::vector<double> time_grid;
	std::vector<std::vector<double>> dyn_mat;
	std::vector<std::vector<double>> dyn_tau;
	std::vector<double> dyn_tau_avg;
	std::vector<wick_base<matrix_t>> obs;
	std::vector<std::string> names;

	event_dynamic_measurement(configuration& config_, Random& rng_,
//		int n_prebin, std::initializer_list<std::string> observables)
		int n_prebin, const std::vector<std::string>& observables)
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
			
			if (observables[i] == "M2")
				add_wick(wick_M2{config, rng});
			else if (observables[i] == "kekule")
				add_wick(wick_kekule{config, rng});
			else if (observables[i] == "epsilon")
				add_wick(wick_epsilon{config, rng});
			else if (observables[i] == "chern")
				add_wick(wick_chern{config, rng});
			else if (observables[i] == "sp")
				add_wick(wick_sp{config, rng});
			else if (observables[i] == "tp")
				add_wick(wick_tp{config, rng});
			
			names.push_back("dyn_"+observables[i]);
			if (config.param.n_matsubara > 0)
				config.measure.add_vectorobservable("dyn_"+observables[i]+"_mat",
					config.param.n_matsubara, n_prebin);
			if (config.param.n_discrete_tau > 0)
				config.measure.add_vectorobservable("dyn_"+observables[i]+"_tau",
					2*config.param.n_discrete_tau + 1, n_prebin);
		}
		dyn_tau_avg.resize(config.param.n_discrete_tau + 1);
	}
	
	template<typename T>
	void add_wick(T&& functor)
	{
		obs.push_back(wick_base<matrix_t>(std::forward<T>(functor)));
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
		//config.M.hirsch_fye_measurement(config.param.n_matsubara, time_grid,
		//	dyn_mat, dyn_tau, obs);

		for (int i = 0; i < dyn_mat.size(); ++i)
		{
			if (config.param.n_matsubara > 0)
				config.measure.add(names[i]+"_mat", dyn_mat[i]);
			if (config.param.n_discrete_tau > 0)
				config.measure.add(names[i]+"_tau", dyn_tau[i]);
		}
	}
};
