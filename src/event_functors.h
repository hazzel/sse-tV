#pragma once
#include "measurements.h"
#include "configuration.h"

struct event_rebuild
{
	configuration* config;
	measurements& measure;

	void trigger()
	{
		config->M.rebuild();
	}
};

struct event_build
{
	configuration* config;
	Random& rng;

	void trigger()
	{
		int n0 = 0.13 * config->params.beta * config->params.V
			* config->l.n_sites();
		std::vector<arg_t> initial_vertices;
		for (int i = 0; i < n0; ++i)
		{
			double tau = rng() * config->params.beta;
			int s1 = rng() * config->l.n_sites();
			int s2 =  config->l.neighbors(s1, "nearest neighbors")
				[rng() * config->l.neighbors(s1, "nearest neighbors").size()];
			initial_vertices.push_back({tau, s1, nn_int});
			initial_vertices.push_back({tau, s2, nn_int});
		}
		config->M.build(initial_vertices, {2*n0, 0});
	}
};
