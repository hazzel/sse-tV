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
		if (config.M.non_ident() >= 0.75 * config.M.max_order())
			config.M.max_order(std::max(10., 4./3. * config.M.max_order()));
		else if (config.M.non_ident() <= 0.5 * config.M.max_order())
			config.M.max_order(std::max(10., 0.8 * config.M.max_order()));
	}
};
