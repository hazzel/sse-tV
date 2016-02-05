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
	}
};
