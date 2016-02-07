#pragma once
#include <vector>
#include <algorithm>
#include "measurements.h"
#include "fast_update.h"
#include "lattice.h"
#include "greens_function.h"
#include "configuration.h"
#include "Random.h"

// ------------ QMC move : update a vertex ------------------

struct move_update_vertex
{
	configuration& config;
	Random& rng;
	int bond_type;

	double attempt()
	{
		return config.M.try_update_vertex(bond_type);
	}

	double accept()
	{
		config.measure.add("update type "+std::to_string(bond_type), 1.0);
		config.M.finish_update_vertex(bond_type);
		return 1.0;
	}

	void reject()
	{
		config.measure.add("update type "+std::to_string(bond_type), 0.0);
	}
};
