#pragma once
#include <vector>
#include <algorithm>
#include "measurements.h"
#include "fast_update.h"
#include "lattice.h"
#include "greens_function.h"
#include "configuration.h"
#include "Random.h"

// ------------ QMC move : inserting a vertex ------------------

template<int N>
struct move_insert
{
	configuration* config;
	Random& rng;

	double attempt()
	{
		return config->M.try_insert_bond();
	}

	double accept()
	{
		config->measure.add("insertion n="+std::to_string(N), 1.0);
		config->M.finish_insert_bond();
		return 1.0;
	}

	void reject()
	{
		config->measure.add("insertion n="+std::to_string(N), 0.0);
	}
};

// ------------ QMC move : deleting a vertex ------------------

template<int N>
struct move_remove
{
	configuration* config;
	Random& rng;

	double attempt()
	{
		return config->M.try_remove_bond();
	}

	double accept()
	{
		config->measure.add("removal n="+std::to_string(N), 1.0);
		config->M.finish_remove_bond();
		return 1.0;
	}

	void reject()
	{
		config->measure.add("removal n="+std::to_string(N), 0.0);
	}
};
