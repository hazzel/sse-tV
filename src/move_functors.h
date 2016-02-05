#pragma once
#include <vector>
#include <algorithm>
#include "measurements.h"
#include "fast_update.h"
#include "lattice.h"
#include "greens_function.h"
#include "configuration.h"
#include "Random.h"

// k! / (k+n)!
template<typename T>
double factorial_ratio(T k, T n)
{
	if (k <= T(0) || n <= T(0))
		return 1.0;
	double result = 1.0;
	for (int i = 1; i <= n; ++i)
		result /= static_cast<double>(k + i);
	return result;
}

// ------------ QMC move : inserting a vertex ------------------

template<int N>
struct move_insert
{
	configuration* config;
	Random& rng;

	double attempt()
	{
		return 1.0;
	}

	double accept()
	{
		config->measure.add("insertion n="+std::to_string(N), 1.0);
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
		return 1.0;
	}

	double accept()
	{
		config->measure.add("removal n="+std::to_string(N), 1.0);
		return 1.0;
	}

	void reject()
	{
		config->measure.add("removal n="+std::to_string(N), 0.0);
	}
};
