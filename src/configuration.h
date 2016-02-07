#pragma once
#include <vector>
#include <algorithm>
#include "measurements.h"
#include "dump.h"
#include "fast_update.h"
#include "lattice.h"
#include "greens_function.h"
#include "parameters.h"
#include "Random.h"

// The Monte Carlo configuration
struct configuration
{
	lattice l;
	parameters param;
	fast_update M;
	measurements measure;
	std::vector<int> shellsize;

	configuration(Random& rng, const greens_function& g0)
		: l{}, param{}, M{rng, l, param}
	{}

	void initialize()
	{
		shellsize.resize(l.max_distance() + 1, 0);
		for (int d = 0; d <= l.max_distance(); ++d)
		{
			int site = 0; //PBC used here
			for (int j = 0; j < l.n_sites(); ++j)
				if (l.distance(site, j) == d)
					shellsize[d] += 1;
		}
		M.initialize();
	}

	void serialize(odump& out)
	{
		M.serialize(out);
	}

	void serialize(idump& in)
	{
		M.serialize(in);
	}
};
