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
	const lattice& l;
	fast_update M;
	const parameters& param;
	measurements& measure;
	std::vector<int> shellsize;

	configuration(const lattice& l_, const greens_function& g0, 
		const parameters& param_, measurements& measure_)
		: l(l_), M{l_, param_}, param(param_), measure(measure_)
	{
		shellsize.resize(l.max_distance() + 1, 0);
		for (int d = 0; d <= l.max_distance(); ++d)
		{
			int site = 0; //PBC used here
			for (int j = 0; j < l.n_sites(); ++j)
				if (l.distance(site, j) == d)
					shellsize[d] += 1;
		}
	}

	int perturbation_order() const { return M.perturbation_order(); }

	void serialize(odump& out)
	{
		M.serialize(out);
	}

	void serialize(idump& in)
	{
		M.serialize(in);
	}
};
