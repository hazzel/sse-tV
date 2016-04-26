#pragma once
#include <vector>
#include <algorithm>
#include "measurements.h"
#include "dump.h"
#include "svd_stabilizer.h"
#include "qr_stabilizer.h"
#include "fast_update.h"
#include "lattice.h"
#include "parameters.h"
#include "Random.h"
#include "spline.h"

// The Monte Carlo configuration
struct configuration
{
	lattice l;
	parameters param;
	measurements& measure;
	spline trig_spline;
	fast_update<qr_stabilizer> M;
	std::vector<int> shellsize;

	configuration(Random& rng, measurements& measure_)
		: l(), param(), measure(measure_), trig_spline(1000),
			M{rng, l, param, measure, trig_spline}
	{}

	void initialize(int max_order)
	{
		shellsize.resize(l.max_distance() + 1, 0);
		for (int d = 0; d <= l.max_distance(); ++d)
		{
			int site = 0; //PBC used here
			for (int j = 0; j < l.n_sites(); ++j)
				if (l.distance(site, j) == d)
					shellsize[d] += 1;
		}
		M.initialize(max_order);
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
