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

// Argument type
struct arg_t
{
	std::pair<int, int> bond;
	bool unity;

	void serialize(odump& out)
	{
		out.write(bond.first);
		out.write(bond.second);
		out.write(unity);
	}

	void serialize(idump& in)
	{
		int s; in.read(s); bond.first = s;
		in.read(s); bond.second = s;
		bool u; in.read(u); unity = u;
	}
};

struct vertex_matrix
{
	const parameters& param;

	double operator()(const arg_t& x, const arg_t& y) const
	{
		return param.lambda; 
	}
};

// The Monte Carlo configuration
struct configuration
{
	const lattice& l;
	fast_update<vertex_matrix, arg_t> M;
	const parameters& param;
	measurements& measure;
	std::vector<int> shellsize;

	configuration(const lattice& l_, const greens_function& g0, 
		const parameters& param_, measurements& measure_)
		: l(l_), M{vertex_matrix{param_}, l_, param_}, param(param_),
			measure(measure_)
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
