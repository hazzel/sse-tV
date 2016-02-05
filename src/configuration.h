#pragma once
#include <vector>
#include <algorithm>
#include "measurements.h"
#include "dump.h"
#include "fast_update.h"
#include "lattice.h"
#include "greens_function.h"
#include "Random.h"

// Argument type
struct arg_t
{
	double tau;
	int site;
	bool worm;

	void serialize(odump& out)
	{
		out.write(tau);
		out.write(site);
		out.write(worm);
	}

	void serialize(idump& in)
	{
		double t; in.read(t); tau = t;
		int s; in.read(s); site = s;
		bool w; in.read(w); worm = w;
	}
};

enum { nn_int, worm };

// The function that appears in the calculation of the determinant
struct full_g_entry
{
	const greens_function& g0;

	double operator()(const arg_t& x, const arg_t& y) const
	{
		return g0.imaginary_time(x.tau - y.tau, x.site, y.site);
	}

	std::complex<double> matsubara_frequency(int omega_n, const arg_t& x,
		const arg_t& y) const
	{
		return g0.matsubara_frequency(omega_n, x.site, y.site)
			* std::exp(g0.matsubara_decimal(omega_n) * (x.tau - y.tau));
	}
};

struct parameters
{	
	double beta, V, zeta2, zeta4;
	int worm_nhood_dist;
	double ratio_w2, ratio_w4;
	//Proposal probabilities
	std::vector<double> add;
	std::vector<double> rem;
	double W2toZ, ZtoW2, ZtoW4, W4toZ, W2toW4, W4toW2, worm_shift;
};

// The Monte Carlo configuration
struct configuration
{
	const lattice& l;
	fast_update<full_g_entry, arg_t> M;
	const parameters& params;
	measurements& measure;
	std::vector<int> shellsize;

	configuration(const lattice& l_, const greens_function& g0, 
		const parameters& params_, measurements& measure_)
		: l(l_), M{full_g_entry{g0}, l_, 2}, params(params_), measure(measure_)
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

	int perturbation_order() const { return M.perturbation_order(nn_int); }
	int worms() const { return M.perturbation_order(worm); }

	void serialize(odump& out)
	{
		M.serialize(out);
	}

	void serialize(idump& in)
	{
		M.serialize(in);
	}
};
