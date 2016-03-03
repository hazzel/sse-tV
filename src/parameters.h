#pragma once
#include <complex>

struct parameters
{
	double beta, t, V, V1, V2, lambda;
	double prop_V1, prop_V2;
	int n_delta, n_discrete_tau, n_matsubara;
};
