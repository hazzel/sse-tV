#pragma once
#include <string>
#include <triqs/h5.hpp>

struct configuration;

// ------------ The main class of the solver -----------------------

class ctint_solver
{
	unsigned int base_seed;
	unsigned int seed;
	double beta;
	int n_slices;
	int worm_nhood_dist = 4;
	double zeta2 = 20.0;
	double zeta4 = 150.0;
	double ratio_w2;
	double ratio_w4;

	public:
		ctint_solver(long unsigned int base_seed, double beta_,
			int n_slices_ = 100);
		
		void write_parameters(triqs::h5::group& group, configuration& config);

		// The method that runs the qmc
		void solve(int L, double V, int n_cycles, int length_cycle = 50,
			int n_warmup_cycles = 5000, std::string random_name = "",
			int max_time = -1);
 
};

