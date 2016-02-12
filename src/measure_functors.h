#pragma once
#include <ostream>
#include <vector>
#include "measurements.h"
#include "parser.h"
#include "configuration.h"

struct measure_estimator
{
	configuration& config;
	measurements& measure;
	parser& pars;
	std::vector<double> density_Corr;
	std::vector<double> matsubara_G;

	void perform()
	{
		measure.add("M2", config.M.measure_M2());
		measure.add("<n1>", config.M.non_ident());
		std::fill(density_Corr.begin(), density_Corr.end(), 0.0);
		config.M.measure_density_correlations(density_Corr);
		for (int i = 0; i < density_Corr.size(); ++i)
			density_Corr[i] /= config.shellsize[i];
		measure.add("<n_r n_0>", density_Corr);
	}

	void collect(std::ostream& os)
	{
		os << "PARAMETERS" << std::endl;
		pars.get_all(os);
		measure.get_statistics(os);
	}
};
