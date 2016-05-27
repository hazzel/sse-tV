#pragma once
#include <ostream>
#include <vector>
#include <Eigen/Dense>
#include "measurements.h"
#include "parser.h"
#include "configuration.h"

struct measure_estimator
{
	configuration& config;
	parser& pars;
	std::vector<double> density_Corr;
	std::vector<double> matsubara_G;

	void perform()
	{
		config.measure.add("M2", config.M.measure_M2());
		config.measure.add("epsilon", config.M.measure_epsilon());
		config.measure.add("kekule", config.M.measure_kekule());
		config.measure.add("<k>_1", config.M.non_ident(0));
		config.measure.add("<k>_2", config.M.non_ident(1));
		config.measure.add("energy", -(config.M.non_ident(0) - config.param.beta
			* config.l.n_bonds() * config.param.t * config.param.t
			/ config.param.V1) / config.param.beta);
		config.M.measure_density_correlations(density_Corr);
		for (int i = 0; i < density_Corr.size(); ++i)
			density_Corr[i] /= config.shellsize[i];
		config.measure.add("<n_r n_0>", density_Corr);
	}

	void collect(std::ostream& os)
	{
		os << "PARAMETERS" << std::endl;
		pars.get_all(os);
		config.measure.get_statistics(os);
	}
};
