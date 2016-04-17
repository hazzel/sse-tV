#pragma once
#include <ostream>
#include <vector>
#include "measurements.h"
#include "parser.h"
#include "configuration.h"

void eval_m2(std::valarray<double>& out,
	std::vector<std::valarray<double>*>& o, double* p)
{
	out.resize((*o[0]).size(), 0.);
	for (int t = 0; t < out.size(); ++t)
		for (int i = 0; i < p[0]; ++i)
			for (int j = 0; j < p[0]; ++j)
			{
//				out[t] += (*o[i*p[0]+j])[t] * (*o[i*p[0]+j])[t] / p[0] / p[0];
				double p_i = (i % 2 == 0) ? 1. : -1;
				double p_j = (j % 2 == 0) ? 1. : -1;
				out[t] += -p_i*p_j*(*o[j*p[0]+i])[out.size()-1-t]*(*o[i*p[0]+j])[t]
					/ p[0] / p[0];
			}
}

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
		measure.add("<k>_1", config.M.non_ident(0));
		measure.add("<k>_2", config.M.non_ident(1));
		measure.add("energy", -(config.M.non_ident(0) - config.param.beta
			* config.l.n_bonds() * config.param.t * config.param.t
			/ config.param.V1) / config.param.beta);
		config.M.measure_density_correlations(density_Corr);
		for (int i = 0; i < density_Corr.size(); ++i)
			density_Corr[i] /= config.shellsize[i];
		measure.add("<n_r n_0>", density_Corr);
	}

	void collect(std::ostream& os)
	{
		double eval_param[] = {static_cast<double>(config.l.n_sites())};
		std::vector<std::string> names;
		for (int i = 0; i < config.l.n_sites(); ++i)
			for (int j = 0; j < config.l.n_sites(); ++j)
				names.push_back("td_gf_" + std::to_string(i) + "_"
					+ std::to_string(j));
		measure.add_vectorevalable("m2_jack", names, eval_m2, eval_param);
		os << "PARAMETERS" << std::endl;
		pars.get_all(os);
		measure.get_statistics(os);
	}
};
