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
	std::vector<double> matsubara_G;

	void perform()
	{
		measure.add("M2", config.M.measure_M2());
		measure.add("<n1>", config.M.non_ident());
	}

	void collect(std::ostream& os)
	{
		os << "PARAMETERS" << std::endl;
		pars.get_all(os);
		measure.get_statistics(os);
	}
};
