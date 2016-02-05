#include <string>
#include <fstream>
#include "mc.h"
#include "move_functors.h"
#include "measure_functors.h"
#include "event_functors.h"

mc::mc(const std::string& dir)
	: rng(Random()), qmc(rng)
{
	//Read parameters
	pars.read_file(dir);
	sweep = 0;
	n_cycles = pars.value_or_default<int>("cycles", 300);
	n_warmup = pars.value_or_default<int>("warmup", 100000);
	n_prebin = pars.value_or_default<int>("prebin", 500);
	n_rebuild = pars.value_or_default<int>("rebuild", 1000);
	n_matsubara = pars.value_or_default<int>("matsubara_freqs", 10);
	hc.L = pars.value_or_default<int>("L", 9);
	param.beta = 1./pars.value_or_default<double>("T", 0.2);
	param.t = pars.value_or_default<double>("t", 1.0);
	param.V = pars.value_or_default<double>("V", 1.355);
	param.lambda = std::log((2.*param.t + param.V)/(2.*param.t - param.V));

	//Proposal probabilites
	param.add = pars.value_or_default<double>("insert", 1.0);
	param.rem = pars.value_or_default<double>("remove", 1.0);

	//Initialize lattice
	lat.generate_graph(hc);
	lat.generate_neighbor_map("nearest neighbors", [this]
		(lattice::vertex_t i, lattice::vertex_t j) {
		return lat.distance(i, j) == 1; });

	//Create configuration
	config = new configuration(lat, g0, param, measure);

	//Set up Monte Carlo moves
	qmc.add_move(move_insert<1>{config, rng}, "insertion", param.add);
	qmc.add_move(move_remove<1>{config, rng}, "removal", param.rem);

	//Set up measurements

	//Measure acceptance probabilities
	measure.add_observable("sign", n_prebin * n_cycles);
	
	qmc.add_measure(measure_estimator{config, measure, pars,
		std::vector<double>(lat.max_distance() + 1, 0.0)}, "measurement");
	
	//Set up events
	qmc.add_event(event_rebuild{config, measure}, "rebuild");
	qmc.add_event(event_build{config, rng}, "initial build");
	//Initialize vertex list to reduce warm up time
	qmc.trigger_event("initial build");
}

mc::~mc()
{
	delete config;
}

void mc::random_write(odump& d)
{
	rng.RngHandle()->write(d);
}
void mc::seed_write(const std::string& fn)
{
	std::ofstream s;
	s.open(fn.c_str());
	s << rng.Seed() << std::endl;
	s.close();
}
void mc::random_read(idump& d)
{
	rng.NewRng();
	rng.RngHandle()->read(d);
}
void mc::init() {}

void mc::write(const std::string& dir)
{
	odump d(dir+"dump");
	random_write(d);
	d.write(sweep);
	config->serialize(d);
	d.close();
	seed_write(dir+"seed");
	std::ofstream f(dir+"bins");
	if (is_thermalized())
	{
		f << "Thermalization: Done." << std::endl;
		f << "Sweeps: " << (sweep - n_warmup) << std::endl;
		f << "Bins: " << static_cast<int>((sweep - n_warmup) / n_prebin)
			<< std::endl;
	}
	else
	{
		f << "Thermalization: " << sweep << std::endl;
		f << "Sweeps: 0" << std::endl;
		f << "Bins: 0" << std::endl;
	}
	f.close();
}
bool mc::read(const std::string& dir)
{
	idump d(dir+"dump");
	if (!d)
	{
		std::cout << "read fail" << std::endl;
		return false;
	}
	else
	{
		random_read(d);
		d.read(sweep);
		config->serialize(d);
		d.close();
		return true;
	}
}

void mc::write_output(const std::string& dir)
{
	std::ofstream f(dir);
	qmc.collect_results(f);
	f.close();
	/*
	const std::vector<std::pair<std::string, double>>& acc =
		qmc.acceptance_rates();
	for (auto a : acc)
		std::cout << a.first << " : " << a.second << std::endl;
	std::cout << "Average sign: " << qmc.average_sign() << std::endl;
	*/
}

bool mc::is_thermalized()
{
	return sweep >= n_warmup;
}

void mc::do_update()
{
	if (!is_thermalized())
		qmc.do_update(measure);
	else
		for (int i = 0; i < n_cycles; ++i)
			qmc.do_update(measure);
	++sweep;
	if (sweep % n_rebuild == 0)
		qmc.trigger_event("rebuild");
	status();
}

void mc::do_measurement()
{
	qmc.do_measurement();
}

void mc::status()
{
//	if (sweep == n_warmup)
//		std::cout << "Thermalization done." << std::endl;
//	if (is_thermalized() && sweep % (10000) == 0)
//	{
//		std::cout << "sweep: " << sweep << std::endl;
//		std::cout << "pert order: " << config->perturbation_order() << std::endl;
//	}
}
