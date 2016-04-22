#include <string>
#include <fstream>
#include "mc.h"
#include "move_functors.h"
#include "measure_functors.h"
#include "event_functors.h"

mc::mc(const std::string& dir)
	: rng(Random()), qmc(rng), config(rng, measure)
{
	//Read parameters
	pars.read_file(dir);
	sweep = 0;
	measure_cnt = 0;
	n_cycles = pars.value_or_default<int>("cycles", 300);
	n_warmup = pars.value_or_default<int>("warmup", 100000);
	n_prebin = pars.value_or_default<int>("prebin", 500);
	n_rebuild = pars.value_or_default<int>("rebuild", 1000);
	config.param.n_discrete_tau = pars.value_or_default<int>("discrete_tau", 10);
	config.param.n_matsubara = pars.value_or_default<int>("matsubara_freqs", 10);
	hc.L = pars.value_or_default<int>("L", 9);
	config.param.beta = 1./pars.value_or_default<double>("T", 0.2);
	config.param.t = pars.value_or_default<double>("t", 1.0);
	config.param.V = pars.value_or_default<double>("V", 1.355);
	double v1_cutoff = pars.value_or_default<double>("v1_cutoff", 0.5);
	if (config.param.V > v1_cutoff)
	{
		config.param.V1 = v1_cutoff;
		config.param.V2 = config.param.V - v1_cutoff;
	}
	else
	{
		config.param.V1 = config.param.V;
		config.param.V2 = 0;
	}
	config.param.lambda = std::log((2.*config.param.t + config.param.V1)
		/ (2.*config.param.t - config.param.V1));
	config.param.n_delta = pars.value_or_default<int>("stabilization", 10);
	int max_order = pars.value_or_default<int>("max_order", 5000);

	//Proposal probabilites
	config.param.prop_V1 = pars.value_or_default<double>("prop_V1", 1.0);
	config.param.prop_V2 = pars.value_or_default<double>("prop_V2", 0.0);

	//Initialize lattice
	config.l.generate_graph(hc);
	hc.generate_maps(config.l);

	//Set up Monte Carlo moves
	qmc.add_move(move_update_vertex{config, rng, 0}, "update type 0",
		config.param.prop_V1);
	qmc.add_move(move_update_vertex{config, rng, 1}, "update type 1",
		config.param.prop_V2);
	
	//Measure acceptance probabilities
	config.measure.add_observable("update type 0", n_prebin * n_cycles);
	config.measure.add_observable("update type 1", n_prebin * n_cycles);
	config.measure.add_observable("sign", n_prebin * n_cycles);

	//Set up measurements
	config.measure.add_observable("M2", n_prebin);
	config.measure.add_observable("epsilon", n_prebin);
	config.measure.add_observable("<k>_1", n_prebin);
	config.measure.add_observable("<k>_2", n_prebin);
	config.measure.add_observable("energy", n_prebin);
	config.measure.add_observable("norm error", n_prebin);
	config.measure.add_observable("max error", n_prebin);
	config.measure.add_observable("avg error", n_prebin);
	config.measure.add_vectorobservable("<n_r n_0>", config.l.max_distance() + 1,
		n_prebin);
	
	qmc.add_measure(measure_estimator{config, pars,
		std::vector<double>(config.l.max_distance() + 1, 0.0)}, "measurement");
	
	//Initialize configuration class
	config.initialize(max_order);
	
	//Set up events
	qmc.add_event(event_rebuild{config, config.measure}, "rebuild");
	qmc.add_event(event_build{config, rng}, "initial build");
	qmc.add_event(event_max_order{config, rng}, "max_order");
	qmc.add_event(event_dynamic_measurement{config, rng, n_prebin,
		{"M2", "sp"}}, "dyn_measure");
//		{"M2", "epsilon", "sp", "tp"}}, "dyn_measure");
//		{"M2", "kekule", "epsilon", "sp", "tp"}}, "dyn_measure");
	//Initialize vertex list to reduce warm up time
	qmc.trigger_event("initial build");
}

mc::~mc()
{}

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
	config.serialize(d);
	d.close();
	seed_write(dir+"seed");
	std::ofstream f(dir+"bins");
	if (is_thermalized())
	{
		f << "Thermalization: Done." << std::endl
			<< "Sweeps: " << (sweep - n_warmup) << std::endl
			<< "Static bins: " << static_cast<int>(measure_static_cnt / n_prebin)
			<< std::endl
			<< "Dynamic bins: " << static_cast<int>(measure_dyn_cnt / n_prebin)
			<< std::endl;
	}
	else
	{
		f << "Thermalization: " << sweep << std::endl
			<< "Sweeps: 0" << std::endl
			<< "Static bins: 0" << std::endl
			<< "Dynamic bins: 0" << std::endl;
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
		config.serialize(d);
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
	for (int n = 0; n < config.M.max_order(); ++n)
	{
		qmc.do_update(config.measure);
		if (is_thermalized())
		{
			++measure_cnt;
			if (n_cycles == measure_cnt)
			{
				++measure_static_cnt;
				qmc.do_measurement();
				measure_cnt = 0;
			}
		}
		config.M.advance_backward();
		config.M.stabilize_backward();
	}
	if (is_thermalized())
	{
		++measure_dyn_cnt;
		qmc.trigger_event("dyn_measure");
	}
	for (int n = 0; n < config.M.max_order(); ++n)
	{
		config.M.advance_forward();
		qmc.do_update(config.measure);
		config.M.stabilize_forward();
		if (is_thermalized())
		{
			++measure_cnt;
			if (n_cycles == measure_cnt)
			{
				++measure_static_cnt;
				qmc.do_measurement();
				measure_cnt = 0;
			}
		}
	}
	if (is_thermalized())
	{
		++measure_dyn_cnt;
		qmc.trigger_event("dyn_measure");
	}
	++sweep;
	if (!is_thermalized())
		qmc.trigger_event("max_order");
	if (sweep == n_warmup)
	{
		qmc.set_proposal_rates({static_cast<double>(config.M.non_ident(0)),
			static_cast<double>(std::max(5, config.M.non_ident(1)))});
		std::cout << "Max order set to " << config.M.max_order() << "."
			<< std::endl;
	}
	status();
}

void mc::do_measurement()
{}

void mc::status()
{
//	if (sweep % 500 == 0)
//		std::cout << "sweep: " << sweep << ", static_meas=" << measure_static_cnt
//			<< ", dyn_meas=" << measure_dyn_cnt << std::endl;
}
