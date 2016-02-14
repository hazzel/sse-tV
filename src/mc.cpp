#include <string>
#include <fstream>
#include "mc.h"
#include "move_functors.h"
#include "measure_functors.h"
#include "event_functors.h"

mc::mc(const std::string& dir)
	: rng(Random()), qmc(rng), g0{}, config(rng, g0)
{
	//Read parameters
	pars.read_file(dir);
	sweep = 0;
	measure_cnt = 0;
	n_cycles = pars.value_or_default<int>("cycles", 300);
	n_warmup = pars.value_or_default<int>("warmup", 100000);
	n_prebin = pars.value_or_default<int>("prebin", 500);
	n_rebuild = pars.value_or_default<int>("rebuild", 1000);
	n_matsubara = pars.value_or_default<int>("matsubara_freqs", 10);
	hc.L = pars.value_or_default<int>("L", 9);
	config.param.beta = 1./pars.value_or_default<double>("T", 0.2);
	config.param.t = pars.value_or_default<double>("t", 1.0);
	config.param.V = pars.value_or_default<double>("V", 1.355);
	config.param.lambda = std::log((2.*config.param.t + config.param.V)
		/ (2.*config.param.t - config.param.V));
	config.param.n_stab = pars.value_or_default<int>("stabilization", 10);

	//Proposal probabilites
	config.param.V1 = pars.value_or_default<double>("V1", 1.0);
	config.param.V2 = pars.value_or_default<double>("V2", 0.0);

	//Initialize lattice
	config.l.generate_graph(hc);
	config.l.generate_neighbor_map("nearest neighbors", [this]
		(lattice::vertex_t i, lattice::vertex_t j) {
		return config.l.distance(i, j) == 1; });

	//Initialize configuration class
	config.initialize();

	//Set up Monte Carlo moves
	qmc.add_move(move_update_vertex{config, rng, 1}, "update type 1",
		config.param.V1);
	qmc.add_move(move_update_vertex{config, rng, 2}, "update type 2",
		config.param.V2);

	//Set up measurements
	config.measure.add_observable("M2", n_prebin);
	config.measure.add_observable("<n1>", n_prebin);
	config.measure.add_vectorobservable("<n_r n_0>", config.l.max_distance() + 1,
		n_prebin);

	//Measure acceptance probabilities
	config.measure.add_observable("update type 1", n_prebin * n_cycles);
	config.measure.add_observable("sign", n_prebin * n_cycles);
	
	qmc.add_measure(measure_estimator{config, config.measure, pars,
		std::vector<double>(config.l.max_distance() + 1, 0.0)}, "measurement");
	
	//Set up events
	qmc.add_event(event_rebuild{config, config.measure}, "rebuild");
	qmc.add_event(event_build{config, rng}, "initial build");
	qmc.add_event(event_max_order{config, rng}, "max_order");
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
	for (int n = 0; n < 2*config.M.max_order(); ++n)
	{
		qmc.do_update(config.measure);
		//config.M.print_gf_from_scratch();

		if (is_thermalized())
		{
			++measure_cnt;
			if (n_cycles == measure_cnt)
			{
				qmc.do_measurement();
				measure_cnt = 0;
			}
		}
		if (n < config.M.max_order())
		{
			config.M.advance_backward();
		//	std::cout << "moved backward" << std::endl;
		}
		else
		{
			config.M.advance_forward();
		//	std::cout << "moved forward" << std::endl;
		}
		//config.M.print_bonds();
		//std::cout << "//////////////" << std::endl;
	}
	++sweep;
	if (!is_thermalized())
		qmc.trigger_event("max_order");
	//if (sweep % n_rebuild == 0)
	//	qmc.trigger_event("rebuild");
	if (sweep == n_warmup)
		std::cout << "Max order set to " << config.M.max_order() << "."
			<< std::endl;
	status();
}

void mc::do_measurement()
{}

void mc::status()
{
//	if (sweep == n_warmup)
//		std::cout << "Thermalization done." << std::endl;
//	if (is_thermalized() && sweep % (10000) == 0)
//	{
//		std::cout << "sweep: " << sweep << std::endl;
//		std::cout << "pert order: " << config.perturbation_order() << std::endl;
//	}
}
