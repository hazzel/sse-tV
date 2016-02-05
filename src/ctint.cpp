#include <iostream>
#include <cassert>
#include <cmath>
#include <initializer_list>
#include <triqs/mc_tools.hpp>
#include <triqs/det_manip.hpp>
#include <triqs/statistics.hpp>
#include <boost/mpi.hpp>
#include <boost/serialization/complex.hpp>
#include "ctint.hpp"
#include "greens_function.h"
#include "fast_update.h"

#define eigen_impl

// --------------- The QMC configuration ----------------

// Argument type
struct arg_t
{
	double tau;
	int site;
	bool worm;
};

enum { nn_int, worm };

// The function that appears in the calculation of the determinant
struct full_g_entry
{
	const greens_function& g0;

	double operator()(const arg_t& x, const arg_t& y) const
	{
		if ((x.tau == y.tau) && (x.site == y.site))
			return 0.0;
		else
			return g0(x.tau - y.tau, x.site, y.site);
	}
};

struct parameters
{	
	double beta, V, zeta2, zeta4;
	int worm_nhood_dist;
	double ratio_w2, ratio_w4;
};

// The Monte Carlo configuration
struct configuration
{
	const lattice& l;
	triqs::det_manip::det_manip<full_g_entry> Mmatrix;
	fast_update<full_g_entry, arg_t> M;
	parameters params;
	triqs::statistics::observable<double> obs_pert_order;
	triqs::statistics::observable<double> obs_Z;
	triqs::statistics::observable<double> obs_W2;
	triqs::statistics::observable<double> obs_W4;
	
	triqs::statistics::observable<double> acc_ZtoW2;
	triqs::statistics::observable<double> acc_W2toZ;
	triqs::statistics::observable<double> acc_ZtoW4;
	triqs::statistics::observable<double> acc_W4toZ;
	triqs::statistics::observable<double> acc_W2toW4;
	triqs::statistics::observable<double> acc_W4toW2;
	triqs::statistics::observable<double> acc_shift;

	#ifndef eigen_impl
	int perturbation_order() const { return Mmatrix.size() / 2; }
	#else
	int perturbation_order() const { return M.perturbation_order(nn_int); }
	int worms() const { return M.perturbation_order(worm); }
	#endif

	configuration(const lattice& l_, const greens_function& g0, 
		const parameters& params_)
		: l(l_), Mmatrix{full_g_entry{g0}, 100}, M{full_g_entry{g0}, l_, 2},
			params(params_), obs_pert_order()
	{}
};

// ------------ QMC move : inserting a vertex ------------------

struct move_insert
{
	configuration* config;
	triqs::mc_tools::random_generator& rng;

	double attempt()
	{
		double tau = rng(config->params.beta);
		int s1 = rng(config->l.n_sites());
		int s2 = config->l.neighbors(s1, 1)
			[rng(config->l.neighbors(s1, 1).size())];
		int k = config->perturbation_order();
		#ifndef eigen_impl
		double det_ratio = config->Mmatrix.try_insert2(2*k, 2*k+1, 2*k,
			2*k+1, {tau, s1}, {tau, s2}, {tau, s1}, {tau, s2});
		#else
		std::vector<arg_t> vec = {arg_t{tau, s1, false}, arg_t{tau, s2, false}};
		double det_ratio = config->M.try_add<1>(vec, nn_int);
		#endif
		assert(det_ratio == det_ratio && "nan value in det ratio");
		return -config->params.beta * config->params.V * config->l.n_bonds()
			/ (k + 1) * det_ratio;
	}

	double accept()
	{
		#ifndef eigen_impl
		config->Mmatrix.complete_operation(); // Finish insertion
		#else
		config->M.finish_add();
		#endif
		return 1.0;
	}

	void reject() {}
};

// ------------ QMC move : deleting a vertex ------------------

struct move_remove
{
	configuration* config;
	triqs::mc_tools::random_generator& rng;

	double attempt()
	{
		int k = config->perturbation_order();
		if (k <= 0) return 0;
		int p = rng(k); // Choose one of the operators for removal
		#ifndef eigen_impl
		double det_ratio = config->Mmatrix.try_remove2(2*p, 2*p+1, 2*p, 2*p+1);
		#else
		std::vector<int> vec = {p};
		double det_ratio = config->M.try_remove<1>(vec, nn_int);
		#endif
		assert(det_ratio == det_ratio && "nan value in det ratio");
		return -k / (config->params.beta * config->params.V
			* config->l.n_bonds()) * det_ratio;
	}

	double accept()
	{
		#ifndef eigen_impl
		config->Mmatrix.complete_operation(); // Finish removal
		#else
		config->M.finish_remove();
		#endif
		return 1.0;
	}

	void reject() {}
};

// ------------ QMC move : Z -> W2 ------------------

struct move_ZtoW2
{
	configuration* config;
	triqs::mc_tools::random_generator& rng;
	bool save_acc;

	double attempt()
	{
		if (config->worms() != 0)
		{
			save_acc = false;
			return 0.0;
		}
		double tau = rng(config->params.beta);
		int s1 = rng(config->l.n_sites());
		const std::vector<int>& neighbors =
			config->l.neighbors(s1, config->params.worm_nhood_dist);
		int s2 = neighbors[rng(neighbors.size())];
		std::vector<arg_t> vec = {arg_t{tau, s1, true}, arg_t{tau, s2, true}};
		double det_ratio = config->M.try_add<1>(vec, worm);
		assert(det_ratio == det_ratio && "nan value in det ratio");
		save_acc = true;
		return config->l.parity(s1) * config->l.parity(s2)
			* config->params.zeta2 * config->params.ratio_w2 * det_ratio;
	}

	double accept()
	{
		config->M.finish_add();
		if (save_acc)
			config->acc_ZtoW2 << 1.0;
		return 1.0;
	}

	void reject()
	{
		if (save_acc)
			config->acc_ZtoW2 << 0.0;
	}
};

// ------------ QMC move : W2 -> Z ------------------

struct move_W2toZ
{
	configuration* config;
	triqs::mc_tools::random_generator& rng;
	bool save_acc;

	double attempt()
	{
		if (config->worms()	!= 1)
		{
			save_acc = false;
			return 0.0;
		}
		int p = 0; //only one worm exists	
		int sites[] = {config->M.vertex(p, worm).site,
							config->M.vertex(p+1, worm).site};
		const std::vector<int>& neighbors =
			config->l.neighbors(sites[0], config->params.worm_nhood_dist);
		if (std::find(neighbors.begin(), neighbors.end(), sites[1])
			== neighbors.end())
		{
			save_acc = false;
			return 0.0;
		}
		std::vector<int> vec = {p};
		double det_ratio = config->M.try_remove<1>(vec, worm);
		assert(det_ratio == det_ratio && "nan value in det ratio");
		save_acc = true;
		return config->l.parity(sites[0]) * config->l.parity(sites[1])
			/ config->params.zeta2 / config->params.ratio_w2 * det_ratio;
	}

	double accept()
	{
		config->M.finish_remove();
		if (save_acc)
			config->acc_W2toZ << 1.0;
		return 1.0;
	}

	void reject()
	{
		if (save_acc)
			config->acc_W2toZ << 0.0;
	}
};

// ------------ QMC move : W2 -> W4 ------------------

struct move_W2toW4
{
	configuration* config;
	triqs::mc_tools::random_generator& rng;
	bool save_acc;

	double attempt()
	{
		if (config->worms() != 1)
		{
			save_acc = false;
			return 0.0;
		}
		double tau = config->M.vertex(0, worm).tau;
		int s1 = rng(config->l.n_sites());
		const std::vector<int>& neighbors =
			config->l.neighbors(s1, config->params.worm_nhood_dist);
		int s2 = neighbors[rng(neighbors.size())];
		std::vector<arg_t> vec = {arg_t{tau, s1, true}, arg_t{tau, s2, true}};
		double det_ratio = config->M.try_add<1>(vec, worm);
		assert(det_ratio == det_ratio && "nan value in det ratio");
		save_acc = true;
		return config->l.parity(s1) * config->l.parity(s2)
			* config->params.zeta4 / config->params.zeta2
			* config->params.ratio_w2 * det_ratio;
	}

	double accept()
	{
		config->M.finish_add();
		if (save_acc)
			config->acc_W2toW4 << 1.0;
		return 1.0;
	}

	void reject()
	{
		if (save_acc)
			config->acc_W2toW4 << 0.0;
	}
};

// ------------ QMC move : W4 -> W2 ------------------

struct move_W4toW2
{
	configuration* config;
	triqs::mc_tools::random_generator& rng;
	bool save_acc;

	double attempt()
	{
		if (config->worms()	!= 2)
		{
			save_acc = false;
			return 0.0;
		}
		int p = rng(config->worms()); //only one worm exists	
		int sites[] = {config->M.vertex(2*p, worm).site,
							config->M.vertex(2*p+1, worm).site};
		const std::vector<int>& neighbors =
			config->l.neighbors(sites[0], config->params.worm_nhood_dist);
		if (std::find(neighbors.begin(), neighbors.end(), sites[1])
			== neighbors.end())
		{
			save_acc = false;
			return 0.0;
		}
		std::vector<int> vec = {p};
		double det_ratio = config->M.try_remove<1>(vec, worm);
		assert(det_ratio == det_ratio && "nan value in det ratio");
		save_acc = true;
		return config->l.parity(sites[0]) * config->l.parity(sites[1])
			* config->params.zeta2 / config->params.zeta4
			/ config->params.ratio_w2 * det_ratio;
	}

	double accept()
	{
		config->M.finish_remove();
		if (save_acc)
			config->acc_W4toW2 << 1.0;
		return 1.0;
	}

	void reject()
	{
		if (save_acc)
			config->acc_W4toW2 << 0.0;
	}
};

// ------------ QMC move : Z -> W4 ------------------

struct move_ZtoW4
{
	configuration* config;
	triqs::mc_tools::random_generator& rng;
	bool save_acc;

	double attempt()
	{
		if (config->worms() != 0)
		{
			save_acc = false;
			return 0.0;
		}
		double tau = rng(config->params.beta);
		int s1 = rng(config->l.n_sites());
		const std::vector<int>& neighbors =
			config->l.neighbors(s1, config->params.worm_nhood_dist);
		int s2 = neighbors[rng(neighbors.size())];
		int s3 = neighbors[rng(neighbors.size())];
		int s4 = neighbors[rng(neighbors.size())];
		std::vector<arg_t> vec = {arg_t{tau, s1, true}, arg_t{tau, s2, true},
			arg_t{tau, s3, true}, arg_t{tau, s4, true}};
		double det_ratio = config->M.try_add<2>(vec, worm);
		assert(det_ratio == det_ratio && "nan value in det ratio");
		save_acc = true;
		return config->l.parity(s1) * config->l.parity(s2) * config->l.parity(s3)
			* config->l.parity(s4) * config->params.zeta4
			* config->params.ratio_w4 * det_ratio;
	}

	double accept()
	{
		config->M.finish_add();
		if (save_acc)
			config->acc_ZtoW4 << 1.0;
		return 1.0;
	}

	void reject()
	{	
		if (save_acc)
			config->acc_ZtoW4 << 0.0;
	}
};

// ------------ QMC move : W4 -> Z ------------------

struct move_W4toZ
{
	configuration* config;
	triqs::mc_tools::random_generator& rng;
	bool save_acc;

	double attempt()
	{
		if (config->worms() != 2)
		{
			save_acc = false;
			return 0.0;
		}
		int p = 0;
		int sites[] = {config->M.vertex(p, worm).site,
							config->M.vertex(p+1, worm).site,
							config->M.vertex(p+2, worm).site,
							config->M.vertex(p+3, worm).site};
		const std::vector<int>& neighbors =
			config->l.neighbors(sites[0], config->params.worm_nhood_dist);
		for (int i = 1; i < 4; ++i)	
			if (std::find(neighbors.begin(), neighbors.end(), sites[i])
				== neighbors.end())
			{
				save_acc = false;
				return 0.0;
			}
		std::vector<int> vec = {p, p+1};
		double det_ratio = config->M.try_remove<2>(vec, worm);
		assert(det_ratio == det_ratio && "nan value in det ratio");
		save_acc = true;
		return config->l.parity(config->M.vertex(p, worm).site)
			* config->l.parity(config->M.vertex(p+1, worm).site)
			* config->l.parity(config->M.vertex(p+2, worm).site)
			* config->l.parity(config->M.vertex(p+3, worm).site)
			/ config->params.zeta4 / config->params.ratio_w4 * det_ratio;
	}

	double accept()
	{
		config->M.finish_remove();
		if (save_acc)
			config->acc_W4toZ << 1.0;
		return 1.0;
	}

	void reject()
	{
		if (save_acc)
			config->acc_W4toZ << 0.0;
	}
};

// ------------ QMC move : worm shift ------------------

struct move_shift
{
	configuration* config;
	triqs::mc_tools::random_generator& rng;
	bool save_acc;

	double attempt()
	{
		if (config->worms() == 0)
		{
			save_acc = false;
			return 0.0;
		}
		std::vector<arg_t> worm_vert(2*config->worms());
		for (int i = 0; i < worm_vert.size(); ++i)
			worm_vert[i] = config->M.vertex(i, worm);
		double tau_shift = -0.05*config->params.beta + 0.1*config->params.beta
			* rng(1.0);
		if (worm_vert[0].tau + tau_shift > config->params.beta)
			tau_shift -= config->params.beta;
		else if(worm_vert[0].tau + tau_shift < 0.0)
			tau_shift += config->params.beta;
		for (int i = 0; i < worm_vert.size(); ++i)
			worm_vert[i].tau += tau_shift;
		int p = rng(worm_vert.size());
		const std::vector<int>& neighbors = config->l.neighbors(
			worm_vert[p].site, config->params.worm_nhood_dist);
		int old_site = worm_vert[p].site;
		worm_vert[p].site = neighbors[rng(neighbors.size())];
		double det_ratio = config->M.try_shift(worm_vert);
		assert(det_ratio == det_ratio && "nan value in det ratio");
		save_acc = true;
		//std::cout << "worm shift try:" << std::endl;
		//std::cout << "shift buffer: " << tau_shift << std::endl;
		//std::cout << worm_vert[0].tau << " , " << worm_vert[0].site << std::endl;
		//std::cout << worm_vert[1].tau << " , " << worm_vert[1].site << std::endl;
		//config->M.print_vertices();
		return config->l.parity(old_site) * config->l.parity(worm_vert[p].site)
			* det_ratio;
	}

	double accept()
	{
		config->M.finish_shift();
		//std::cout << "worm shift done:" << std::endl;
		//config->M.print_vertices();
		//std::cout << std::endl;
		//std::cin.get();
		if (save_acc)
			config->acc_shift << 1.0;
		return 1.0;
	}

	void reject()
	{
		if (save_acc)
			config->acc_shift << 0.0;
	}
};

//  -------------- QMC measurement ----------------

struct measure_M
{
	configuration* config;
	double Z = 0, k = 0;

	measure_M(configuration* config_)
		: config(config_)
	{}

	void accumulate(double sign)
	{
		Z += sign;
		if (config->worms() == 0) //measure Z
		{
			k += config->perturbation_order();
			config->obs_pert_order << config->perturbation_order();
			config->obs_Z << 1.0;
			config->obs_W2 << 0.0;
			config->obs_W4 << 0.0;
		}
		else if (config->worms() == 1) //measure W2
		{
			config->obs_Z << 0.0;
			config->obs_W2 << 1.0;
			config->obs_W4 << 0.0;
		}
		else if (config->worms() == 2) //measure W4
		{
			config->obs_Z << 0.0;
			config->obs_W2 << 0.0;
			config->obs_W4 << 1.0;
		}
	};

	void collect_results(const boost::mpi::communicator& c)
	{
		boost::mpi::all_reduce(c, Z, Z, std14::plus<>());
		boost::mpi::all_reduce(c, k, k, std14::plus<>());
		//boost::mpi::all_reduce(c, config->obs_pert_order, config->obs_pert_order,
		//	std14::plus<triqs::statistics::observable<double>>());
		
		int bin_size = 1000;
//		triqs::h5::file file("test_store.h5", H5F_ACC_RDWR);
//		triqs::h5::group root(file);
//		triqs::h5::group run = root.create_group("run "
//			+ std::to_string(c.rank()));
//		triqs::h5::group gr_measure = run.create_group("observables");
//		triqs::h5::h5_write(gr_measure, "pert_order", triqs::statistics::
//			make_binned_series(config->obs_pert_order, bin_size).data());

		if (c.rank() == 0)
		{
			try
			{	
				std::cout << "Average perturbation order = "
					<< triqs::statistics::average_and_error(config->obs_pert_order,
							bin_size)
					//<< " (Autocorrelationtime = "
					//<< triqs::statistics::autocorrelation_time_from_binning(
					//	config->obs_pert_order)
					<< std::endl;
				std::cout << "Z = "
					<< triqs::statistics::average_and_error(config->obs_Z, bin_size)
					//<< " (Autocorrelationtime = "
					//<< triqs::statistics::autocorrelation_time_from_binning(
					//	config->obs_Z)
					<< std::endl;
				std::cout << "W2 = "
					<< triqs::statistics::average_and_error(config->obs_W2, bin_size)
					//<< " (Autocorrelationtime = "
					//<< triqs::statistics::autocorrelation_time_from_binning(
					//	config->obs_W2)
					<< std::endl;
				std::cout << "W4 = "
					<< triqs::statistics::average_and_error(config->obs_W4, bin_size)
					//<< " (Autocorrelationtime = "
					//<< triqs::statistics::autocorrelation_time_from_binning(
					//	config->obs_W4)
					<< std::endl;
				std::cout << "M_2 = "
					<< triqs::statistics::average_and_error(config->obs_W2
						/ config->obs_Z / config->params.zeta2)
					<< std::endl;
				std::cout << "M_4 = "
					<< triqs::statistics::average_and_error(config->obs_W4
						/ config->obs_Z / config->params.zeta4)
					<< std::endl;
				std::cout << "Binder = "
					<< triqs::statistics::average_and_error(config->obs_W4
						* config->obs_Z / (config->obs_W2 * config->obs_W2)
						* config->params.zeta2 * config->params.zeta2
						/ config->params.zeta4)
					<< std::endl;
			
				std::cout << "ZtoW2 = "
					<< triqs::statistics::average_and_error(config->acc_ZtoW2)
					<< std::endl;
				std::cout << "W2toZ = "
					<< triqs::statistics::average_and_error(config->acc_W2toZ)
					<< std::endl;
				std::cout << "ZtoW4 = "
					<< triqs::statistics::average_and_error(config->acc_ZtoW4)
					<< std::endl;
				std::cout << "W4toZ = "
					<< triqs::statistics::average_and_error(config->acc_W4toZ)
					<< std::endl;
				std::cout << "W2toW4 = "
					<< triqs::statistics::average_and_error(config->acc_W2toW4)
					<< std::endl;
				std::cout << "W4toW2 = "
					<< triqs::statistics::average_and_error(config->acc_W4toW2)
					<< std::endl;
				std::cout << "worm shift = "
					<< triqs::statistics::average_and_error(config->acc_shift)
					<< std::endl;
			}
			catch(const triqs::runtime_error& c)
			{
				std::cout << c.what() << std::endl;
			}
		}
	}
};

// ------------ The main class of the solver ------------------------

ctint_solver::ctint_solver(long unsigned int base_seed_, double beta_,
	int n_slices_)
	: base_seed(base_seed_), beta(beta_), n_slices(n_slices_)
{}

void ctint_solver::write_parameters(triqs::h5::group& group,
	configuration& config)
{
	triqs::h5::h5_write(group, "seed", std::to_string(seed));
	triqs::h5::h5_write(group, "L", std::sqrt(config.l.n_sites()/2.0));
	triqs::h5::h5_write(group, "beta", config.params.beta);
	triqs::h5::h5_write(group, "V", config.params.V);
	triqs::h5::h5_write(group, "zeta2", config.params.zeta2);
	triqs::h5::h5_write(group, "zeta4", config.params.zeta4);
	triqs::h5::h5_write(group, "worm neighborhood distance",
	config.params.worm_nhood_dist);
}

// The method that runs the qmc
void ctint_solver::solve(int L, double V, int n_cycles, int length_cycle,
	int n_warmup_cycles, std::string random_name, int max_time)
{
	boost::mpi::communicator world;

	// Rank-specific variables
	int verbosity = (world.rank() == 0 ? 3 : 0);
	seed = base_seed + 928374 * world.rank();
	if (world.rank() == 0)
	{
		std::cout << "CT-INT with L=" << L << ", V=" << V << ", beta=" << beta
			<< std::endl;
		#ifndef eigen_impl
		std::cout << "Using 'triqs' for fast updates." << std::endl;
		#else
		std::cout << "Using 'eigen' for fast updates." << std::endl;
		#endif
	}
	std::cout << "Rank " << world.rank() << " using seed: " << seed << std::endl;

	// Construct a Monte Carlo loop
	triqs::mc_tools::mc_generic<double> CTQMC(n_cycles, length_cycle,
		n_warmup_cycles, random_name, seed, verbosity);

	// Prepare the configuration
	honeycomb h(L);
	lattice l;
	l.generate_graph(h);
	if (worm_nhood_dist == -1)
		worm_nhood_dist = l.max_distance();
	l.generate_neighbor_map(worm_nhood_dist);
	double ratio_w2 = static_cast<double>(l.neighbors(0, worm_nhood_dist).size())
	  	/ static_cast<double>(l.n_sites());
	double ratio_w4 = std::pow(static_cast<double>(l.neighbors(
		0, worm_nhood_dist).size()) / static_cast<double>(l.n_sites()), 3.0);
	greens_function g0;
	g0.generate_mesh(&l, beta, n_slices);
	auto config = configuration{l, g0,
		{beta, V, zeta2, zeta4, worm_nhood_dist, ratio_w2, ratio_w4}};

	triqs::h5::file file("test_store.h5", H5F_ACC_TRUNC);
	triqs::h5::group root(file);
	triqs::h5::group run = root.create_group("run "
		+ std::to_string(world.rank()));
	triqs::h5::group gr_parameters = run.create_group("parameters");
	if (world.rank() == 0)
		write_parameters(gr_parameters, config);
	triqs::h5::group gr_measure = run.create_group("observables");

	// Register moves and measurements
	CTQMC.add_move(move_insert{&config, CTQMC.rng()}, "insertion");
	CTQMC.add_move(move_remove{&config, CTQMC.rng()}, "removal");
	CTQMC.add_move(move_ZtoW2{&config, CTQMC.rng(), false}, "Z -> W2");
	CTQMC.add_move(move_W2toZ{&config, CTQMC.rng(), false}, "W2 -> Z");
	CTQMC.add_move(move_ZtoW4{&config, CTQMC.rng(), false}, "Z -> W4");
	CTQMC.add_move(move_W4toZ{&config, CTQMC.rng(), false}, "W4 -> Z");
	CTQMC.add_move(move_W2toW4{&config, CTQMC.rng(), false}, "W2 -> W4");
	CTQMC.add_move(move_W4toW2{&config, CTQMC.rng(), false}, "W4 -> W2");
	CTQMC.add_move(move_shift{&config, CTQMC.rng(), false}, "worm shift");
	CTQMC.add_measure(measure_M{&config}, "M measurement");

	// Run and collect results
	CTQMC.start(1.0, triqs::utility::clock_callback(max_time));
	CTQMC.collect_results(world);
}
