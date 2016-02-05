#pragma once
#include <vector>
#include <algorithm>
#include "measurements.h"
#include "fast_update.h"
#include "lattice.h"
#include "greens_function.h"
#include "configuration.h"
#include "Random.h"

// k! / (k+n)!
template<typename T>
double factorial_ratio(T k, T n)
{
	if (k <= T(0) || n <= T(0))
		return 1.0;
	double result = 1.0;
	for (int i = 1; i <= n; ++i)
		result /= static_cast<double>(k + i);
	return result;
}

// ------------ QMC move : inserting a vertex ------------------

template<int N>
struct move_insert
{
	configuration* config;
	Random& rng;

	double attempt()
	{
		std::vector<arg_t> vec; vec.reserve(2*N);
		for (int i = 0; i < N; ++i)
		{
			double tau = config->params.beta * rng();
			int s1 = config->l.n_sites() * rng();
			int s2 = config->l.neighbors(s1, "nearest neighbors")
				[config->l.neighbors(s1, "nearest neighbors").size() * rng()];
			vec.push_back(arg_t{tau, s1, false});
			vec.push_back(arg_t{tau, s2, false});
		}
		int k = config->perturbation_order();
		double det_ratio = config->M.try_add<N>(vec, nn_int);
		assert(det_ratio == det_ratio && "nan value in det ratio");
		return config->params.rem[N-1] / config->params.add[N-1]
			* std::pow(-config->params.beta * config->params.V *
			config->l.n_bonds(), N) * factorial_ratio(k, N) * det_ratio;
	}

	double accept()
	{
		config->M.finish_add<N>();
		config->measure.add("insertion n="+std::to_string(N), 1.0);
		return 1.0;
	}

	void reject()
	{
		config->measure.add("insertion n="+std::to_string(N), 0.0);
	}
};

// ------------ QMC move : deleting a vertex ------------------

template<int N>
struct move_remove
{
	configuration* config;
	Random& rng;

	double attempt()
	{
		int k = config->perturbation_order();
		if (k < N) return 0.0;
		std::vector<int> vec; vec.reserve(N);
		for (int i = 0; i < N; ++i)
		{
			int p = k * rng();
			while (std::find(vec.begin(), vec.end(), p) != vec.end())
				p = k * rng();
			vec.push_back(p);
		}
		std::sort(vec.begin(), vec.end());
		double det_ratio = config->M.try_remove<N>(vec, nn_int);
		assert(det_ratio == det_ratio && "nan value in det ratio");
		return config->params.add[N-1] / config->params.rem[N-1]
			* std::pow(-config->params.beta * config->params.V *
			config->l.n_bonds(), -N) / factorial_ratio(k-N, N) * det_ratio;
	}

	double accept()
	{
		config->M.finish_remove<N>();
		config->measure.add("removal n="+std::to_string(N), 1.0);
		return 1.0;
	}

	void reject()
	{
		config->measure.add("removal n="+std::to_string(N), 0.0);
	}
};

// ------------ QMC move : Z -> W2 ------------------

struct move_ZtoW2
{
	configuration* config;
	Random& rng;
	bool save_acc;

	double attempt()
	{
		if (config->worms() != 0)
		{
			save_acc = false;
			return 0.0;
		}
		double tau = config->params.beta * rng();
		int s1 = config->l.n_sites() * rng();
		const std::vector<int>& neighbors =
			config->l.neighbors(s1, "worm nhood");
		int s2 = neighbors[neighbors.size() * rng()];
		std::vector<arg_t> vec = {arg_t{tau, s1, true}, arg_t{tau, s2, true}};
		double det_ratio = config->M.try_add<1>(vec, worm);
		assert(det_ratio == det_ratio && "nan value in det ratio");
		save_acc = true; 
		return config->params.W2toZ / config->params.ZtoW2
			* config->l.parity(s1) * config->l.parity(s2) * config->params.zeta2
			* config->params.ratio_w2 * det_ratio;
	}

	double accept()
	{
		config->M.finish_add<1>();
		if (save_acc)
			config->measure.add("Z -> W2", 1.0);
		return 1.0;
	}

	void reject()
	{
		if (save_acc)
			config->measure.add("Z -> W2", 0.0);
	}
};

// ------------ QMC move : W2 -> Z ------------------

struct move_W2toZ
{
	configuration* config;
	Random& rng;
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
			config->l.neighbors(sites[0], "worm nhood");
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
		return config->params.ZtoW2 / config->params.W2toZ
			* config->l.parity(sites[0]) * config->l.parity(sites[1])
			/ config->params.zeta2 / config->params.ratio_w2 * det_ratio;
	}

	double accept()
	{
		config->M.finish_remove<1>();
		if (save_acc)
			config->measure.add("W2 -> Z", 1.0);
		return 1.0;
	}

	void reject()
	{
		if (save_acc)
			config->measure.add("W2 -> Z", 0.0);
	}
};

// ------------ QMC move : W2 -> W4 ------------------

struct move_W2toW4
{
	configuration* config;
	Random& rng;
	bool save_acc;

	double attempt()
	{
		if (config->worms() != 1)
		{
			save_acc = false;
			return 0.0;
		}
		double tau = config->M.vertex(0, worm).tau;
		int s1 = config->l.n_sites() * rng();
		const std::vector<int>& neighbors =
			config->l.neighbors(s1, "worm nhood");
		int s2 = neighbors[neighbors.size() * rng()];
		std::vector<arg_t> vec = {arg_t{tau, s1, true}, arg_t{tau, s2, true}};
		double det_ratio = config->M.try_add<1>(vec, worm);
		assert(det_ratio == det_ratio && "nan value in det ratio");
		save_acc = true;
		return config->params.W4toW2 / config->params.W2toW4
			* config->l.parity(s1) * config->l.parity(s2) * config->params.zeta4
			/ config->params.zeta2 * config->params.ratio_w2 * det_ratio;
	}

	double accept()
	{
		config->M.finish_add<1>();
		if (save_acc)
			config->measure.add("W2 -> W4", 1.0);
		return 1.0;
	}

	void reject()
	{
		if (save_acc)
			config->measure.add("W2 -> W4", 0.0);
	}
};

// ------------ QMC move : W4 -> W2 ------------------

struct move_W4toW2
{
	configuration* config;
	Random& rng;
	bool save_acc;

	double attempt()
	{
		if (config->worms()	!= 2)
		{
			save_acc = false;
			return 0.0;
		}
		int p = config->worms() * rng();	
		int sites[] = {config->M.vertex(2*p, worm).site,
							config->M.vertex(2*p+1, worm).site};
		const std::vector<int>& neighbors =
			config->l.neighbors(sites[0], "worm nhood");
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
		return config->params.W2toW4 / config->params.W4toW2
			* config->l.parity(sites[0]) * config->l.parity(sites[1])
			* config->params.zeta2 / config->params.zeta4
			/ config->params.ratio_w2 * det_ratio;
	}

	double accept()
	{
		config->M.finish_remove<1>();
		if (save_acc)
			config->measure.add("W4 -> W2", 1.0);
		return 1.0;
	}

	void reject()
	{
		if (save_acc)
			config->measure.add("W4 -> W2", 0.0);
	}
};

// ------------ QMC move : Z -> W4 ------------------

struct move_ZtoW4
{
	configuration* config;
	Random& rng;
	bool save_acc;

	double attempt()
	{
		if (config->worms() != 0)
		{
			save_acc = false;
			return 0.0;
		}
		double tau = config->params.beta * rng();
		int s1 = config->l.n_sites() * rng();
		const std::vector<int>& neighbors =
			config->l.neighbors(s1, "worm nhood");
		int s2 = neighbors[neighbors.size() * rng()];
		int s3 = neighbors[neighbors.size() * rng()];
		int s4 = neighbors[neighbors.size() * rng()];
		std::vector<arg_t> vec = {arg_t{tau, s1, true}, arg_t{tau, s2, true},
			arg_t{tau, s3, true}, arg_t{tau, s4, true}};
		double det_ratio = config->M.try_add<2>(vec, worm);
		assert(det_ratio == det_ratio && "nan value in det ratio");
		save_acc = true;
		return config->params.W4toZ / config->params.ZtoW4 * config->l.parity(s1)
			* config->l.parity(s2) * config->l.parity(s3) * config->l.parity(s4)
			* config->params.zeta4 * config->params.ratio_w4 * det_ratio;
	}

	double accept()
	{
		config->M.finish_add<2>();
		if (save_acc)
			config->measure.add("Z -> W4", 1.0);
		return 1.0;
	}

	void reject()
	{	
		if (save_acc)
			config->measure.add("Z -> W4", 0.0);
	}
};

// ------------ QMC move : W4 -> Z ------------------

struct move_W4toZ
{
	configuration* config;
	Random& rng;
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
			config->l.neighbors(sites[0], "worm nhood");
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
		return config->params.ZtoW4 / config->params.W4toZ
			* config->l.parity(config->M.vertex(p, worm).site)
			* config->l.parity(config->M.vertex(p+1, worm).site)
			* config->l.parity(config->M.vertex(p+2, worm).site)
			* config->l.parity(config->M.vertex(p+3, worm).site)
			/ config->params.zeta4 / config->params.ratio_w4 * det_ratio;
	}

	double accept()
	{
		config->M.finish_remove<2>();
		if (save_acc)
			config->measure.add("W4 -> Z", 1.0);
		return 1.0;
	}

	void reject()
	{
		if (save_acc)
			config->measure.add("W4 -> Z", 0.0);
	}
};

// ------------ QMC move : worm shift ------------------

struct move_shift
{
	configuration* config;
	Random& rng;
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
			* rng();
		if (worm_vert[0].tau + tau_shift > config->params.beta)
			tau_shift -= config->params.beta;
		else if(worm_vert[0].tau + tau_shift < 0.0)
			tau_shift += config->params.beta;
		for (int i = 0; i < worm_vert.size(); ++i)
			worm_vert[i].tau += tau_shift;
		int p = worm_vert.size() * rng();
		const std::vector<int>& neighbors = config->l.neighbors(
			worm_vert[p].site, "shift nhood");
		int old_site = worm_vert[p].site;
		worm_vert[p].site = neighbors[neighbors.size() * rng()];
		int new_site = worm_vert[p].site;
		std::random_shuffle(worm_vert.begin(), worm_vert.end());
		double det_ratio;
		if (config->worms() == 1)
			det_ratio = config->M.try_shift<1>(worm_vert);
		else if (config->worms() == 2)
			det_ratio = config->M.try_shift<2>(worm_vert);
		assert(det_ratio == det_ratio && "nan value in det ratio");
		save_acc = true;
		return config->l.parity(old_site) * config->l.parity(new_site)
			* det_ratio;
	}

	double accept()
	{
		if (config->worms() == 1)
			config->M.finish_shift<1>();
		else if (config->worms() == 2)
			config->M.finish_shift<2>();
		if (save_acc)
			config->measure.add("worm shift", 1.0);
		return 1.0;
	}

	void reject()
	{
		if (save_acc)
			config->measure.add("worm shift", 0.0);
	}

};
// ------------ QMC move : worm shift ------------------

struct move_shift_2
{
	configuration* config;
	Random& rng;
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
			* rng();
		if (worm_vert[0].tau + tau_shift > config->params.beta)
			tau_shift -= config->params.beta;
		else if(worm_vert[0].tau + tau_shift < 0.0)
			tau_shift += config->params.beta;
		for (int i = 0; i < worm_vert.size(); ++i)
			worm_vert[i].tau += tau_shift;
		int p = worm_vert.size() * rng();
		const std::vector<int>& neighbors = config->l.neighbors(
			worm_vert[p].site, "shift nhood");
		int old_site = worm_vert[p].site;
		worm_vert[p].site = neighbors[neighbors.size() * rng()];
		int new_site = worm_vert[p].site;
		std::random_shuffle(worm_vert.begin(), worm_vert.end());
		double det_ratio;
		if (config->worms() == 1)
			det_ratio = config->M.try_shift<1>(worm_vert);
		else if (config->worms() == 2)
			det_ratio = config->M.try_shift<2>(worm_vert);
		assert(det_ratio == det_ratio && "nan value in det ratio");
		save_acc = true;
		return config->l.parity(old_site) * config->l.parity(new_site)
			* det_ratio;
	}

	double accept()
	{
		if (config->worms() == 1)
			config->M.finish_shift<1>();
		else if (config->worms() == 2)
			config->M.finish_shift<2>();
		if (save_acc)
			config->measure.add("worm shift 2", 1.0);
		return 1.0;
	}

	void reject()
	{
		if (save_acc)
			config->measure.add("worm shift 2", 0.0);
	}
};
