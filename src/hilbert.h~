#pragma once
#include <iostream>
#include <map>
#include <functional>
#include <armadillo>
#include "lattice.h"
#include "honeycomb.h"

typedef long long int int_t;

struct state
{
	int_t sign;
	int_t id;
};

class hilbert
{
	public:
		hilbert(lattice& lat_)
			: lat(lat_)
		{
			dim = std::pow(2, lat.n_sites());
		}

		int_t n_el(const state& psi)
		{
			int_t n = 0;
			for (int_t i = 0; i < lat.n_sites(); ++i)
				n += n_i(psi, i);
			return n;
		}
	
		int_t n_i(const state& psi, int_t i)
		{
			return psi.sign * test_bit(psi.id, i);
		}

		state c_i(const state& psi, int_t i)
		{
			if (!test_bit(psi.id, i))
				return state{0, 0};
			int_t sign = psi.sign;
			for (int_t j=i+1; j < base_index.size(); ++j)
				if (test_bit(psi.id, j))
					sign *= -1;
			return state{sign, clear_bit(psi.id, i)};
		}

		state c_dag_i(const state& psi, int_t i)
		{
			if (test_bit(psi.id, i))
				return state{0, 0};
			int_t sign = psi.sign;
			for (int_t j=i+1; j < base_index.size(); ++j)
				if (test_bit(psi.id, j))
					sign *= -1;
			return state{sign, set_bit(psi.id, i)};
		}

		int_t sub_dimension() { return base_index.size(); }
		int_t dimension() { return dim; }
		int_t index(int_t state_id) { return base_index[state_id]; }

		void build_basis(std::function<bool(int_t)> keep)
		{
			base_index.clear();
			int_t cnt = 0;
			for (int_t i = 0; i < dim; ++i)
			{
				if (keep(i))
				{
					base_index[i] = cnt;
					++cnt;
				}
			}
		}

		void build_operator(std::function<void(const std::pair<int_t, int_t>&)> op)
		{
			std::for_each(base_index.begin(), base_index.end(), op);
		}
	private:
		int_t set_bit(int_t integer, int_t offset)
		{
			return integer | (1 << offset);
		}
		int_t clear_bit(int_t integer, int_t offset)
		{
			return integer & (~(1 << offset));
		}
		int_t invert_bit(int_t integer, int_t offset)
		{
			return integer ^ (1 << offset);
		}
		int_t test_bit(int_t integer, int_t offset)
		{
			return (integer & (1 << offset)) >> offset;
		}
	private:
		lattice& lat;
		std::map<int_t, int_t> base_index;
		int_t dim;
};
