#pragma once
#include <cmath>
#include "interpolation.h"

struct spline
{
	int n_mesh;
	double pi = 4. * std::atan(1);
	alglib::real_1d_array mesh_x;
	alglib::real_1d_array cos_mesh_y;
	alglib::real_1d_array sin_mesh_y;
	alglib::spline1dinterpolant cos_spline;
	alglib::spline1dinterpolant sin_spline;

	spline(int n_mesh_)
		: n_mesh(n_mesh_)
	{
		int n_mesh = 1000;
		mesh_x.setlength(n_mesh);
		cos_mesh_y.setlength(n_mesh);
		sin_mesh_y.setlength(n_mesh);
		for (int i = 0; i < n_mesh; ++i)
		{
			mesh_x[i] = static_cast<double>(i) / static_cast<double>(n_mesh-1)
				* 2. * 4. * std::atan(1.);
			cos_mesh_y[i] = std::cos(mesh_x[i]);
			sin_mesh_y[i] = std::sin(mesh_x[i]);
		}
		alglib::spline1dbuildakima(mesh_x, cos_mesh_y, cos_spline);
		alglib::spline1dbuildakima(mesh_x, sin_mesh_y, sin_spline);
	}

	double cos(double x) const
	{
		return alglib::spline1dcalc(cos_spline, std::fmod(std::abs(x), 2.*pi));
	}
	
	double sin(double x) const
	{
		if (x > 0.)
			return alglib::spline1dcalc(sin_spline, std::fmod(x, 2.*pi));
		else if (x < 0.)
			return -alglib::spline1dcalc(sin_spline, std::fmod(-x, 2.*pi));
		else
			return 0.;
	}
};
