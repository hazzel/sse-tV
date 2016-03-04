import re
import os
import glob
import sys
sys.path.append('/home/stephan/mc/ctqmc')
sys.path.append("/net/home/lxtsfs1/tpc/hesselmann/mc/ctqmc")
import numpy as np
from decimal import *
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
from ParseDataOutput import *
sys.path.append("/net/home/lxtsfs1/tpc/hesselmann/mc/qising-SSE")
from Fit import *
from texify import *

def FitFunction(x, a, b, c):
	return a+b*np.exp(c*x)

def combinatorial_factor(n, k):
	prod = Decimal(1)
	for j in range(n + 1):
		if k != j:
			kDec = Decimal(str(k))
			jDec = Decimal(str(j))
			prod *= (kDec + jDec) * (kDec - jDec)
	return Decimal(1) / prod

def estimator(n, beta, C):
	omega1 = Decimal(str(2. * np.pi / beta))
	sum1 = Decimal(0)
	sum2 = Decimal(0)
	for k in range(n + 1):
		sum1 += Decimal(str(k**2.)) * combinatorial_factor(n, k) * Decimal(str(C[k]))
		sum2 += combinatorial_factor(n, k) * Decimal(str(C[k]))
	return float(omega1 * abs(- sum1 / sum2).sqrt())

latexify()
color_cycle = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'darkgreen']
marker_cycle = ['o', 'D', '<', 'p', '>', 'v', '*', '^', 's']

filelist = []
filelist.append(glob.glob("../bin/job/*.out"))
filelist.append(glob.glob("../data/dyn_M2/*.out"))
filelist.sort()

filelist = [item for sublist in filelist for item in sublist]
for f in filelist:
	figure, (ax1, ax2, ax3) = plt.subplots(1, 3)
	plist = ParseParameters(f)
	elist = ParseEvalables(f)
		
	for i in range(len(plist)):
		n_matsubara = int(plist[i]["matsubara_freqs"])
		n_discrete_tau = int(plist[i]["discrete_tau"]) + 1
		h = float(plist[i]["V"])
		T = float(plist[i]["T"])
		L = float(plist[i]["L"])
		figure.suptitle(r"$L = " + str(L) + ",\ V = " + str(h) + ",\ T = " + str(T) + "$")
	
		x_mat = (np.array(range(0, n_matsubara)) * 2.) * np.pi * T
		y_mat = np.array(ArrangePlot(elist[i], "dynamical_M2_mat")[0])
		err_mat = np.array(ArrangePlot(elist[i], "dynamical_M2_mat")[1])
		x_tau = np.array(range(0, n_discrete_tau-1)) / float(n_discrete_tau) / T
		y_tau = np.array(ArrangePlot(elist[i], "dynamical_M2_tau")[0])[:-1]
		err_tau = np.array(ArrangePlot(elist[i], "dynamical_M2_tau")[1])[:-1]
		#n_matsubara = 100
		x_delta = np.array(range(1, n_matsubara))
		y_delta = []
		for n in range(1, n_matsubara):
			y_delta.append(estimator(n, 1./T, y_mat))

		c = 0
		ax1.set_xlabel(r"$\omega_n$")
		ax1.set_ylabel(r"$M_2(\omega_n) \cdot \omega_n^2$")
		ax1.plot(x_mat, y_mat * x_mat**2., marker=marker_cycle[c%len(marker_cycle)], color=color_cycle[c%len(color_cycle)], markersize=10.0, linewidth=2.0, label=r'$L='+str(int(L))+'$')
		(_, caps, _) = ax1.errorbar(x_mat, y_mat * x_mat**2., yerr=err_mat* x_mat**2., marker='None', capsize=8, color=color_cycle[c%len(color_cycle)])
		for cap in caps:
			cap.set_markeredgewidth(1.4)
			
		c = 1
		ax2.set_xlabel(r"$\tau$")
		ax2.set_ylabel(r"$M_2(\tau)$")
		ax2.set_yscale("log")
		ax2.plot(x_tau, y_tau, marker=marker_cycle[c%len(marker_cycle)], color=color_cycle[c%len(color_cycle)], markersize=10.0, linewidth=2.0, label=r'$L='+str(int(L))+'$')
		(_, caps, _) = ax2.errorbar(x_tau, y_tau, yerr=err_tau, marker='None', capsize=8, color=color_cycle[c%len(color_cycle)])
		for cap in caps:
			cap.set_markeredgewidth(1.4)

		'''
		parameter, perr = fit_function( [0.0, 0.05, 9.], x_tau, y_tau, FitFunction, datayerrors=err_tau)
		px = np.linspace(x_tau[0], x_tau[-1], 1000)
		ax2.plot(px, FitFunction(px, *parameter), 'k-', linewidth=3.0)
		d = -int(np.log10(abs(perr[2])))+2
		ax2.text(0.05, 0.98, r"$\Delta_{FIT} = " + ("{:."+str(d)+"f}").format(parameter[2]) + "(" + str(round(perr[2], d)*10.**d).partition('.')[0] + ")$", transform=ax2.transAxes, fontsize=20, va='top')
		ax2.text(0.05, 0.92, r"$\Delta_{ED} = 0.9264$", transform=ax2.transAxes, fontsize=20, va='top')
		print parameter
		print perr
		'''
		
		c = 2
		ax3.set_xlabel(r"$n$")
		ax3.set_ylabel(r"$\Delta_n$")
		ax3.plot(x_delta, y_delta, marker=marker_cycle[c%len(marker_cycle)], color=color_cycle[c%len(color_cycle)], markersize=10.0, linewidth=2.0, label=r'$L='+str(int(L))+'$')
		
	plt.tight_layout()
plt.show()
