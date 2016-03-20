import re
import os
import glob
import sys
sys.path.append('/home/stephan/mc/ctqmc')
sys.path.append("/net/home/lxtsfs1/tpc/hesselmann/mc/ctqmc")
import numpy as np
from cdecimal import *
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import pylab
from ParseDataOutput import *
sys.path.append("/net/home/lxtsfs1/tpc/hesselmann/mc/qising-SSE")
sys.path.append("/home/stephan/mc/qising-SSE")
from Fit import *
from texify import *
import scipy.integrate

def FitFunction(x, a, b, c):
	return a+b*np.exp(-c*x)

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
#filelist.append(glob.glob("../data/dyn_M2/*.out"))
#filelist.append(glob.glob("/net/home/lxtsfs1/tpc/hesselmann/cluster_work/code/ctint/jobs/spectroscopy/L2-V1.355-T0.02/*task*.out"))
filelist.append(glob.glob("/net/home/lxtsfs1/tpc/hesselmann/cluster_work/code/ctint/jobs/spectroscopy/job-L2-V1.355-T0.20/*task*.out"))
filelist.sort()
ed_data = pylab.loadtxt(glob.glob("../../ctint/data/ed*")[0])

filelist = [item for sublist in filelist for item in sublist]
for f in filelist:
	figure, (ax1, ax2, ax3) = plt.subplots(1, 3)
	plist = ParseParameters(f)
	elist = ParseEvalables(f)
		
	for i in range(len(plist)):
		n_matsubara = int(plist[i]["matsubara_freqs"])
		n_discrete_tau = int(plist[i]["discrete_tau"])
		h = float(plist[i]["V"])
		T = float(plist[i]["T"])
		L = float(plist[i]["L"])
		for j in range(len(ed_data)):
			if h == ed_data[j,2] and T == ed_data[j,3] and L == int(ed_data[j,1]):
				n_ed_tau = int(ed_data[j,8])
				n_ed_mat = int(ed_data[j,10+n_ed_tau])
		figure.suptitle(r"$L = " + str(L) + ",\ V = " + str(h) + ",\ T = " + str(T) + "$")
	
		x_tau = np.array(range(0, n_discrete_tau + 1)) / float(n_discrete_tau) / T
		y_tau = np.array(ArrangePlot(elist[i], "dynamical_M2_tau")[0])
		err_tau = np.array(ArrangePlot(elist[i], "dynamical_M2_tau")[1])

		ax1.set_xlabel(r"$\omega_n$")
		ax1.set_ylabel(r"$M_2(\omega_n) \cdot \omega_n^2$")
		for j in range(len(ed_data)):
			if h == ed_data[j,2] and T == ed_data[j,3] and L == int(ed_data[j,1]):
				ax1.plot(np.array(range(0, n_ed_mat)) * 2. * np.pi * T, ed_data[j,11+n_ed_tau:] * ((np.array(range(0, n_ed_mat)) * 2.) * np.pi * T)**2., marker='o', color="r", markersize=10.0, linewidth=2.0)
			
		ax2.set_xlabel(r"$\tau$")
		ax2.set_ylabel(r"$M_2(\tau)$")
		ax2.set_yscale("log")
		ax2.plot(x_tau, y_tau, marker="o", color="green", markersize=10.0, linewidth=2.0, label=r'$L='+str(int(L))+'$')
		(_, caps, _) = ax2.errorbar(x_tau, y_tau, yerr=err_tau, marker='None', capsize=8, color="green")
		for cap in caps:
			cap.set_markeredgewidth(1.4)
		for j in range(len(ed_data)):
			if h == ed_data[j,2] and T == ed_data[j,3] and L == int(ed_data[j,1]):
				ax2.plot(np.linspace(0., 1./T/2., n_ed_tau + 1), ed_data[j,9:10+n_ed_tau], marker='o', color="r", markersize=10.0, linewidth=2.0, label=r'$L='+str(int(L))+'$')
		
		try:
			nmin = 0; nmax = len(x_tau)/2
			parameter, perr = fit_function( [0.0, 1., 1.], x_tau[nmin:nmax], y_tau[nmin:nmax], FitFunction, datayerrors=err_tau[nmin:nmax])
			px = np.linspace(x_tau[nmin], x_tau[nmax], 1000)
			ax2.plot(px, FitFunction(px, *parameter), 'k-', linewidth=3.0)
			d = -int(np.log10(abs(perr[2])))+2
			ax2.text(0.05, 0.98, r"$\Delta_{FIT} = " + ("{:."+str(d)+"f}").format(parameter[2]) + "(" + str(round(perr[2], d)*10.**d).partition('.')[0] + ")$", transform=ax2.transAxes, fontsize=20, va='top')
			ax2.text(0.05, 0.92, r"$\Delta_{ED} = 0.9264$", transform=ax2.transAxes, fontsize=20, va='top')
			print parameter
			print perr
		except:
			print "runtime error"
		
		ax3.set_xlabel(r"$n$")
		ax3.set_ylabel(r"$\Delta_n$")
		for j in range(len(ed_data)):
			if h == ed_data[j,2] and T == ed_data[j,3] and L == int(ed_data[j,1]):
				x_delta = np.array(range(1, n_ed_mat))
				y_delta = np.zeros(n_ed_mat - 1)
				for n in range(1, n_ed_mat):
					for j in range(len(ed_data)):
						if h == ed_data[j,2] and T == ed_data[j,3] and L == int(ed_data[j,1]):
							y_delta[n-1] = estimator(n, 1./T, ed_data[j,11+n_ed_tau:])
				ax3.plot(x_delta, y_delta, marker="o", color="red", markersize=10.0, linewidth=2.0, label=r'$L='+str(int(L))+'$')
		
	plt.tight_layout()
plt.show()
