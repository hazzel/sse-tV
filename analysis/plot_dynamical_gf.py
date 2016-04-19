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
	return a + b*np.exp(-c*x)

def combinatorial_factor(n, k):
	prod = Decimal(1)
	for j in range(n + 1):
		if k != j:
			kDec = Decimal(str(k))
			jDec = Decimal(str(j))
			prod *= (kDec + jDec) * (kDec - jDec)
	return Decimal(1) / prod

def estimator(n, beta, C, parity):
	omega1 = Decimal(str((2. + (1.-parity)/2.) * np.pi / beta))
	sum1 = Decimal(0)
	sum2 = Decimal(0)
	for k in range(n + 1):
		sum1 += Decimal(str(k**2.)) * combinatorial_factor(n, k) * Decimal(str(C[k]))
		sum2 += combinatorial_factor(n, k) * Decimal(str(C[k]))
	return float(omega1 * abs(- sum1 / sum2).sqrt())

def parse_ed_file(filename):
	ed_data = []
	ed_lines = ed_file.read().splitlines()
	for line in ed_lines:
		ed_data.append([])
		values = filter(None, line.split("\t"))
		for x in values:
			ed_data[-1].append(float(x))
		ed_data[-1] = np.array(ed_data[-1])
	return ed_data

latexify()
color_cycle = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'darkgreen']
marker_cycle = ['o', 'D', '<', 'p', '>', 'v', '*', '^', 's']

filelist = []
filelist.append(glob.glob("../bin/job/*.out"))
#filelist.append(glob.glob("/net/home/lxtsfs1/tpc/hesselmann/cluster_work/code/sse-tV/jobs/spectroscopy/job-L2-V1.0-T0.15/*task*.out"))
#filelist.append(glob.glob("/net/home/lxtsfs1/tpc/hesselmann/cluster_work/code/sse-tV/jobs/spectroscopy/job-L3-V1.0-T0.10/*task*.out"))
#filelist.append(glob.glob("/net/home/lxtsfs1/tpc/hesselmann/cluster_work/code/sse-tV/jobs/spectroscopy/job-L4-V1.0-T0.10/*task*.out"))
#filelist.append(glob.glob("/net/home/lxtsfs1/tpc/hesselmann/cluster_work/code/sse-tV/jobs/spectroscopy/job-L2-V0.5-T0.5/*task*.out"))
#filelist.append(glob.glob("/net/home/lxtsfs1/tpc/hesselmann/cluster_work/code/sse-tV/jobs/spectroscopy/job-L2-V0.5-T1.0/*task*.out"))
filelist.sort()

filelist = [item for sublist in filelist for item in sublist]
for f in filelist:
	figure, (ax1, ax2, ax3) = plt.subplots(1, 3)
	plist = ParseParameters(f)
	elist = ParseEvalables(f)

	obs = "M2"
	if obs == "M2":
		ed_n = 1
		parity = 1.
	if obs == "kekule":
		ed_n = 3
		parity = 1.
	elif obs == "epsilon":
		ed_n = 5
		parity = 1.
	elif obs == "epsilon_nn":
		ed_n = 7
		parity = 1.
		obs = "epsilon"
	elif obs == "sp":
		ed_n = 11
		parity = 1.
	elif obs == "tp":
		ed_n = 13
		parity = 1.
		
	for i in range(len(plist)):
		n_matsubara = int(plist[i]["matsubara_freqs"])
		n_discrete_tau = 2*int(plist[i]["discrete_tau"])
		h = float(plist[i]["V"])
		T = float(plist[i]["T"])
		L = float(plist[i]["L"])
		ed_glob = glob.glob("../../ctint/data/ed*" + "L_" + str(int(L)) + "*V_" + format(h, '.6f') + "*T_" + format(T, '.6f') + "*")
		if len(ed_glob) > 0:
			ed_file = open(ed_glob[0])
			ed_data = parse_ed_file(ed_file)
			n_ed_tau = int(ed_data[0][8])
			n_ed_mat = int(ed_data[0][9])

		figure.suptitle(r"$L = " + str(L) + ",\ V = " + str(h) + ",\ T = " + str(T) + "$")
		
		x_mat = (np.array(range(0, n_matsubara)) * 2. + (1.-parity)/2.) * np.pi * T
		y_mat = np.array(ArrangePlot(elist[i], "dyn_"+obs+"_mat")[0])
		err_mat = np.array(ArrangePlot(elist[i], "dyn_"+obs+"_mat")[1])
		x_tau = np.array(range(0, n_discrete_tau + 1)) / float(n_discrete_tau) / T
		y_tau = np.abs(np.array(ArrangePlot(elist[i], "dyn_"+obs+"_tau")[0]))
		err_tau = np.array(ArrangePlot(elist[i], "dyn_"+obs+"_tau")[1])

		N_bootstrap = 25
		x_delta = np.array(range(1, n_matsubara))
		y_delta = []
		for j in range(N_bootstrap):
			y_delta.append(np.zeros(n_matsubara - 1))
			y_boot = np.zeros(n_matsubara)
			for k in range(len(y_boot)):
				y_boot[k] = y_mat[k] + np.random.normal(0., 0.01 * abs(y_mat[k]))
			for n in range(1, n_matsubara):
				y_delta[j][n-1] = estimator(n, 1./T, y_boot, parity)

		ax1.set_xlabel(r"$\omega_n$")
		ax1.set_ylabel(r"$M_2(\omega_n) \cdot \omega_n^2$")
		xscale = np.copy(x_mat)
		xscale[0] = 1.
		xscale = xscale**2.
		ax1.plot(x_mat, y_mat * xscale, marker="o", color="green", markersize=10.0, linestyle="None", linewidth=0.0)
		(_, caps, _) = ax1.errorbar(x_mat, y_mat * xscale, yerr=err_mat * xscale, marker='None', capsize=8, color="green")
		for cap in caps:
			cap.set_markeredgewidth(1.4)
		if len(ed_glob) > 0:
			x_mat = (np.array(range(0, n_ed_mat)) * 2. + (1.-parity)/2.) * np.pi * T
			xscale = np.copy(x_mat)
			xscale[0] = 1.
			xscale = xscale**2.
			ax1.plot(x_mat, ed_data[ed_n+1] * xscale, marker='o', color="r", markersize=10.0, linewidth=2.0)
			
			ed_tau = np.linspace(0., 1./T, n_ed_tau + 1)
			y_num_int = np.zeros(n_ed_mat)
			for n in range(n_ed_mat):
				y_num_int[n] = scipy.integrate.simps(ed_data[ed_n] * np.cos(ed_tau * x_mat[n]), ed_tau)
			ax1.plot(x_mat, y_num_int * xscale, marker='o', color="b", markersize=10.0, linewidth=2.0)

		ax2.set_xlabel(r"$\tau$")
		ax2.set_ylabel(r"$M_2(\tau)$")
		ax2.set_yscale("log")
		ax2.plot(x_tau, y_tau, marker="o", color="green", markersize=10.0, linewidth=2.0, label=r'$L='+str(int(L))+'$')
		(_, caps, _) = ax2.errorbar(x_tau, y_tau, yerr=err_tau, marker='None', capsize=8, color="green")
		for cap in caps:
			cap.set_markeredgewidth(1.4)
		if len(ed_glob) > 0:
			ax2.plot(np.linspace(0., 1./T, n_ed_tau + 1), ed_data[ed_n], marker='o', color="r", markersize=10.0, linewidth=2.0, label=r'$L='+str(int(L))+'$')
			ax2.plot(np.linspace(0., 1./T, n_ed_tau + 1), np.flipud(ed_data[ed_n]), marker='o', color="orange", markersize=10.0, linewidth=2.0, label=r'$L='+str(int(L))+'$')
		
		try:
			nmin = len(x_tau)*3/32; nmax = len(x_tau)*10/32
			#nmin = len(x_tau)*10/16; nmax = len(x_tau)*14/16
			#nmin = len(x_tau)*17/32; nmax = len(x_tau)*31/32
			#nmin = 0; nmax = len(x_tau)*2/16
			#parameter, perr = fit_function( [0.1, 0.1, 1.], x_tau[nmin:nmax], y_tau[nmin:nmax], FitFunction, datayerrors=err_tau[nmin:nmax])
			parameter, perr = curve_fit(FitFunction, x_tau[nmin:nmax], y_tau[nmin:nmax], p0=[0.01, 0.01, 1.])
			px = np.linspace(x_tau[nmin], x_tau[nmax], 1000)
			ax2.plot(px, FitFunction(px, *parameter), 'k-', linewidth=3.0)
			d = -int(np.log10(abs(perr[2])))+2
			ax2.text(0.10, 0.98, r"$\Delta_{FIT} = " + ("{:."+str(d)+"f}").format(parameter[2]) + "(" + str(round(perr[2], d)*10.**d).partition('.')[0] + ")$", transform=ax2.transAxes, fontsize=20, va='top')
			#ax2.text(0.05, 0.92, r"$\Delta_{ED} = 0.9264$", transform=ax2.transAxes, fontsize=20, va='top')
			print parameter
			print perr
		except:
			print "runtime error"
		
		ax3.set_xlabel(r"$n$")
		ax3.set_ylabel(r"$\Delta_n$")
		ax3.plot(x_delta, np.mean(y_delta, axis=0), marker="o", color="green", markersize=10.0, linewidth=2.0)
		(_, caps, _) = ax3.errorbar(x_delta, np.mean(y_delta, axis=0), yerr=np.std(y_delta, axis=0), marker="None", color="green", markersize=10.0, linewidth=2.0)
		for cap in caps:
			cap.set_markeredgewidth(1.4)

		if len(ed_glob) > 0:
			x_delta = np.array(range(1, n_ed_mat))
			y_delta = np.zeros(n_ed_mat - 1)
			for n in range(1, n_ed_mat):
				y_delta[n-1] = estimator(n, 1./T, ed_data[ed_n+1], parity)
			ax3.plot(x_delta, y_delta, marker="o", color="red", markersize=10.0, linewidth=2.0, label=r'$L='+str(int(L))+'$')
		
	plt.tight_layout()
plt.show()
