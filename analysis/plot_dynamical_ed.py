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

def FitFunctionL(x, a, b, c):
	return a + b*np.exp(-c*x)
	#return 2.85*10.**-9. + b*np.exp(-c*x)

def FitFunctionR(x, a, b, c):
	return a + b*np.exp(c*x)

def ExpSumFunction(x, a, b, c, d):
	return a + b*np.exp(-d*x) + c*np.exp(d*x)

def LinearFunction(x, a, b):
	return a - b*x

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


obs = "epsilon"
if obs == "M2":
	ed_n = 1
	parity = 1.
if obs == "kekule":
	ed_n = 3
	parity = 1.
elif obs == "chern":
	ed_n = 5
	parity = 1.
elif obs == "epsilon":
	ed_n = 7
	parity = 1.
elif obs == "epsilon_nn":
	ed_n = 9
	parity = 1.
	obs = "epsilon"
elif obs == "sp":
	ed_n = 13
	parity = 1.
elif obs == "tp":
	ed_n = 15
	parity = 1.

T = 0.01
h = 1.355
L = 2
ed_glob = glob.glob("../../ctint/data/ed*" + "L_" + str(int(L)) + "*V_" + format(h, '.6f') + "*T_" + format(T, '.6f') + "*")
if len(ed_glob) > 0:
	ed_file = open(ed_glob[0])
	ed_data = parse_ed_file(ed_file)
	n_ed_tau = int(ed_data[0][8])
	n_ed_mat = int(ed_data[0][9])
	ed_tau = np.linspace(0., 1./T, n_ed_tau + 1)
	
	if obs == "epsilon":
		if T == 0.001:
			e0 = -0.44975132
		elif T == 0.01:
			e0 = -0.44975132
		elif T == 0.025:
			e0 = -0.44975132
		elif T == 0.075:
			e0 = -0.44975067
		elif T == 0.10:
			e0 = -0.4497369
		elif T == 0.20:
			e0 = -0.44821574
		elif T == 0.20:
			e0 = -0.44203009
		else:
			e0 = 0.

		#ed_y = np.abs(ed_data[ed_n] - e0**2.)
		ed_y = np.abs(ed_data[ed_n])

figure, ax2 = plt.subplots(1, 1)
figure.suptitle(r"$L = " + str(L) + ",\ V = " + str(h) + ",\ T = " + str(T) + "$")

ax2.set_xlabel(r"$\tau$")
ax2.set_ylabel(r"$"+obs+r"(\tau)$")
ax2.set_yscale("log")
if len(ed_glob) > 0 and len(ed_data) > ed_n:
	ax2.plot(ed_tau, ed_y, marker='o', color="r", markersize=10.0, linewidth=0.0, label=r'$L='+str(int(L))+'$')
	#ax2.plot(ed_tau, np.flipud(ed_data[ed_n]), marker='o', color="r", markersize=10.0, linewidth=2.0, label=r'$L='+str(int(L))+'$')

d = 3
if len(ed_glob) > 0 and len(ed_data) > ed_n:
	#parameter_ed, perr_ed = scipy.optimize.curve_fit( ExpSumFunction, ed_tau, ed_data[ed_n], p0=[0.1, 0.1, 0.1, 1.0])
	#px = np.linspace(ed_tau[0], ed_tau[len(ed_data[ed_n])-1], 1000)
	#ax2.plot(px, ExpSumFunction(px, *parameter_ed), 'r-', linewidth=3.0)
	#ax2.text(0.10, 0.93, r"$\Delta_{FIT\ ED} = " + ("{:."+str(d)+"f}").format(abs(parameter_ed[3])) + "$", transform=ax2.transAxes, fontsize=20, va='top')
	
	if obs == "epsilon":
		n_min = 1
		parameter_ed, perr_ed = scipy.optimize.curve_fit( FitFunctionL, ed_tau[n_min:len(ed_data[ed_n])*16/32], ed_y[n_min:len(ed_data[ed_n])*16/32], p0=[28., 2.5, 3.])
	elif obs == "sp":
		parameter_ed, perr_ed = scipy.optimize.curve_fit( FitFunctionL, ed_tau[len(ed_data[ed_n])/4:len(ed_data[ed_n])/2], ed_data[ed_n][len(ed_data[ed_n])/4:len(ed_data[ed_n])/2], p0=[0.1, 0.1, 1.])
		#parameter_ed, perr_ed = scipy.optimize.curve_fit( FitFunctionR, ed_tau[len(ed_data[ed_n])/2:len(ed_data[ed_n])-15], ed_data[ed_n][len(ed_data[ed_n])/2:len(ed_data[ed_n])-15], p0=[0.1, 0.1, 1.])
	else:
		parameter_ed, perr_ed = scipy.optimize.curve_fit( FitFunctionL, ed_tau[3:len(ed_data[ed_n])*12/32], ed_data[ed_n][3:len(ed_data[ed_n])*12/32], p0=[0.1, 0.1, 1.])
	#px = np.linspace(ed_tau[0], ed_tau[len(ed_data[ed_n])/2], 1000)
	px = np.linspace(ed_tau[n_min], ed_tau[len(ed_data[ed_n])*16/32], 1000)
	#px = np.linspace(ed_tau[len(ed_data[ed_n])/2], ed_tau[len(ed_data[ed_n])-1], 1000)
	ax2.plot(px, FitFunctionL(px, *parameter_ed), 'r-', linewidth=3.0)
	ax2.text(0.10, 0.93, r"$\Delta_{FIT\ ED} = " + ("{:."+str(d)+"f}").format(abs(parameter_ed[2])) + "$", transform=ax2.transAxes, fontsize=20, va='top')
	print parameter_ed

plt.show()
