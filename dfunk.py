import numpy as np
#from scipy.integrate import quad
from matplotlib import pyplot as plt
import os,inspect
import matplotlib.ticker as ticker
from matplotlib import gridspec
import ROOT
import rootpy
import rootpy.plotting as rplot
import rootpy.plotting.root2matplotlib as rplt

def plot_error_region(yvalues,yerrors,bins,color='orange'):

	for i in range(len(bins)-1):
		plt.fill_between(x=[bins[i],bins[i+1]],
			y1=yvalues[i]-yerrors[i],
			y2=yvalues[i]+yerrors[i],
			alpha=0.7,
			color=color)

def plot_error_region2(yvalues,yerrors,bins,color='orange'):
	bins_a = bins[:-1]
	bins_e = bins[1:]
	bins_m = (bins_a+bins_e)/2
	plt.errorbar(bins_m,yvalues,yerr=yerrors,drawstyle = 'steps-mid',color=color,lw=0.8)

def plot_error_region2ax(ax,yvalues,yerrors,bins,color='orange'):
	bins_a = bins[:-1]
	bins_e = bins[1:]
	bins_m = (bins_a+bins_e)/2
	ax.errorbar(bins_m,yvalues,yerr=yerrors,drawstyle = 'steps-mid',color=color,lw=0.8)

def parameters_kin(sam,var,par):
	topmassup=400
	topmassdown=0
	bmassup=60
	bmassdown=0
	umassup=1
	umassdown=0
	wmassup=100
	wmassdown=0
	jmassup=60
	jmassdown=0

	default_range_up=800
	default_range_down=0

	default_nbins = 64
	mass_nbins = 64
	nbin = default_nbins
	lower_range=default_range_down
	upper_range=default_range_up
	up_l = 0.1
	#plt.style.use("atlas.mplstyle")
	if((par == "Jet") and (var=="PT")):
		upper_range=500
		lower_range=default_range_down
		nbin=default_nbins
	if((sam=="pro") | (sam=="int")):
		upper_range=800
	if(par == "TopQuark" and var=="PT"):
		upper_range=default_range_up
		lower_range=default_range_down
		nbin=64
	if (var == "Eta" or var == "Phi"):
		lower_range=-5
		upper_range=5
		nbin = 32
		up_l = 0
	if (var == "Eta" and (par=="Photon" or par=="bJet" or par=="Jet")):
		lower_range=-2.5
		upper_range=2.5
		nbin = 32
		up_l = 0
	if(var=="M"):
		up_l=0
		nbin=mass_nbins
		if("W" in par):
			lower_range=wmassdown
			upper_range=wmassup
		if("Top" in par):
			lower_range=topmassdown
			upper_range=topmassup
		if("Jet" in par):
				lower_range=jmassdown
				upper_range=jmassup
		if("bJet" in par):
			lower_range=bmassdown
			upper_range=bmassup
		if("UQuark" in par):
			lower_range=umassdown
			upper_range=umassup
	if(var=="PT"):
		if("bJet" in par):
			upper_range=400
		if("W" in par):
			nbin = 32
		if("Top" in par):
			nbin = 32
		if("Jet" in par):
			nbin = 32
		if("bJet" in par):
			nbin = 32
	return [lower_range,upper_range,nbin,up_l]



def parameters_RM(p1,p2,va):
	up_l=0.001
	if va == "R":
		lower_range = 0
		upper_range = 5
		nbin = 32
		if p1 == "TopQuark" and p2 == "WBoson":
			lower_range = 0
			upper_range = 4
	
	if va == "M":
		nbin = 64
		lower_range = 0
		if p1 =="Photon":
			lower_range = 0
			upper_range = 1400
			if p2 == "TopQuark":
				lower_range = 0
				upper_range = 1400
				nbin = 32
			if p2 == "LeadingJet":
				lower_range = 0
				upper_range = 700
			if p2 == "bJet":
				lower_range = 0
				upper_range = 500
				nbin = 32
			if p2 == "WBoson":
				nbin = 32
		if p1 =="TopQuark":
			lower_range = 0
			upper_range = 800
			if p2 == "WBoson":
				lower_range = 50
				nbin = 32
			if p2 == "Photon":
				lower_range = 0
				upper_range = 1400
				nbin = 32
			if p2 == "LeadingJet":
				upper_range = 1200
				nbin = 32
	return [lower_range,upper_range,nbin,up_l]


def labelkin(va,par):
	label = r"$"
	if va == "PT":
		label += r"p_T"
	if va == "Eta":
		label += r"\eta"
	if va == "Phi":
		label += r"\Phi"
	if va == "M":
		label += r"m"
	if par == "Photon":
		label += r"(\gamma)$"
	if par == "TopQuark":
		label += r"(t)$"
	if par == "bJet":
		label += r"(bJet)$"
	if par == "WBoson":
		label += r"(W)$"
	if par == "Jet":
		label += r"(Jet)$"
	if va == "PT" or va == "M":
		label += "/GeV"
	plt.xlabel(label)
	plt.ylabel("number of events (normalized)/GeV",rotation=90)
	ax = plt.gca()
	ax.yaxis.set_label_coords(-0.08,0.5)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
	plt.grid(alpha=0.5)
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

	if va == "M":
		plt.legend(loc="upper right")
	else:
		plt.legend(loc="best")
	if va == "Eta":
		plt.legend(loc="best",fontsize="small")	

def labelkinax(ax,ax2,va,par):
	label = r"$"
	if va == "PT":
		label += r"p_T"
	if va == "Eta":
		label += r"\eta"
	if va == "Phi":
		label += r"\Phi"
	if va == "M":
		label += r"m"
	if par == "Photon":
		label += r"(\gamma)$"
	if par == "TopQuark":
		label += r"(t)$"
	if par == "bJet":
		label += r"(bJet)$"
	if par == "WBoson":
		label += r"(W)$"
	if par == "Jet":
		label += r"(Jet)$"
	if va == "PT" or va == "M":
		label += "/GeV"
	#ax.xlabel(label)
	ax.set_ylabel("number of events (normalized)/GeV",rotation=90)
	ax.yaxis.set_label_coords(-0.08,0.5)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
	ax.grid(alpha=0.5)
	ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	ax2.set_xlabel(label)
	ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator())
	ax2.grid(alpha=0.5)
	ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

	if va == "M":
		ax.legend(loc="upper right")
	else:
		ax.legend(loc="best")
	if va == "Eta":
		ax.legend(loc="best",fontsize="small")		


def labelRM(va,p1,p2):
	vari = r"$"
	if va == "R":
		vari += r"\Delta R"
	if va == "M":
		vari += r"m"

	parts = ""

	if p1 == "Photon":
		parts += r"\gamma,"
	if p1 == "TopQuark":
		parts += r"t,"

	if p2 == "Photon":
		parts += r"\gamma"
	if p2 == "TopQuark":
		parts += r"t"
	if p2 == "bJet":
		parts += r"bJet"
	if p2 == "WBoson":
		parts += r"W"
	if p2 == "LeadingJet":
		parts += r"Jet"

	einh = ""
	if va == "M":
		einh += "/GeV"

	plt.xlabel(vari+"("+parts+")"+r"$"+einh)
	plt.ylabel("number of events (normalized)/GeV",rotation=90)
	ax = plt.gca()
	ax.yaxis.set_label_coords(-0.08,0.5)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
	plt.grid(alpha=0.5)
	#if p1 == "Photon" and p2 == "LeadingJet" and va=="M":
	#	pass
	#else:
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")


def labelRMax(ax,ax2,va,p1,p2):
	vari = r"$"
	if va == "R":
		vari += r"\Delta R"
	if va == "M":
		vari += r"m"

	parts = ""

	if p1 == "Photon":
		parts += r"\gamma,"
	if p1 == "TopQuark":
		parts += r"t,"

	if p2 == "Photon":
		parts += r"\gamma"
	if p2 == "TopQuark":
		parts += r"t"
	if p2 == "bJet":
		parts += r"bJet"
	if p2 == "WBoson":
		parts += r"W"
	if p2 == "LeadingJet":
		parts += r"Jet"

	einh = ""
	if va == "M":
		einh += "/GeV"

	#ax.xlabel(vari+"("+parts+")"+r"$"+einh)
	ax.set_ylabel("number of events (normalized)/GeV",rotation=90)
	ax.yaxis.set_label_coords(-0.08,0.5)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
	ax.grid(alpha=0.5)
	#if p1 == "Photon" and p2 == "LeadingJet" and va=="M":
	#	pass
	#else:
	ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	ax.legend(loc="best")
	ax2.set_xlabel(vari+"("+parts+")"+r"$"+einh)
	ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator())
	ax2.grid(alpha=0.5)
	ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))



def set_dyn_binning(va,lo,up,n,err=0.15):
	#eps=0.0001
	#va = va[(np.abs(va)>up_l) & (va!=999.9)]
	#va = np.clip(va,lo+eps,up-eps)
	binning = np.linspace(lo,up,n)
	h1 = rplot.Hist(binning)
	map(h1.Fill,va)
	h1.Scale(1/(h1.Integral(0,h1.GetNbinsX()+1)))
	#print(binning)

	for i in range(n-1,1,-1):
		#print(i," ",binning[i])
		#if (h1.GetBinContent(i)==0):
			#print(0)
		#	continue
		#print(h1.GetBinError(i)/h1.GetBinContent(i))

		if (h1.GetBinContent(i)==0):
			binning = np.delete(binning, i-1)
			h1 = rplot.Hist(binning)
			map(h1.Fill,va)
			h1.Scale(1/(h1.Integral(0,h1.GetNbinsX()+1)))
			continue

		if (h1.GetBinError(i)/h1.GetBinContent(i)) > err:
			binning = np.delete(binning, i-1)
			h1 = rplot.Hist(binning)
			map(h1.Fill,va)
			h1.Scale(1/(h1.Integral(0,h1.GetNbinsX()+1)))
			continue

	return binning




def KMS(a,b):
	na = len(a)
	nb = len(b)
	a = np.sort(a)
	b = np.sort(b)
	h = tm.KolmogorovTest(na,a,nb,b,"")
	return h