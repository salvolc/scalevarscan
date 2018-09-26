from tfunk import *

samples = ["dec","int","pro"]
vari = ["PT","Eta","Phi","M"]
part = ["TopQuark","Photon","BQuark","WBoson","UQuark"]


PREFIX = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+"/"

datas = os.listdir(PREFIX+"data")

ntruth = 0
for i in range(len(datas)):
	if("truth" in datas[i]):
		ntruth += 1

for par in part:
	if(par+"_truth" not in os.listdir(PREFIX+"plots/")):
		os.chdir(PREFIX+"plots/")
		os.mkdir(par+"_truth/")
		os.chdir(PREFIX)



ratiox=6
ratioy=4
eps=0.001


for par in part:
	for var in vari:
		if(par == "Photon" and var=="M"):
			continue
		para = parameters_kin("dec", var, par)
		lower_range=para[0];upper_range=para[1];nbin=para[2];up_l=para[3]

		evdyn 		= np.genfromtxt("data/"+"dyn"+"_"+par+"_"+var+"_truth.txt")
		evtop 		= np.genfromtxt("data/"+"top"+"_"+par+"_"+var+"_truth.txt")
		evtwotop	= np.genfromtxt("data/"+"ttop"+"_"+par+"_"+var+"_truth.txt")
		evhalftop	= np.genfromtxt("data/"+"htop"+"_"+par+"_"+var+"_truth.txt")
		evtripeltop	= np.genfromtxt("data/"+"tritop"+"_"+par+"_"+var+"_truth.txt")
		evowntop	= np.genfromtxt("data/"+"owntop"+"_"+par+"_"+var+"_truth.txt")

		#evi = np.genfromtxt("data/"+"int"+"_"+par+"_"+var+"_truth.txt")
		#evd = np.genfromtxt("data/"+"dec"+"_"+par+"_"+var+"_truth.txt")
		#evp = np.genfromtxt("data/"+"pro"+"_"+par+"_"+var+"_truth.txt")

		evdyn 		= evdyn[(np.abs(evdyn)>up_l) & (evdyn!=999.9)];
		evtop 		= evtop[(np.abs(evtop)>up_l) & (evtop!=999.9)];
		evtwotop	= evtwotop[(np.abs(evtwotop)>up_l) & (evtwotop!=999.9)];
		evhalftop	= evhalftop[(np.abs(evhalftop)>up_l) & (evhalftop!=999.9)];
		evtripeltop	= evtripeltop[(np.abs(evtripeltop)>up_l) & (evtripeltop!=999.9)];
		evowntop	= evowntop[(np.abs(evowntop)>up_l) & (evowntop!=999.9)]

		evdyn 		 = np.clip(evdyn, lower_range+eps, upper_range-eps);
		evtop 		 = np.clip(evtop, lower_range+eps, upper_range-eps);
		evtwotop	 = np.clip(evtwotop, lower_range+eps, upper_range-eps);
		evhalftop	 = np.clip(evhalftop, lower_range+eps, upper_range-eps);
		evtripeltop	 = np.clip(evtripeltop, lower_range+eps, upper_range-eps);
		evowntop	 = np.clip(evowntop, lower_range+eps, upper_range-eps)

		gs = gridspec.GridSpec(2, 1, width_ratios=[1],height_ratios=[3,1]) 
		f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
		ax1 = plt.subplot(gs[0])
		ax2 = plt.subplot(gs[1])
		#plt.subplots_adjust(wspace=-1)
		#nbin=128
		#binning = np.linspace(lower_range, upper_range, nbin)
		binning = set_dyn_binning(evdyn, lower_range, upper_range, nbin)

		hdyn 		= rplot.Hist(binning);map(hdyn.Fill, evdyn)
		htop 		= rplot.Hist(binning);map(htop.Fill, evtop)
		htwotop		= rplot.Hist(binning);map(htwotop.Fill, evtwotop)
		hhalftop	= rplot.Hist(binning);map(hhalftop.Fill, evhalftop)
		htripeltop	= rplot.Hist(binning);map(htripeltop.Fill, evtripeltop)
		howntop		= rplot.Hist(binning);map(howntop.Fill, evowntop)

		hdyn.Sumw2(True);htop.Sumw2(True);htwotop.Sumw2(True);hhalftop.Sumw2(True);htripeltop.Sumw2(True);howntop.Sumw2(True)

		hdyn.Scale(1/(hdyn.Integral(0,hdyn.GetNbinsX()+1)),"width")
		hdyn.linecolor = "black";hdyn.linewidth = 1
		rplt.hist(hdyn,fmt="none",axes=ax1,lw=0.4,color="black",label="Dynamic scale")
		rplt.errorbar(hdyn,fmt="none",axes=ax1,lw=0.4,color="black",label="_nolegend_")

		htop.Scale(1/(htop.Integral(0,htop.GetNbinsX()+1)),"width")
		htop.linecolor = "green";htop.linewidth = 1
		rplt.hist(htop,fmt="none",axes=ax1,lw=0.4,color="green",label="Top Mass")
		rplt.errorbar(htop,fmt="none",axes=ax1,lw=0.4,color="green",label="_nolegend_")

		htwotop.Scale(1/(htwotop.Integral(0,htwotop.GetNbinsX()+1)),"width")
		htwotop.linecolor = "red";htwotop.linewidth = 1
		rplt.hist(htwotop,fmt="none",axes=ax1,lw=0.4,color="red",label="Double Top Mass")
		rplt.errorbar(htwotop,fmt="none",axes=ax1,lw=0.4,color="red",label="_nolegend_")

		hhalftop.Scale(1/(hhalftop.Integral(0,hhalftop.GetNbinsX()+1)),"width")
		hhalftop.linecolor = "blue";hhalftop.linewidth = 1
		rplt.hist(hhalftop,fmt="none",axes=ax1,lw=0.4,color="blue",label="Half Top Mass")
		rplt.errorbar(hhalftop,fmt="none",axes=ax1,lw=0.4,color="blue",label="_nolegend_")

		htripeltop.Scale(1/(htripeltop.Integral(0,htripeltop.GetNbinsX()+1)),"width")
		htripeltop.linecolor = "magenta";htripeltop.linewidth = 1
		#rplt.hist(htripeltop,fmt="none",axes=ax1,lw=0.4,color="magenta",label="Tripel Top Mass")
		#rplt.errorbar(htripeltop,fmt="none",axes=ax1,lw=0.4,color="magenta",label="_nolegend_")

		howntop.Scale(1/(howntop.Integral(0,howntop.GetNbinsX()+1)),"width")
		howntop.linecolor = "gray";howntop.linewidth = 1
		#rplt.hist(howntop,fmt="none",axes=ax1,lw=0.4,color="gray",label="Own Scale (top mass)")
		#rplt.errorbar(howntop,fmt="none",axes=ax1,lw=0.4,color="gray",label="_nolegend_")

		#print(par + " " + var + " KS Test")
		ksval = hdyn.KolmogorovTest(htop)
		#print(hi.KolmogorovTest(hdp))
		#print()

		f = open("scalevarexp/scale_"+var+"_"+par+"_truth.txt","w")

		for i in range(0,htop.GetNbinsX()+2):
			f.write(str(htop.GetBinLowEdge(i))+"\n")
			if htop.GetBinContent(i) == 0:
				f.write(str(0)+" ")
				f.write(str(0)+"\n")
				continue
			f.write(str(htwotop.GetBinContent(i)/htop.GetBinContent(i))+" ")
			f.write(str(hhalftop.GetBinContent(i)/htop.GetBinContent(i))+"\n")
		f.close()

		ax2.hist(binning[1:]-np.diff(binning)/2,bins=binning,weights=np.ones_like(binning[1:]),histtype="step",lw=0.8,color="black")
		htop.Divide(hdyn);htwotop.Divide(hdyn);hhalftop.Divide(hdyn);#htripeltop.Divide(hdyn);
		howntop.Divide(hdyn)

		rplt.errorbar(htop,fmt="none",axes=ax2,lw=0.6,color="green",label="_nolegend_")
		rplt.errorbar(htwotop,fmt="none",axes=ax2,lw=0.6,color="red",label="_nolegend_")
		rplt.errorbar(hhalftop,fmt="none",axes=ax2,lw=0.6,color="blue",label="_nolegend_")
		#rplt.errorbar(htripeltop,fmt="none",axes=ax2,lw=0.6,color="magenta",label="_nolegend_")
		#rplt.errorbar(howntop,fmt="none",axes=ax2,lw=0.6,color="gray",label="_nolegend_")

		#hi.Divide(hi)
		#rplt.errorbar(hi,fmt="none",axes=ax2,lw=0.6,color="black",label="_nolegend_")

		ax2.set_ylabel("Ratio to dyn. scale")
		ax2.grid(alpha=0.6)

		ax1.set_xlim(lower_range,upper_range)
		ax2.set_xlim(lower_range,upper_range)
		ax2.set_ylim(0.5,1.5)

		if par == "Photon" and var=="PT":
			ax1.set_ylim(ymin=0)
			ax2.set_ylim(0.5,2)

		labelkinax(ax1,ax2,var,par)
		leg = ax1.legend()#title="KS Test: "+"{:.7f}".format(ksval),fontsize="small")
		#leg.get_title().set_fontsize('small')
		
		left = lower_range
		right = upper_range
		bottom, top = ax1.get_ylim()
		ax1.text(0.7*(left+right),0.4*(bottom+top),r"$\sqrt{s}=13\,$ TeV"+"\n"+"truth level")

		plt.savefig("plots/"+par+"_truth/"+"scalevar"+"_"+par+"_"+var+".pdf",bbox_inches='tight')
		plt.close()

"""
wdyn=np.genfromtxt("data/pro_weight.txt");hwdyn=rplot.Hist(np.linspace(0,np.max(wdyn)*1.1,64));map(hwdyn.Fill, wdyn)
rplt.hist(hwdyn,fmt="none",lw=0.4,color="blue",label="_nolegend_")
rplt.errorbar(hwdyn,fmt="none",lw=0.6,color="blue",label="_nolegend_")
plt.xlabel("weights")
plt.ylabel("number of entries",rotation=90)
plt.grid(alpha=0.7)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.savefig("weights_dynscale.pdf",bbox_inches='tight')
plt.close()

wtop=np.genfromtxt("data/dec_weight.txt");hwtop=rplot.Hist(np.linspace(0,np.max(wtop)*1.1,64));map(hwtop.Fill, wtop)
rplt.hist(hwtop,fmt="none",lw=0.4,color="blue",label="_nolegend_")
rplt.errorbar(hwtop,fmt="none",lw=0.6,color="blue",label="_nolegend_")
plt.xlabel("weights")
plt.ylabel("number of entries",rotation=90)
plt.grid(alpha=0.7)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.savefig("weights_topscale.pdf",bbox_inches='tight')
plt.close()

wtwotop=np.genfromtxt("data/int_weight.txt");hwtwotop=rplot.Hist(np.linspace(0,np.max(wtwotop)*1.1,64));map(hwtwotop.Fill, wtwotop)
rplt.hist(hwtwotop,fmt="none",lw=0.4,color="blue",label="_nolegend_")
rplt.errorbar(hwtwotop,fmt="none",lw=0.6,color="blue",label="_nolegend_")
plt.xlabel("weights")
plt.ylabel("number of entries",rotation=90)
plt.grid(alpha=0.7)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.savefig("weights_twotopscale.pdf",bbox_inches='tight')
plt.close()

whalftop=np.genfromtxt("data/int_weight.txt");hwhalftop=rplot.Hist(np.linspace(0,np.max(whalftop)*1.1,64));map(hwhalftop.Fill, whalftop)
rplt.hist(hwhalftop,fmt="none",lw=0.4,color="blue",label="_nolegend_")
rplt.errorbar(hwhalftop,fmt="none",lw=0.6,color="blue",label="_nolegend_")
plt.xlabel("weights")
plt.ylabel("number of entries",rotation=90)
plt.grid(alpha=0.7)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.savefig("weights_halftopscale.pdf",bbox_inches='tight')
plt.close()"""




if("R_truth" not in os.listdir(PREFIX+"plots/")):
	os.chdir(PREFIX+"plots/")
	os.mkdir("R_truth/")
	os.chdir(PREFIX)

if("M_truth" not in os.listdir(PREFIX+"plots/")):
	os.chdir(PREFIX+"plots/")
	os.mkdir("M_truth/")
	os.chdir(PREFIX)

		
var=["R","M"]
RMPart=["Photon","TopQuark"]
RMPartP=["Photon","TopQuark","BQuark","WBoson","UQuark"]

up_l=0.01

for va in var:
	for p1 in RMPart:
		for p2 in RMPartP:
			if(p1==p2):
				continue

			para = parameters_RM(p1, p2, va)
			lower_range=para[0];upper_range=para[1];nbin=para[2];up_l=para[3]




			evdyn 		= np.genfromtxt("data/"+"dyn_"+p1+"_"+p2+"_"+va+"_truth.txt")
			evtop 		= np.genfromtxt("data/"+"top_"+p1+"_"+p2+"_"+va+"_truth.txt")
			evtwotop	= np.genfromtxt("data/"+"ttop_"+p1+"_"+p2+"_"+va+"_truth.txt")
			evhalftop	= np.genfromtxt("data/"+"htop_"+p1+"_"+p2+"_"+va+"_truth.txt")
			evtripeltop	= np.genfromtxt("data/"+"tritop_"+p1+"_"+p2+"_"+va+"_truth.txt")
			evowntop	= np.genfromtxt("data/"+"owntop_"+p1+"_"+p2+"_"+va+"_truth.txt")

			evdyn 		= evdyn[(np.abs(evdyn)>up_l) & (evdyn!=999.9)];
			evtop 		= evtop[(np.abs(evtop)>up_l) & (evtop!=999.9)];
			evtwotop	= evtwotop[(np.abs(evtwotop)>up_l) & (evtwotop!=999.9)];
			evhalftop	= evhalftop[(np.abs(evhalftop)>up_l) & (evhalftop!=999.9)];
			evtripeltop	= evtripeltop[(np.abs(evtripeltop)>up_l) & (evtripeltop!=999.9)]
			evowntop	= evowntop[(np.abs(evowntop)>up_l) & (evowntop!=999.9)]

			evdyn 		 = np.clip(evdyn, lower_range+eps, upper_range-eps);
			evtop 		 = np.clip(evtop, lower_range+eps, upper_range-eps);
			evtwotop	 = np.clip(evtwotop, lower_range+eps, upper_range-eps);
			evhalftop	 = np.clip(evhalftop, lower_range+eps, upper_range-eps);			
			evtripeltop	 = np.clip(evtripeltop, lower_range+eps, upper_range-eps)
			evowntop	 = np.clip(evowntop, lower_range+eps, upper_range-eps)


			binning = set_dyn_binning(evdyn, lower_range, upper_range, nbin)
			hdyn 		= rplot.Hist(binning);map(hdyn.Fill, evdyn)
			htop 		= rplot.Hist(binning);map(htop.Fill, evtop)
			htwotop		= rplot.Hist(binning);map(htwotop.Fill, evtwotop)
			hhalftop	= rplot.Hist(binning);map(hhalftop.Fill, evhalftop)
			htripeltop	= rplot.Hist(binning);map(htripeltop.Fill, evtripeltop)
			howntop		= rplot.Hist(binning);map(howntop.Fill, evowntop)


			gs = gridspec.GridSpec(2, 1, width_ratios=[1],height_ratios=[3,1]) 
			f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
			ax1 = plt.subplot(gs[0])
			ax2 = plt.subplot(gs[1])
			gs.update(wspace=0.03)

			hdyn.Sumw2(True);htop.Sumw2(True);htwotop.Sumw2(True);hhalftop.Sumw2(True);howntop.Sumw2(True)


			hdyn.Scale(1/(hdyn.Integral(0,hdyn.GetNbinsX()+1)),"width")
			hdyn.linecolor = "black";hdyn.linewidth = 1
			rplt.hist(hdyn,fmt="none",axes=ax1,lw=0.4,color="black",label="Dynamic scale")
			rplt.errorbar(hdyn,fmt="none",axes=ax1,lw=0.4,color="black",label="_nolegend_")

			htop.Scale(1/(htop.Integral(0,htop.GetNbinsX()+1)),"width")
			htop.linecolor = "green";htop.linewidth = 1
			rplt.hist(htop,fmt="none",axes=ax1,lw=0.4,color="green",label="Top Mass")
			rplt.errorbar(htop,fmt="none",axes=ax1,lw=0.4,color="green",label="_nolegend_")

			htwotop.Scale(1/(htwotop.Integral(0,htwotop.GetNbinsX()+1)),"width")
			htwotop.linecolor = "red";htwotop.linewidth = 1
			rplt.hist(htwotop,fmt="none",axes=ax1,lw=0.4,color="red",label="Double Top Mass")
			rplt.errorbar(htwotop,fmt="none",axes=ax1,lw=0.4,color="red",label="_nolegend_")

			hhalftop.Scale(1/(hhalftop.Integral(0,hhalftop.GetNbinsX()+1)),"width")
			hhalftop.linecolor = "blue";hhalftop.linewidth = 1
			rplt.hist(hhalftop,fmt="none",axes=ax1,lw=0.4,color="blue",label="Half Top Mass")
			rplt.errorbar(hhalftop,fmt="none",axes=ax1,lw=0.4,color="blue",label="_nolegend_")

			#htripeltop.Scale(1/(htripeltop.Integral(0,htripeltop.GetNbinsX()+1)),"width")
			#htripeltop.linecolor = "magenta";htripeltop.linewidth = 1
			#rplt.hist(htripeltop,fmt="none",axes=ax1,lw=0.4,color="magenta",label="Tripel Top Mass")
			#rplt.errorbar(htripeltop,fmt="none",axes=ax1,lw=0.4,color="magenta",label="_nolegend_")

			howntop.Scale(1/(howntop.Integral(0,howntop.GetNbinsX()+1)),"width")
			howntop.linecolor = "gray";howntop.linewidth = 1
			#rplt.hist(howntop,fmt="none",axes=ax1,lw=0.4,color="gray",label="Own Scale (top mass)")
			#rplt.errorbar(howntop,fmt="none",axes=ax1,lw=0.4,color="gray",label="_nolegend_")

			#print(p1 + " " + p2 + " " + va + " KS Test")
			ksval = hdyn.KolmogorovTest(htop)
			#print(hi.KolmogorovTest(hdp))
			#print()
			f = open("scalevarexp/scale_"+p1+"_"+p2+"_"+va+"_truth.txt","w")

			for i in range(0,htop.GetNbinsX()+2):
				f.write(str(htop.GetBinLowEdge(i))+"\n")
				if htop.GetBinContent(i) == 0:
					f.write(str(0)+" ")
					f.write(str(0)+"\n")
					continue
				f.write(str(htwotop.GetBinContent(i)/htop.GetBinContent(i))+" ")
				f.write(str(hhalftop.GetBinContent(i)/htop.GetBinContent(i))+"\n")
			f.close()

			ax2.hist(binning[1:]-np.diff(binning)/2,bins=binning,weights=np.ones_like(binning[1:]),histtype="step",lw=0.8,color="black")
			htop.Divide(hdyn);htwotop.Divide(hdyn);hhalftop.Divide(hdyn);#htripeltop.Divide(hdyn);
			howntop.Divide(hdyn)
		
			rplt.errorbar(htop,fmt="none",axes=ax2,lw=0.6,color="green",label="_nolegend_")
			rplt.errorbar(htwotop,fmt="none",axes=ax2,lw=0.6,color="red",label="_nolegend_")
			rplt.errorbar(hhalftop,fmt="none",axes=ax2,lw=0.6,color="blue",label="_nolegend_")
			#rplt.errorbar(htripeltop,fmt="none",axes=ax2,lw=0.6,color="magenta",label="_nolegend_")
			#rplt.errorbar(howntop,fmt="none",axes=ax2,lw=0.6,color="gray",label="_nolegend_")

			ax2.set_ylabel("Ratio to dyn. scale")
			#ax2.legend(loc="best")
			ax2.grid(alpha=0.6)

			ax1.set_xlim(lower_range,upper_range)
			ax2.set_xlim(lower_range,upper_range)
			ax2.set_ylim(0.5,1.5)

			labelRMax(ax1,ax2,va,p1,p2)
			leg = ax1.legend()#title="KS Test: "+"{:.7f}".format(ksval),fontsize="small")
			#leg.get_title().set_fontsize('small')
			left = lower_range
			right = upper_range
			bottom, top = ax1.get_ylim()
			ax1.text(0.7*(left+right),0.4*(bottom+top),r"$\sqrt{s}=13\,$ TeV"+"\n"+"truth level")

			plt.savefig("plots/"+va+"_truth/"+"scalevar_"+p1+"_"+p2+"_"+va+".pdf",bbox_inches='tight')
			plt.close()


scale = [87.75,172.5,345,517.5]
xsec_pro 		=	np.array([0.004596,0.004036,0.003581,0.003371])
xsec_pro_err	=	np.array([5.507e-06,4.992e-06,5.01e-06,4.229e-06])

plt.figure(figsize=(5,3.75))

plt.errorbar(scale,xsec_pro,yerr=xsec_pro_err,label="Production",fmt="+",alpha=0.7)

plt.fill_between([0,1000], 0.003735+4.721e-06, 0.003735-4.721e-06,color='black',alpha=0.5,label="Dynamic",lw=0.1)
plt.fill_between([0,1000], 0.00403+5.538e-06, 0.00403-5.538e-06,color='gray',alpha=0.5,label="OwnScale",lw=0.1)
plt.xlim(80,530)
plt.ylabel("Crosssection in pb")
plt.xlabel(r"$\alpha_s$")
#plt.yscale("log")
#plt.xscale("log")
#plt.grid(alpha=0.7)
plt.legend(loc="best")
plt.savefig("plots/crosssec.pdf",bbox_inches='tight')



"""		hd = rplot.Hist(binning);map(hd.Fill, evd,np.ones_like(evd)*decay_fraction)
		hp = rplot.Hist(binning);map(hp.Fill, evp,np.ones_like(evp)*prod_fraction)

		hp.Scale(1/(hp.Integral(0,hp.GetNbinsX()+1)))
		hp.linecolor = "blue";hp.linewidth = 1
		rplt.hist(hp,fmt="none",lw=0.4,color="blue",label="production sample")
		rplt.errorbar(hp,fmt="none",lw=0.4,color="blue",label="_nolegend_")

		hd.Scale(1/(hd.Integral(0,hd.GetNbinsX()+1)))
		hd.linecolor = "red";hd.linewidth = 1
		rplt.hist(hd,fmt="none",lw=0.4,color="red",label="decay sample")
		rplt.errorbar(hd,fmt="none",lw=0.4,color="red",label="_nolegend_")

		plt.xlim(lower_range,upper_range)
		labelkin(var,par)

		plt.savefig("plots/"+par+"_truth/"+"decpro"+"_"+par+"_"+var+".pdf",bbox_inches='tight')
		plt.close()

		hp.Scale(1/(hp.Integral(0,hp.GetNbinsX()+1)))
		hp.linecolor = "blue";hp.linewidth = 1
		rplt.hist(hp,fmt="none",lw=0.4,color="blue",label="production sample")
		rplt.errorbar(hp,fmt="none",lw=0.4,color="blue",label="_nolegend_")

		hd.Scale(1/(hd.Integral(0,hd.GetNbinsX()+1)))
		hd.linecolor = "red";hd.linewidth = 1
		rplt.hist(hd,fmt="none",lw=0.4,color="red",label="decay sample")
		rplt.errorbar(hd,fmt="none",lw=0.4,color="red",label="_nolegend_")

		plt.xlim(lower_range,upper_range)
		labelRM(va,p1,p2)
			
		plt.savefig("plots/"+va+"_truth/"+"decpro"+p1+"_"+p2+"_"+va+".pdf",bbox_inches='tight')
		plt.close()
		"""