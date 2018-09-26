from dfunk import *

samples = ["dec","int","pro"]
vari = ["PT","Eta","Phi","M"]
part = ["TopQuark","Photon","bJet","Jet","WBoson"]

PREFIX = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+"/"

datas = os.listdir(PREFIX+"data")

ntruth = 0
for i in range(len(datas)):
	if("truth" in datas[i]):
		ntruth += 1

for par in part:
	if(par not in os.listdir(PREFIX+"plots/")):
		os.chdir(PREFIX+"plots/")
		os.mkdir(par)
		os.chdir(PREFIX)

decay_fraction = 2266./7084.
prod_fraction = 4818./7084.


ratiox=6
ratioy=4
eps=0.001

for par in part:
	for var in vari:
		if(par == "Photon" and var=="M"):
			continue
		para = parameters_kin("dec", var, par)
		lower_range=para[0]
		upper_range=para[1]
		nbin=para[2]
		up_l=para[3]



		evdyn 		= np.genfromtxt("data/"+"dyn"+"_"+par+"_"+var+".txt")
		evtop 		= np.genfromtxt("data/"+"top"+"_"+par+"_"+var+".txt")
		evtwotop	= np.genfromtxt("data/"+"ttop"+"_"+par+"_"+var+".txt")
		evhalftop	= np.genfromtxt("data/"+"htop"+"_"+par+"_"+var+".txt")
		evtripeltop	= np.genfromtxt("data/"+"tritop"+"_"+par+"_"+var+".txt")
		evowntop	= np.genfromtxt("data/"+"owntop"+"_"+par+"_"+var+".txt")

		#evi = np.genfromtxt("data/"+"int"+"_"+par+"_"+var+"_truth.txt")
		#evd = np.genfromtxt("data/"+"dec"+"_"+par+"_"+var+"_truth.txt")
		#evp = np.genfromtxt("data/"+"pro"+"_"+par+"_"+var+"_truth.txt")

		evdyn 		= evdyn[(np.abs(evdyn)>up_l) & (evdyn!=999.9)];
		evtop 		= evtop[(np.abs(evtop)>up_l) & (evtop!=999.9)];
		evtwotop	= evtwotop[(np.abs(evtwotop)>up_l) & (evtwotop!=999.9)];
		evhalftop	= evhalftop[(np.abs(evhalftop)>up_l) & (evhalftop!=999.9)]
		evtripeltop	= evtripeltop[(np.abs(evtripeltop)>up_l) & (evtripeltop!=999.9)]
		evowntop	= evowntop[(np.abs(evowntop)>up_l) & (evowntop!=999.9)]

		evdyn 		 = np.clip(evdyn, lower_range+eps, upper_range-eps);
		evtop 		 = np.clip(evtop, lower_range+eps, upper_range-eps);
		evtwotop	 = np.clip(evtwotop, lower_range+eps, upper_range-eps);
		evhalftop	 = np.clip(evhalftop, lower_range+eps, upper_range-eps)
		evtripeltop	 = np.clip(evtripeltop, lower_range+eps, upper_range-eps)
		evowntop	 = np.clip(evowntop, lower_range+eps, upper_range-eps)
	
	
		

		gs = gridspec.GridSpec(2, 1, width_ratios=[1],height_ratios=[3,1]) 
		f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
		ax1 = plt.subplot(gs[0])
		ax2 = plt.subplot(gs[1])
		#plt.subplots_adjust(wspace=-1)
		#nbin=128
		#binning = np.linspace(lower_range, upper_range, nbin)
		binning = set_dyn_binning(evdyn, lower_range, upper_range, nbin, 0.07)

		hdyn 		= rplot.Hist(binning);map(hdyn.Fill, evdyn)
		htop 		= rplot.Hist(binning);map(htop.Fill, evtop)
		htwotop		= rplot.Hist(binning);map(htwotop.Fill, evtwotop)
		hhalftop	= rplot.Hist(binning);map(hhalftop.Fill, evhalftop)
		htripeltop	= rplot.Hist(binning);map(htripeltop.Fill, evtripeltop)
		howntop		= rplot.Hist(binning);map(howntop.Fill, evowntop)

		hdyn.Sumw2(True);htop.Sumw2(True);htwotop.Sumw2(True);hhalftop.Sumw2(True);#htripeltop.Sumw2(True);
		howntop.Sumw2(True)

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
		#rplt.hist(howntop,fmt="none",axes=ax1,lw=0.4,color="gray",label="Own scale (top mass)")
		#rplt.errorbar(howntop,fmt="none",axes=ax1,lw=0.4,color="gray",label="_nolegend_")

		#print(par + " " + var + " KS Test")
		ksval = hdyn.KolmogorovTest(htop)
		#print(hi.KolmogorovTest(hdp))
		#print()

		f = open("scalevarexp/scale_"+var+"_"+par+"_delph.txt","w")

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
		htop.Divide(hdyn);htwotop.Divide(hdyn);hhalftop.Divide(hdyn);htripeltop.Divide(hdyn);howntop.Divide(hdyn)
		

		rplt.errorbar(htop,fmt="none",axes=ax2,lw=0.6,color="green",label="_nolegend_")
		rplt.errorbar(htwotop,fmt="none",axes=ax2,lw=0.6,color="red",label="_nolegend_")
		rplt.errorbar(hhalftop,fmt="none",axes=ax2,lw=0.6,color="blue",label="_nolegend_")
		#rplt.errorbar(htripeltop,fmt="none",axes=ax2,lw=0.6,color="magenta",label="_nolegend_")
		#rplt.errorbar(howntop,fmt="none",axes=ax2,lw=0.6,color="gray",label="_nolegend_")

		#hi.Divide(hi)
		#rplt.errorbar(hi,fmt="none",axes=ax2,lw=0.6,color="black",label="_nolegend_")

		ax2.set_ylabel("Ratio to dyn. scale(MG)")
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
		ax1.text(0.7*(left+right),0.4*(bottom+top),r"$\sqrt{s}=13\,$ TeV"+"\n"+"reco level")

		plt.savefig("plots/"+par+"/"+"scalevar"+"_"+par+"_"+var+".pdf",bbox_inches='tight')
		plt.close()


"""for par in part:
	for var in vari:
		if(par == "Photon" and var=="M"):
			continue

		para = parameters_kin("dec", var, par)
		lower_range=para[0]
		upper_range=para[1]
		nbin=para[2]
		up_l=para[3]

		evp = np.genfromtxt("data/pro_"+par+"_"+var+".txt")
		evd = np.genfromtxt("data/dec_"+par+"_"+var+".txt")
		evi = np.genfromtxt("data/int_"+par+"_"+var+".txt")
		evp = evp[(np.abs(evp)>up_l) & (evp!=999.9)];evd = evd[(np.abs(evd)>up_l) & (evd!=999.9)];evi = evi[(np.abs(evi)>up_l) & (evi!=999.9)]
		ev = np.concatenate((evp, evd), axis=0)
		#print(par+var)
		ev = np.clip(ev, lower_range, upper_range)
		evi = np.clip(evi, lower_range, upper_range)
		evp = np.clip(evp, lower_range, upper_range)
		evd = np.clip(evd, lower_range, upper_range)


		nbpV,binspV,apV = plt.hist(evp,bins=nbin,lw=0.5,alpha=0.0,color="blue",fill=False,weights=np.ones_like(evp)*prod_fraction,range=(lower_range,upper_range),histtype='step')
		nbdV,binsdV,adV = plt.hist(evd,bins=binspV,lw=0.5,alpha=0.0,color="blue",fill=False,weights=np.ones_like(evd)*decay_fraction,range=(lower_range,upper_range),histtype='step')
		nV,binsV,aV 	= plt.hist(binspV[1:]-np.diff(binspV)/2,label=r"production+decay",bins=binspV,lw=0.5,color="blue",fill=False,weights=nbpV+nbdV,range=(lower_range,upper_range),histtype='step')

		vnI,vbinsI,vaI = plt.hist(evi,label=r"interference",bins=nbin,lw=0.5,color="red",fill=False,normed=False,range=(lower_range,upper_range),histtype='step')
		plt.close()

		plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')

		nbp,binsp,ap = plt.hist(evp,bins=nbin,lw=0.5,alpha=0.0,color="blue",fill=False,weights=np.ones_like(evp)*prod_fraction/len(evp),range=(lower_range,upper_range),histtype='step')
		nbd,binsd,ad = plt.hist(evd,bins=binsp,lw=0.5,alpha=0.0,color="blue",fill=False,weights=np.ones_like(evd)*decay_fraction/len(evd),range=(lower_range,upper_range),histtype='step')
		n,bins,a 	 = plt.hist(binsp[1:]-np.diff(binsp)/2,label=r"production + decay mode",bins=binsp,lw=0.8,color="blue",fill=False,weights=nbp+nbd,range=(lower_range,upper_range),histtype='step')

		plot_error_region2(n,1/np.sqrt(nV)*n, bins,"blue")
		nI,binsI,aI = plt.hist(evi,label=r"interference sample",bins=bins,lw=0.8,color="red",fill=False,weights=np.ones_like(evi)/float(len(evi)),range=(lower_range,upper_range),histtype='step')
		plot_error_region2(nI,1/np.sqrt(vnI)*nI, binsI,"red")
		

		plt.xlim(lower_range,upper_range)
		labelkin(var,par)

		plt.savefig("plots/"+par+"/"+"all"+"_"+par+"_"+var+".pdf",bbox_inches='tight')
		plt.close()"""



if("R" not in os.listdir(PREFIX+"plots/")):
	os.chdir(PREFIX+"plots/")
	os.mkdir("R/")
	os.chdir(PREFIX)

if("M" not in os.listdir(PREFIX+"plots/")):
	os.chdir(PREFIX+"plots/")
	os.mkdir("M")
	os.chdir(PREFIX)

		
var=["R","M"]
RMPart=["Photon","TopQuark"]
RMPartP=["Photon","TopQuark","bJet","WBoson","LeadingJet"]



for va in var:
	for p1 in RMPart:
		for p2 in RMPartP:
			if(p1==p2):
				continue

			para = parameters_RM(p1, p2, va)
			lower_range=para[0];upper_range=para[1];nbin=para[2];up_l=para[3]


			evdyn 		= np.genfromtxt("data/"+"dyn_"+p1+"_"+p2+"_"+va+".txt")
			evtop 		= np.genfromtxt("data/"+"top_"+p1+"_"+p2+"_"+va+".txt")
			evtwotop	= np.genfromtxt("data/"+"ttop_"+p1+"_"+p2+"_"+va+".txt")
			evhalftop	= np.genfromtxt("data/"+"htop_"+p1+"_"+p2+"_"+va+".txt")
			evtripeltop	= np.genfromtxt("data/"+"tritop_"+p1+"_"+p2+"_"+va+".txt")
			evowntop	= np.genfromtxt("data/"+"owntop_"+p1+"_"+p2+"_"+va+".txt")

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


			hdyn 		= rplot.Hist(binning);map(hdyn.Fill, evdyn)
			htop 		= rplot.Hist(binning);map(htop.Fill, evtop)
			htwotop		= rplot.Hist(binning);map(htwotop.Fill, evtwotop)
			hhalftop	= rplot.Hist(binning);map(hhalftop.Fill, evhalftop)

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
			f = open("scalevarexp/scale_"+p1+"_"+p2+"_"+va+"_delph.txt","w")

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



			ax2.set_ylabel("Ratio to dyn. scale(MG)")
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
			ax1.text(0.7*(left+right),0.4*(bottom+top),r"$\sqrt{s}=13\,$ TeV"+"\n"+"reco level")
			
			plt.savefig("plots/"+va+"/"+"scalevar_"+p1+"_"+p2+"_"+va+".pdf",bbox_inches='tight')
			plt.close()




"""for va in var:
	for p1 in RMPart:
		for p2 in RMPartP:
			if(p1==p2):
				continue

			para = parameters_RM(p1, p2, va)
			lower_range=para[0];upper_range=para[1];nbin=para[2];up_l=para[3]

			evp = np.genfromtxt("data/pro_"+p1+"_"+p2+"_"+va+".txt")
			evd = np.genfromtxt("data/dec_"+p1+"_"+p2+"_"+va+".txt")
			evi = np.genfromtxt("data/int_"+p1+"_"+p2+"_"+va+".txt")
			evp = evp[(np.abs(evp)>up_l) & (evp!=999.9)];evd = evd[(np.abs(evd)>up_l) & (evd!=999.9)];evi = evi[(np.abs(evi)>up_l) & (evi!=999.9)]
			evi = np.clip(evi, lower_range, upper_range);evp = np.clip(evp, lower_range, upper_range);evd = np.clip(evd, lower_range, upper_range)
			#ev = np.concatenate((evp, evd), axis=0)
			#ev = np.clip(ev, lower_range, upper_range)

			fig = plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')


			nbpV2,binspV2,apV2 = plt.hist(evp,bins=nbin,lw=0.5,alpha=0.0,color="blue",fill=False,weights=np.ones_like(evp)*prod_fraction,range=(lower_range,upper_range),histtype='step')
			nbdV2,binsdV2,adV2 = plt.hist(evd,bins=binspV2,lw=0.5,alpha=0.0,color="blue",fill=False,weights=np.ones_like(evd)*decay_fraction,range=(lower_range,upper_range),histtype='step')
			nV2,binsV2,aV2 	= plt.hist(binspV2[1:]-np.diff(binspV2)/2,label=r"production+decay",bins=binspV2,lw=0.5,color="blue",fill=False,weights=nbpV2+nbdV2,range=(lower_range,upper_range),histtype='step')

			vnI2,vbinsI2,vaI2 = plt.hist(evi,label=r"interference",bins=nbin,lw=0.8,color="red",fill=False,normed=False,histtype='step')#,range=(lower_range,upper_range))
			
			plt.close()
			plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')
			
			nbp2,binsp2,ap2 = plt.hist(evp,bins=binspV2,lw=0.5,alpha=0.0,color="blue",fill=False,weights=np.ones_like(evp)*prod_fraction/len(evp),range=(lower_range,upper_range),histtype='step')
			nbd2,binsd2,ad2 = plt.hist(evd,bins=binspV2,lw=0.5,alpha=0.0,color="blue",fill=False,weights=np.ones_like(evd)*decay_fraction/len(evd),range=(lower_range,upper_range),histtype='step')
			n2,bins2,a2 	= plt.hist(binspV2[1:]-np.diff(binspV2)/2,label=r"production + decay mode",bins=binspV2,lw=0.8,color="blue",fill=False,weights=nbp2+nbd2,range=(lower_range,upper_range),histtype='step')
			
			nI,binsI,aI = plt.hist(evi,label=r"interference sample",bins=bins2,lw=0.8,color="red",fill=False,weights=np.ones_like(evi)/float(len(evi)),histtype='step')#,range=(lower_range,upper_range))
			
			plot_error_region2(n2,1/np.sqrt(nV2)*n2, bins2,"blue")
			plot_error_region2(nI,1/np.sqrt(vnI2)*nI, binsI,"red")

			plt.xlim(lower_range,upper_range)
			labelRM(va,p1,p2)
			plt.savefig("plots/"+va+"/"+"all_"+p1+"_"+p2+"_"+va+".pdf",bbox_inches='tight')
			plt.close()"""





		#plt.errorbar(binsp[:-1]+np.diff(binsp)/2,nbp+nbd,drawstyle = 'steps-mid',color="blue",lw=0.5)
		#plt.bar(binsp[:-1],(nbp+nbd),align="edge",label=r"production+decay",width=np.diff(binsp),lw=0.5,fill=False,edgecolor="blue")
