from coolPlot import coolPlot
import collections
from ROOT import *
import gc
gc.disable()

gROOT.LoadMacro("convertL1Jet.C+")
gRandom.SetSeed(123)
Variable = collections.namedtuple('Variable', 'name thevar binning L1jet L1threshold HLTjet HLTthreshold')

sampleslist = [
    ("MinbiasWithPU50","Event.Weight",False),
    ("MinbiasWithPU55","Event.Weight",False),
    ("MinbiasWithPU50","Event.Weight",True),
    ("MinbiasWithPU55","Event.Weight",True),
    ("n3n4_1TeVsquark_showerLHE_withPU55.20k","Event.Weight",False),
    ("n3n4_1TeVsquark_showerLHE_zeroPU.20k","Event.Weight",False),
#    ("MinbiasInelasticWithPU50","Event.Weight",False),
#    ("minbias_withPU.10k","Event.Weight",False),
#    ("minbiasnoelastic_withPU.10k","Event.Weight",False),
#    ("minbiasnondiffractive_withPU.10k","Event.Weight",False),
#    ("minbiasnondiffractive_withPU55.10k","Event.Weight",False),

    #("MinbiasInelasticWithPU50","Event.Weight",True),
    #("minbias_withPU.10k","Event.Weight",False),
    #("minbiasnoelastic_withPU.10k","Event.Weight",False),
    #("minbiasnondiffractive_withPU.10k","Event.Weight",False),
    #("minbiasnondiffractive_withPU55.10k","Event.Weight",False),
    #("dijet_withPUnoelastic.10k","Event.Weight",False),
    #("dijetpthat1_withPUnoelastic.1k","Event.Weight",False),
    #("dijetpthat4bias2_withPUnoelastic.1k","Event.Weight",False),
    #("dijetpthat8nobias_withPU50.100","Event.Weight",False),
    #("dijetpthat8nobias_withPU50_pujets.100","Event.Weight",False),
    #("dijetpthat2nobias_withPU50.10k","Event.Weight",False),
    #("dijetpthat4nobias_withPU50.10k","Event.Weight",False),
    ]

addsf = 1

pass_J    = "pass_J20+pass_J50+pass_J75+pass_J100"
reference = "pass_J50"

variables = [
    Variable("pass_J50","pass_J50",(10,-0.5,9.5),4,60,6,60),
    Variable("pass_J",pass_J,(10,-0.5,9.5),4,60,6,60),
    Variable("pass_HT21","pass_HT21",(10,-0.5,9.5),4,60,6,60),
    Variable("pass_HT31","pass_HT31",(10,-0.5,9.5),4,60,6,60),
    Variable("pass_4J15","pass_4J15",(2,-0.5,1.5),4,60,6,60),
    Variable("count_J15","count_J15",(7,-0.5,7.5),4,60,6,60),
    Variable("pass_3J40","pass_3J40",(2,-0.5,1.5),4,60,6,60),
    Variable("pass_J85_3J30","pass_J85_3J30",(2,-0.5,1.5),4,60,6,60),
    Variable("n_jet","Jet_size",(10,-0.5,9.5),1,200,1,450),
    Variable("ht21","myHT21",(10,0,200),1,200,1,450),
    Variable("ht31","myHT31",(10,0,200),1,200,1,450),
    Variable("jet_pt1","Jet.PT[0]",(10,0,200),1,200,1,450),
    Variable("jet_pt4","Jet.PT[3]",(10,0,200),4,60,4,120),
    Variable("jet_pt6","Jet.PT[5]",(5,0,100),4,60,6,60),
    ]
assert variables[0].name == reference
href_J = TH1D("href_J","href_J",10,-0.5,9.5)
href_J.SetBinContent(1,28550-553.-38.2-10.-3.7)
href_J.SetBinContent(2,553.-38.2-10.-3.7)
href_J.SetBinContent(3,38.2-10.-3.7)
href_J.SetBinContent(4,10.-3.7)
href_J.SetBinContent(5,3.7)
href_J.SetLineStyle(2)
href_J.SetLineColor(kGray)
href_4J15 = TH1D("href_4J15","href_4J15",2,-0.5,1.5)
href_4J15.SetBinContent(1,28550-4.5)
href_4J15.SetBinContent(2,4.5)
href_4J15.SetLineStyle(2)
href_4J15.SetLineColor(kGray)
href_3J40 = TH1D("href_3J40","href_3J40",2,-0.5,1.5)
href_3J40.SetBinContent(2,1.1)
href_3J40.SetLineStyle(2)
href_3J40.SetLineColor(kGray)
hists = collections.OrderedDict()
sfs = collections.OrderedDict()
theline = "Jet.PT[0],Jet.Eta[0],rng," + ",".join(["Alt$(Jet.PT[%d],0),Alt$(Jet.Eta[%d],0),rng"%(i,i) for i in range(1,8)])
thepuline = theline.replace("Jet","JetNoPUcorr")
thepuline = thepuline.replace(",JetNoPUcorr.Eta","-12,JetNoPUcorr.Eta").replace(",0),Alt$(JetNoPUcorr.Eta","-13,0),Alt$(JetNoPUcorr.Eta")
print thepuline
for (sample,weight,pujets) in sampleslist:
    f = TFile.Open("../Delphes-3.4.2/LLPtrigger_samples/%s.root"%sample)
    tree = f.Get("Delphes")
    tree.SetAlias("rng","sin(2*TMath::Pi()*rndm)*sqrt(-2*log(rndm))")
    if pujets:
        tree.SetAlias("L1jet0pt","convertL1Jet(0,%s)+0"%thepuline)
        tree.SetAlias("L1jet1pt","convertL1Jet(1,%s)+0"%thepuline)
        tree.SetAlias("L1jet2pt","convertL1Jet(2,%s)+0"%thepuline)
        tree.SetAlias("L1jet3pt","convertL1Jet(3,%s)+0"%thepuline)
        tree.SetAlias("L1jet4pt","convertL1Jet(4,%s)+0"%thepuline)
        tree.SetAlias("L1jet5pt","convertL1Jet(5,%s)+0"%thepuline)
        tree.SetAlias("L1jet6pt","convertL1Jet(6,%s)+0"%thepuline)
    else:
        tree.SetAlias("L1jet0pt","convertL1Jet(0,%s)+0"%theline)
        tree.SetAlias("L1jet1pt","convertL1Jet(1,%s)+0"%theline)
        tree.SetAlias("L1jet2pt","convertL1Jet(2,%s)+0"%theline)
        tree.SetAlias("L1jet3pt","convertL1Jet(3,%s)+0"%theline)
        tree.SetAlias("L1jet4pt","convertL1Jet(4,%s)+0"%theline)
        tree.SetAlias("L1jet5pt","convertL1Jet(6,%s)+0"%theline)
        tree.SetAlias("L1jet6pt","convertL1Jet(6,%s)+0"%theline)
    tree.SetAlias("myHT21","convertL1Jet(-1,%s,2.1)+0"%theline)
    tree.SetAlias("myHT31","convertL1Jet(-1,%s,3.1)+0"%theline)
    tree.SetAlias("pass_J100","L1jet0pt>100")
    tree.SetAlias("pass_J75" ,"L1jet0pt>75")
    tree.SetAlias("pass_J50" ,"L1jet0pt>50")
    tree.SetAlias("pass_J20" ,"L1jet0pt>20")
    tree.SetAlias("pass_4J15" ,"L1jet3pt>15")
    tree.SetAlias("count_J15" ,"+".join(["(L1jet%dpt>15)"%j for j in range(7)]))
    tree.SetAlias("pass_3J40" ,"L1jet2pt>40")
    tree.SetAlias("pass_HT21" ,"myHT21>190")
    tree.SetAlias("pass_HT31" ,"myHT31>150")
    tree.SetAlias("pass_J85_3J30" ,"L1jet0pt>85 && L1jet2pt>30")

    for var in variables:
        if not var in hists: hists[var] = collections.OrderedDict()
        name = sample+var.name+str(pujets)
        name = name.replace(".","_")
        h = TH1D(name,name,*var.binning)
        thevar = var.thevar
        #if pujets: thevar = thevar.replace("Jet","JetNoPUcorr")
        tree.Draw("%s >> %s"%(var.thevar,name),weight)
        h.SetDirectory(0)
        hists[var][(sample,pujets)] = h

        if var.name == reference:
            if var.name == "pass_4J15":
                sfs[(sample,pujets)] = (4.8/(h.GetBinContent(2) or 1),  h.GetBinError(2)/(h.GetBinContent(2) or 1))
            if var.name == "pass_J100":
                sfs[(sample,pujets)] = (3.7/(h.GetBinContent(2) or 1),  h.GetBinError(2)/(h.GetBinContent(2) or 1))
            if var.name == "pass_J75":                  
                sfs[(sample,pujets)] = (10./(h.GetBinContent(2) or 1),  h.GetBinError(2)/(h.GetBinContent(2) or 1))
            if var.name == "pass_J50":
                sfs[(sample,pujets)] = (38.2/(h.GetBinContent(2) or 1), h.GetBinError(2)/(h.GetBinContent(2) or 1))
            if var.name == "pass_J20":
                sfs[(sample,pujets)] = (553./(h.GetBinContent(2) or 1), h.GetBinError(2)/(h.GetBinContent(2) or 1))
            if var.name == "pass_J":
                sfs[(sample,pujets)] = ((553.-38.2-10.-3.7)/(h.GetBinContent(2) or 1), h.GetBinError(2)/(h.GetBinContent(2) or 1))
    f.Close()
        
for var in variables:
    for ((s,pujets),h),(sf,sferr) in zip(hists[var].items(),sfs.values()):
        h.Scale(sf*addsf)
        if pujets: h.SetLineStyle(2)
        if var.name == reference:
            print s,"normalized with sf and relative error:",sf,sferr
        #for b in range(h.GetNbinsX()):
        #    h.SetBinError(b+1,max(h.GetBinError(b+1), h.GetBinContent(b+1)*sferr))
        if "pass" in var.name:
            print var.name,s,h.GetBinContent(2),h.GetBinError(2)

    if "ht" in var.name: hists[var].values()[0].Fit("expo","","",60,200)
    extra, extratitle = [],[]
    if var.name == "pass_J": extra,extratitle = [href_J],["reference"]
    if var.name == "pass_4J15": extra,extratitle = [href_4J15],["reference"]
    #if var.name == "pass_3J40": extra,extratitle = [href_3J40],["reference"]
    #c = coolPlot("samples_withPU_"+var.name,hists.values(), titlelist=hists.keys(),yrangeratio=(0,2))
    #c = coolPlot("samples_withPU_"+var.name,hists[var].values(), titlelist=[h[0] for h in hists[var].keys()],yrange=((0,100) if "pass" in var.name else None),yrangeratio=(0,2),additionals=extra)
    c = coolPlot("samples_withPU_"+var.name,extra+hists[var].values(), titlelist=extratitle+[h[0] for h in hists[var].keys()],yrange=((0,40) if "pass" in var.name or "ht" in var.name else None),yrangeratio=(0,2))
    c = coolPlot("samples_withPU_norm_"+var.name,extra+hists[var].values(), titlelist=extratitle+[h[0] for h in hists[var].keys()],normalized=True,yrangeratio=(0,4),yrange=(0,0.01))
