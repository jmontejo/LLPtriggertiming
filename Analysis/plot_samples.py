from coolPlot import coolPlot
import collections
from ROOT import *
import gc
gc.disable()

gRandom.SetSeed(123)
Variable = collections.namedtuple('Variable', 'name thevar binning L1jet L1threshold HLTjet HLTthreshold')

sampleslist = [
    ("dijet_withPUnoelastic.10k","Event.Weight"),
    ("dijetpthat1_withPUnoelastic.1k","Event.Weight"),
    ("dijetpthat4bias2_withPUnoelastic.1k","Event.Weight"),
    ("minbias_withPU.10k","1"),
    ("minbiasnoelastic_withPU.10k","1"),
    ("minbiasnondiffractive_withPU.10k","1"),
    ]
sampleslist = [
    ("dijetpthat10_zeroPU.100k","Event.Weight",False),
    ("dijetpthat8_zeroPU.100k","Event.Weight",False),
#    ("dijetpthat8_zeroPU.100k","Event.Weight",True),
    ("dijetpthat6_zeroPU.10k","Event.Weight",False),
    ("dijetpthat4_zeroPU.100k","Event.Weight",False),
#    ("dijetpthat6_zeroPU.10k","Event.Weight",True),

#    ("dijetpthat8_withminbiasPU.100k","Event.Weight",False),
#    ("dijetpthat8nobias_withPU50.1k","Event.Weight",False),
#    ("dijetpthat8nobias_withPU50.1k","Event.Weight",True),
##    ("dijetpthat4nobias_withPU50.10k","Event.Weight",False),
#    ("dijetpthat4nobias_withPU50.10k","Event.Weight",True),
##    ("dijetpthat2nobias_withPU50.10k","Event.Weight",False),
#    ("dijetpthat2nobias_withPU50.10k","Event.Weight",True),
#    ("dijetpthat4bias2_withPUnoelastic.1k","Event.Weight",True),
#    ("dijet_withPUnoelastic.10k","Event.Weight",True),
#    ("dijetpthat4_zeroPU.100k","Event.Weight"),
#    ("dijetpthat1_zeroPU.100k","Event.Weight"),
#    ("dijet_zeroPU.100k","Event.Weight"),
#    ("minbias_zeroPU.100k","1"),
#    ("minbiasinelastic_zeroPU.100k","1"),
#    ("minbiasnoelastic_zeroPU.100k","1"),
#    ("minbiasnoelasticpthat1_zeroPU.100k","1"),
#    ("minbiasnoelasticpthat8_zeroPU.10k","1"),
#    ("minbiasnondiffractive_zeroPU.100k","1"),
    ]


addsf = 1

samples = collections.OrderedDict()
for s,w,pujets in sampleslist:
    f = TFile.Open("../Delphes-3.4.2/%s.root"%s)
    t = f.Get("Delphes")
    t.GetEntry(0) #WTF, I get an IO error later without this
    t.SetDirectory(0)
    samples[(s,pujets)] = (t,w)

for t,w in samples.values():
    t.SetAlias("rng","sin(2*TMath::Pi()*rndm)*sqrt(-2*log(rndm))")
    t.SetAlias("L1jet0pt","Jet.PT[0]*(1+(0.1+0.1/Jet.PT[0])*rng)") 
    t.SetAlias("pass_J100","fabs(Jet.Eta[0])<3.1 && L1jet0pt>165")
    t.SetAlias("pass_J75" ,"fabs(Jet.Eta[0])<3.1 && L1jet0pt>130")
    t.SetAlias("pass_J50" ,"fabs(Jet.Eta[0])<3.1 && L1jet0pt>90")
    t.SetAlias("pass_J20" ,"fabs(Jet.Eta[0])<3.1 && L1jet0pt>48")
pass_J    = "pass_J20+pass_J50+pass_J75+pass_J100"
pass_4J15 = "Jet.PT[3]*(1+0.15*rng)>45"
reference = "pass_J"

variables = [
    Variable("pass_J",pass_J,(10,-0.5,9.5),4,60,6,60),
    Variable("pass_4J15",pass_4J15,(2,-0.5,1.5),4,60,6,60),
    Variable("jet_pt1","Jet.PT[0]",(10,0,200),1,200,1,450),
    Variable("jet_pt4","Jet.PT[3]",(10,0,200),4,60,4,120),
    Variable("jet_pt6","Jet.PT[5]",(5,0,100),4,60,6,60),
    ]
assert variables[0].name == reference
href_J = TH1F("href_J","href_J",*variables[0].binning)
href_J.SetBinContent(2,553.)
href_J.SetBinContent(3,38.2)
href_J.SetBinContent(4,10.)
href_J.SetBinContent(5,3.7)
href_J.SetLineStyle(2)
href_J.SetLineColor(kGray)
href_4J15 = TH1F("href_4J15","href_4J15",*variables[0].binning)
href_4J15.SetBinContent(2,4.5)
href_4J15.SetLineStyle(2)
href_4J15.SetLineColor(kGray)
for var in variables:
    hists = collections.OrderedDict()
    for (sample,pujets), (tree,weight) in samples.iteritems():
        name = sample+var.name+str(pujets)
        name = name.replace(".","_")
        h = TH1F(name,name,*var.binning)
        thevar = var.thevar
        if pujets: thevar = thevar.replace("Jet","JetNoPUcorr")
        tree.Draw("Alt$(%s,0) >> %s"%(var.thevar,name),weight)
        hists[(sample,pujets)] = h

    if var.name == reference:
        if var.name == "pass_J100":
            sfs = [(3.7/(h.GetBinContent(2) or 1),  h.GetBinError(2)/(h.GetBinContent(2) or 1)) for h in hists.values()]
        if var.name == "pass_J75":                  
            sfs = [(10./(h.GetBinContent(2) or 1),  h.GetBinError(2)/(h.GetBinContent(2) or 1)) for h in hists.values()]
        if var.name == "pass_J50":
            sfs = [(38.2/(h.GetBinContent(2) or 1), h.GetBinError(2)/(h.GetBinContent(2) or 1)) for h in hists.values()]
        if var.name == "pass_J20":
            sfs = [(553./(h.GetBinContent(2) or 1), h.GetBinError(2)/(h.GetBinContent(2) or 1)) for h in hists.values()]
        if var.name == "pass_J":
            sfs = [((553.-38.2-10.-3.7)/(h.GetBinContent(2) or 1), h.GetBinError(2)/(h.GetBinContent(2) or 1)) for h in hists.values()]
        
    for ((s,pujets),h),(sf,sferr) in zip(hists.items(),sfs):
        h.Scale(sf*addsf)
        if pujets: h.SetLineStyle(2)
        if var.name == reference:
            print s,"normalized with sf and relative error:",sf,sferr
        #for b in range(h.GetNbinsX()):
        #    h.SetBinError(b+1,max(h.GetBinError(b+1), h.GetBinContent(b+1)*sferr))
        if "pass" in var.name:
            print var.name,s,h.GetBinContent(2),h.GetBinError(2)

    extra = []
    if var.name == "pass_J": extra = [href_J]
    if var.name == "pass_4J15": extra = [href_4J15]
    #c = coolPlot("samples_zeroPU_"+var.name,hists.values(), titlelist=hists.keys(),yrangeratio=(0,2))
    c = coolPlot("samples_zeroPU_"+var.name,hists.values(), titlelist=[h[0] for h in hists.keys()],yrange=((0,100) if "pass" in var.name else None),yrangeratio=(0,2),additionals=extra)
    #c = coolPlot("samples_zeroPU_norm_"+var.name,hists.values(), titlelist=hists.keys(),normalized=True,yrangeratio=(0,2))
    #c = coolPlot("samples_zeroPU_norm_"+var.name,hists.values(), yrange=(0,0.02), titlelist=hists.keys(),normalized=True,yrangeratio=(0,2))
