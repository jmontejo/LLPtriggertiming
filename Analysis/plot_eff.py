from coolPlot import coolPlot
from ROOT import *
import collections

gROOT.LoadMacro("convertL1Jet.C+")

Variable = collections.namedtuple('Variable', 'name thevar binning L1jet L1threshold HLTjet HLTthreshold')
Efficiency = collections.namedtuple('Efficiency', 'name thevar thecut binning')

samples = [
    ("n3n4_1TeVsquark_150gev_showerLHE_withPU55.20k.root","higgsino 150gev PU 55"),
    ("MinbiasWithPU55.root","minbias PU 55"),
    ("MinbiasWithPU50.root","minbias PU 50"),
]
#sample = ("minbias","10k","1")
#sample = ("minbiasnoelastic","10k","1")


pass_J100 = "L1jet0pt>100"
pass_J75  = "L1jet0pt>75"
pass_J50  = "L1jet0pt>50"
pass_J20  = "L1jet0pt>20"
pass_J15  = "L1jet0pt>15"
pass_4J15 = "L1jet3pt>15"

variables = [
    Variable("jet_pt1","Jet.PT[0]",(30,0,300),1,200,1,450),
    Variable("jet_pt4","Jet.PT[3]",(20,0,200),4,60,4,120),
    Variable("jet_pt6","Jet.PT[5]",(20,0,100),4,60,6,60),
    Variable("pass_J100",pass_J100,(2,-0.5,1.5),4,60,6,60),
    Variable("pass_J15",pass_J15,(2,-0.5,1.5),4,60,6,60),
    Variable("pass_4J15",pass_4J15,(2,-0.5,1.5),4,60,6,60),
    ]

#for var in variables:
#    withpuname = "withpu"+var.name
#    withpu = TH1F(withpuname,withpuname,*var.binning)
#    withputree.Draw("Alt$(%s,0) >> %s"%(var.thevar,withpuname),"Event.Weight")
#    
#    c = coolPlot(sample[0]+"_"+var.name,[nopu,withpu,zeropu])

efficiencies = [
    Efficiency("eff_LLPjet_pt1_vs_L1jet_pt1","MaxIf$(Jet.PT,Jet.Flavor>0)",pass_J100,(30,0,300)),
    #Efficiency("eff_jet_pt1_vs_L1jet_pt1","Jet.PT[0]",pass_J15,(20,0,100)),
    Efficiency("eff_jet_pt1_vs_L1jet_J100_pt1","Jet.PT[0]",pass_J100,(20,100,300)),
    Efficiency("eff_jet_pt1_vs_L1jet_J15_pt1","Jet.PT[0]",pass_J15,(20,0,100)),
    Efficiency("eff_jet_pt4_vs_L1jet_pt4","Jet.PT[3]",pass_4J15,(20,0,100)),
    ]
theline = "Jet.PT[0],Jet.Eta[0],rng," + ",".join(["Alt$(Jet.PT[%d],0),Alt$(Jet.Eta[%d],0),rng"%(i,i) for i in range(1,8)])

tefficiencies = {}
for sample, samplename in samples:
    rfile = TFile.Open("../LLPtrigger_samples/"+sample)
    tree = rfile.Get("Delphes")
    
    tree.SetAlias("rng","sin(2*TMath::Pi()*rndm)*sqrt(-2*log(rndm))")
    tree.SetAlias("L1jet0pt","convertL1Jet(0,%s)+0"%theline)
    tree.SetAlias("L1jet3pt","convertL1Jet(3,%s)+0"%theline)
    tefficiencies[sample]={}

    for eff in efficiencies:
        effname = sample+eff.name
        h2d = TH2F(effname, effname,2,-0.5,1.5,*eff.binning)
        teff = TEfficiency(effname+"eff", effname+"eff",*eff.binning)
        teff.SetDirectory(0)
        tree.Draw("%s:%s >> %s"%(eff.thevar,eff.thecut,effname),"Event.Weight")

        hpass = h2d.ProjectionY("pass",2,-1)
        htotal = h2d.ProjectionY("total",1,-1)
        teff.SetTotalHistogram(htotal,"")
        teff.SetPassedHistogram(hpass,"")
        tefficiencies[sample][eff.name] = teff
    rfile.Close()

print tefficiencies
for eff in efficiencies:
    coolPlot(eff.name,[TH1F("dummy","dummy",*eff.binning)]+[tefficiencies[sample[0]][eff.name] for sample in samples], titlelist=["dummy"]+[sample[1] for sample in samples], overflow=False,yrange=(0,1.2),plotratio=False)
