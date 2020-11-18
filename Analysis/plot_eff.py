from coolPlot import coolPlot
from ROOT import *
import collections

gROOT.LoadMacro("convertL1JetNoPUcorr.C+")

Variable = collections.namedtuple('Variable', 'name thevar binning L1jet L1threshold HLTjet HLTthreshold')
Efficiency = collections.namedtuple('Efficiency', 'name thevar thecut binning')

sample = ("n3n4_1TeVsquark_showerLHE","20k","Event.Weight")
#sample = ("minbias","10k","1")
#sample = ("minbiasnoelastic","10k","1")

withpufile = TFile.Open("../Delphes-3.4.2/LLPtrigger_samples/%s_withPU55.%s.root"%(sample[0],sample[1]))
withputree = withpufile.Get("Delphes")
zeropufile = TFile.Open("../Delphes-3.4.2/LLPtrigger_samples/%s_zeroPU.%s.root"%(sample[0],sample[1]))
zeroputree = zeropufile.Get("Delphes")

withputree.SetAlias("rng","sin(2*TMath::Pi()*rndm)*sqrt(-2*log(rndm))")
zeroputree.SetAlias("rng","sin(2*TMath::Pi()*rndm)*sqrt(-2*log(rndm))")
theline = "JetNoPUcorr.PT[0],JetNoPUcorr.Eta[0],rng," + ",".join(["Alt$(JetNoPUcorr.PT[%d],0),Alt$(JetNoPUcorr.Eta[%d],0),rng"%(i,i) for i in range(1,8)])
withputree.SetAlias("L1jet0pt","convertL1JetNoPUcorr(0,%s)+0"%theline)
zeroputree.SetAlias("L1jet0pt","convertL1JetNoPUcorr(0,%s)+0"%theline)
withputree.SetAlias("L1jet3pt","convertL1JetNoPUcorr(3,%s)+0"%theline)
zeroputree.SetAlias("L1jet3pt","convertL1JetNoPUcorr(3,%s)+0"%theline)

pass_J100 = "L1jet0pt>100"
pass_J75  = "L1jet0pt>75"
pass_J50  = "L1jet0pt>50"
pass_J20  = "L1jet0pt>20"
pass_J15  = "L1jet0pt>15"
pass_4J15 = "L1jet3pt>15"

variables = [
    Variable("jet_pt1","JetNoPUcorr.PT[0]",(30,0,300),1,200,1,450),
    Variable("jet_pt4","JetNoPUcorr.PT[3]",(20,0,200),4,60,4,120),
    Variable("jet_pt6","JetNoPUcorr.PT[5]",(20,0,100),4,60,6,60),
    Variable("pass_J100",pass_J100,(2,-0.5,1.5),4,60,6,60),
    Variable("pass_J15",pass_J15,(2,-0.5,1.5),4,60,6,60),
    Variable("pass_4J15",pass_4J15,(2,-0.5,1.5),4,60,6,60),
    ]

weight=sample[2]
#for var in variables:
#    withpuname = "withpu"+var.name
#    withpu = TH1F(withpuname,withpuname,*var.binning)
#    withputree.Draw("Alt$(%s,0) >> %s"%(var.thevar,withpuname),weight)
#    
#    c = coolPlot(sample[0]+"_"+var.name,[nopu,withpu,zeropu])

efficiencies = [
    Efficiency("eff_LLPjet_pt1_vs_L1jet_pt1","MaxIf$(JetNoPUcorr.PT,JetNoPUcorr.Flavor>0)",pass_J100,(30,0,300)),
    #Efficiency("eff_jet_pt1_vs_L1jet_pt1","JetNoPUcorr.PT[0]",pass_J15,(20,0,100)),
    Efficiency("eff_jet_pt1_vs_L1jet_J100_pt1","JetNoPUcorr.PT[0]",pass_J100,(20,100,300)),
    Efficiency("eff_jet_pt1_vs_L1jet_J15_pt1","JetNoPUcorr.PT[0]",pass_J15,(20,0,100)),
    Efficiency("eff_jet_pt4_vs_L1jet_pt4","JetNoPUcorr.PT[3]",pass_4J15,(20,0,100)),
    ]

for eff in efficiencies:
    withpuname = "withpu"+eff.name
    withpu = TH2F(withpuname, withpuname,2,-0.5,1.5,*eff.binning)
    withpueff = TEfficiency(withpuname+"eff", withpuname+"eff",*eff.binning)
    withputree.Draw("%s:%s >> %s"%(eff.thevar,eff.thecut,withpuname),"(%s)"%(weight))
    hwithpupass = withpu.ProjectionY("pass",2,-1)
    hwithputotal = withpu.ProjectionY("total",1,-1)
    withpueff.SetTotalHistogram(hwithputotal,"")
    withpueff.SetPassedHistogram(hwithpupass,"")

    zeropuname = "zeropu"+eff.name
    zeropu = TH2F(zeropuname, zeropuname,2,-0.5,1.5,*eff.binning)
    zeropueff = TEfficiency(zeropuname+"eff", zeropuname+"eff",*eff.binning)
    zeroputree.Draw("%s:%s >> %s"%(eff.thevar,eff.thecut,zeropuname),"(%s)"%(weight))
    hzeropupass = zeropu.ProjectionY("pass",2,-1)
    hzeroputotal = zeropu.ProjectionY("total",1,-1)
    zeropueff.SetTotalHistogram(hzeroputotal,"")
    zeropueff.SetPassedHistogram(hzeropupass,"")

    coolPlot(sample[0]+"_test_"+eff.name,[TH1F("dummy","dummy",*eff.binning),withpueff,zeropueff],overflow=False,yrange=(0,1.2),plotratio=False)
