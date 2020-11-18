from coolPlot import coolPlot
from ROOT import *

gROOT.LoadMacro("getResiduals.C+")

#rfile = TFile.Open("../Delphes-3.4.2/LLPtrigger_samples/n3n4_1TeVsquark_150gev_showerLHE_withPU55.20k.root")
rfile = TFile.Open("../Delphes-3.4.2/LLPtrigger_samples/n3n4_1TeVsquark_showerLHE_withPU55.20k.root")
#rfile = TFile.Open("../Delphes-3.4.2/LLPtrigger_samples/n3n4_1TeVsquark_showerLHE_zeroPU.20k.root")
tree = rfile.Get("Delphes")

pujets = True
tag = "_pu" if pujets else  ""
ptoffset = 10 if pujets else 0

nbins = 20
maxpt = 200
h2d = TH2D("h2d","h2d;Truth pT;Reco pT - Truth pT",nbins,5,maxpt+5,nbins,-20+ptoffset,20+ptoffset)
canv = TCanvas()

#tree.Draw("getResiduals(GenJet.PT,GenJet.Eta,GenJet.Phi,GenJet.Mass,Alt$(JetNoPUcorr.PT[0],0),Alt$(JetNoPUcorr.Eta[0],0),Alt$(JetNoPUcorr.Phi[0],0),Alt$(JetNoPUcorr.Mass[0],0),Alt$(JetNoPUcorr.PT[1],0),Alt$(JetNoPUcorr.Eta[1],0),Alt$(JetNoPUcorr.Phi[1],0),Alt$(JetNoPUcorr.Mass[1],0),Alt$(JetNoPUcorr.PT[2],0),Alt$(JetNoPUcorr.Eta[2],0),Alt$(JetNoPUcorr.Phi[2],0),Alt$(JetNoPUcorr.Mass[2],0),Alt$(JetNoPUcorr.PT[3],0),Alt$(JetNoPUcorr.Eta[3],0),Alt$(JetNoPUcorr.Phi[3],0),Alt$(JetNoPUcorr.Mass[3],0),Alt$(JetNoPUcorr.PT[4],0),Alt$(JetNoPUcorr.Eta[4],0),Alt$(JetNoPUcorr.Phi[4],0),Alt$(JetNoPUcorr.Mass[4],0),Alt$(JetNoPUcorr.PT[5],0),Alt$(JetNoPUcorr.Eta[5],0),Alt$(JetNoPUcorr.Phi[5],0),Alt$(JetNoPUcorr.Mass[5],0),Alt$(JetNoPUcorr.PT[6],0),Alt$(JetNoPUcorr.Eta[6],0),Alt$(JetNoPUcorr.Phi[6],0),Alt$(JetNoPUcorr.Mass[6],0),Alt$(JetNoPUcorr.PT[7],0),Alt$(JetNoPUcorr.Eta[7],0),Alt$(JetNoPUcorr.Phi[7],0),Alt$(JetNoPUcorr.Mass[7],0),Alt$(JetNoPUcorr.PT[8],0),Alt$(JetNoPUcorr.Eta[8],0),Alt$(JetNoPUcorr.Phi[8],0),Alt$(JetNoPUcorr.Mass[8],0),Alt$(JetNoPUcorr.PT[9],0),Alt$(JetNoPUcorr.Eta[9],0),Alt$(JetNoPUcorr.Phi[9],0),Alt$(JetNoPUcorr.Mass[9],0) ):GenJet.PT >> h2d","Event.Weight","colz")
for j in range(10):
    offset = max(0,j-3)
    tree.Draw("getResiduals(GenJet.PT[{j}],GenJet.Eta[{j}],GenJet.Phi[{j}],GenJet.Mass[{j}],Alt$({jet}.PT[0+{off}],0),Alt$({jet}.Eta[0+{off}],0),Alt$({jet}.Phi[0+{off}],0),Alt$({jet}.Mass[0+{off}],0),Alt$({jet}.PT[1+{off}],0),Alt$({jet}.Eta[1+{off}],0),Alt$({jet}.Phi[1+{off}],0),Alt$({jet}.Mass[1+{off}],0),Alt$({jet}.PT[2+{off}],0),Alt$({jet}.Eta[2+{off}],0),Alt$({jet}.Phi[2+{off}],0),Alt$({jet}.Mass[2+{off}],0), Alt$({jet}.PT[3+{off}],0),Alt$({jet}.Eta[3+{off}],0),Alt$({jet}.Phi[3+{off}],0),Alt$({jet}.Mass[3+{off}],0), Alt$({jet}.PT[4+{off}],0),Alt$({jet}.Eta[4+{off}],0),Alt$({jet}.Phi[4+{off}],0),Alt$({jet}.Mass[4+{off}],0),Alt$({jet}.PT[5+{off}],0),Alt$({jet}.Eta[5+{off}],0),Alt$({jet}.Phi[5+{off}],0),Alt$({jet}.Mass[5+{off}],0),Alt$({jet}.PT[6+{off}],0),Alt$({jet}.Eta[6+{off}],0),Alt$({jet}.Phi[6+{off}],0),Alt$({jet}.Mass[6+{off}],0),   Alt$({jet}.PT[7+{off}],0),Alt$({jet}.Eta[7+{off}],0),Alt$({jet}.Phi[7+{off}],0),Alt$({jet}.Mass[7+{off}],0),Alt$({jet}.PT[8+{off}],0),Alt$({jet}.Eta[8+{off}],0),Alt$({jet}.Phi[8+{off}],0),Alt$({jet}.Mass[8+{off}],0),Alt$({jet}.PT[9+{off}],0),Alt$({jet}.Eta[9+{off}],0),Alt$({jet}.Phi[9+{off}],0),Alt$({jet}.Mass[9+{off}],0)):GenJet.PT[{j}] >> {add}h2d".format(jet="JetNoPUcorr" if pujets else "Jet",j=j,off=offset,add="+" if j!=0 else ""),"Event.Weight","colz")

g = TF1("g","gaus(0)",-10+ptoffset,10+ptoffset)

px = h2d.ProfileX()
h2d.FitSlicesY(g)
means = gROOT.FindObject("h2d_1")
means.Print("range")
stddv = gROOT.FindObject("h2d_2")
stddv.Print("range")

graph = TGraphErrors(nbins)
for i in range(nbins):
    graph.SetPoint(i, px.GetBinCenter(i+1), means.GetBinContent(i+1))
    graph.SetPointError(i, 0, stddv.GetBinContent(i+1))
graph.SetLineColor(2)
graph.SetMarkerColor(2)
graph.SetMarkerStyle(20)
graph.Draw("same")
px.Draw("same")

canv.SaveAs("h2d%s.png"%tag)

stddv.Draw("E")
canv.SaveAs("stddv%spng"%tag)
for i in range(nbins):
    stddv.SetBinContent(i+1, stddv.GetBinContent(i+1)/stddv.GetBinCenter(i+1) )
    stddv.SetBinError(i+1, stddv.GetBinError(i+1)/stddv.GetBinCenter(i+1) )
myf = TF1("myf","[0]/x+[1]",0,maxpt)
stddv.Fit("myf")
stddv.Draw("E")
canv.SaveAs("relstddv%s.png"%tag)
