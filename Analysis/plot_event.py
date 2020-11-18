from coolPlot import coolPlot
from ROOT import *

rfile = TFile.Open("../Delphes-3.4.2/n3n4_1TeVsquark_showerLHE_150gev_zeroPU.32.myflavorlast.root")
tree = rfile.Get("Delphes")

toplot = {}
tree.GetEntry(0)

key = "GenJet"
toplot[key] = []
for j in range(tree.GenJet_size):
    jet = tree.GenJet[j]
    print vars(jet)
    print jet
    toplot[key].append((tree.GenJet[j].PT, tree.GenJet.Eta[j], tree.GenJet.Phi[j], tree.GenJet.Flavor[j]))

key = "Jet"
toplot[key] = []
for j in range(tree.Jet_size):
    toplot[key].append((tree.Jet.PT[j], tree.Jet.Eta[j], tree.Jet.Phi[j], tree.Jet.Flavor[j]))

key = "ParticleIni"
toplot[key] = []
for j in range(tree.ParticleAll_size):
    if tree.ParticleAll.Status != 23: continue
    toplot[key].append((tree.ParticleAll.PT[j], tree.ParticleAll.Eta[j], tree.ParticleAll.Phi[j], tree.ParticleAll.Flavor[j]))

print toplot
