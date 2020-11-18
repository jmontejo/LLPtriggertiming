from math import exp, sin, cos
from ROOT import *
gROOT.SetBatch(1)
import array
import os, sys

gROOT.Macro('$ROOTCOREDIR/scripts/load_packages.C')
# Initialize the xAOD infrastructure
xAOD.Init()

sample = 'higgsino_150'
samples = {
    'higgsino_200':"/afs/cern.ch/work/j/jmontejo/LLPtrigger/SignalProduction/n3n4_1TeVsquark_200gev.DAOD_TRUTH3.root",
    'higgsino_150':"/afs/cern.ch/work/j/jmontejo/LLPtrigger/SignalProduction/n3n4_1TeVsquark_150gev.DAOD_TRUTH3.root",
    'higgsino_150_pythiadecay':"/afs/cern.ch/work/j/jmontejo/LLPtrigger/SignalProduction/n3n4_1TeVsquark_150gev_pythiadecay.DAOD_TRUTH3.root",
}

rfile = TFile.Open(samples[sample])
tree  = xAOD.MakeTransientTree(rfile, "CollectionTree")
recopt = xAOD.Jet.Decorator('float')('recopt')
l1pt = xAOD.Jet.Decorator('float')('l1pt')
delayfactor = xAOD.Jet.Decorator('float')('delayfactor')
Rdisplacementfactor = xAOD.Jet.Decorator('float')('Rdisplacementfactor')
Zdisplacementfactor = xAOD.Jet.Decorator('float')('Zdisplacementfactor')
isN1decay = xAOD.Jet.Decorator('int')('isN1decay')

rnd = TRandom3(123)
truth2reco_res = TF1("truth2reco_res","[0]/x+[1]",0,1000)
truth2reco_res.SetParameters(7.06285e+00, 2.74922e-02)
reco2l1_res = TF1("reco2l1_res","[0]/x+[1]",0,1000)
reco2l1_res.SetParameters(1.79319e+01,  4.24584e-02)
nbins = 25
gev = 1e-3
maxht = 800
maxjetpt = 500
time_cutpoint = 2
time_cutpoint_high = 12
Rspace_cutpoint_high = 1000 #mm end of TRT
Zspace_cutpoint_high = 300 #???
lifetimes = sorted([pow(10,i) for i in range(-1,3)]+[3*pow(10,i) for i in range(-1,3)])
cuts = ('nominal','delayed_highpt','delayed_ISR','delayed_PU')

def truth2reco(truthpt):
    recopt = truthpt - 4 #in GeV
    recopt *= (1+truth2reco_res.Eval(truthpt)*rnd.Gaus())
    return max(1,recopt)

def reco2hlt(recopt):
    hltpt = recopt #in GeV
    hltpt *= (1+0.02)*rnd.Gaus()
    return hltpt

def reco2l1(recopt):
    l1pt = (recopt - 2.15239e+01)/1.38092e+00
    l1pt *= (1+reco2l1_res.Eval(recopt)*rnd.Gaus())
    return l1pt

def getLast(p):
    for c in range(p.nChildren()):
        if p.child(c).pdgId() == p.pdgId():
            return getLast(p.child(c))
    return p

def getDecays(tree):
    toret = []
    for bsm in tree.TruthBSMWithDecayParticles:
      if bsm.status()==62:
        assert bsm.nChildren()==3, bsm.nChildren()
        child1 = bsm.child(0)
        child2 = bsm.child(1)
        child3 = bsm.child(2)
        toret.append([bsm, child1,child2,child3])
        if len(toret)==2: break
    return toret

def decorateJet(jet, dec1, dec2, todebug=None):
    recopt.set(jet, truth2reco(jet.pt()*gev) )
    l1pt.set(jet, reco2l1(recopt(jet)) )
    flav = jet.auxdataConst['int']('TrueFlavor') 
    dr1, dr2 = 9, 9
    ppt = 0
    for p in dec1:
        if p.pdgId() == flav: 
            dr1 = p.p4().DeltaR(jet.p4())
            ppt = p.pt()
            break
    for p in dec2:
        if p.pdgId() == flav: 
            dr2 = p.p4().DeltaR(jet.p4())
            ppt = p.pt()
            break
    n1 = None
    if dr1 < dr2 and dr1<0.6:
        n1 = dec1[0]
    elif dr2 < dr1 and dr2<0.6:
        n1 = dec2[0]
    elif min(dr1,dr2) < 9 and min(dr1,dr2>0.6):
        if todebug.count(flav) <=2:
            #print "No good dR matching",dr1,dr2,flav,jet.pt(),jet.eta(),jet.phi(),[(p.pdgId(), p.pt(), p.eta(),p.phi()) for p in dec1+dec2 if p.pdgId()==flav]
            #print "Bad ratio",jet.pt()/ppt
            #print todebug
            return False
    if n1:
        n14v = n1.p4()
        delayfactor.set(jet, n14v.Gamma()*(1-n14v.Beta()))
        Rdisplacementfactor.set(jet, n14v.Beta()* n14v.Gamma()*300*sin(n14v.Theta()) )
        Zdisplacementfactor.set(jet, n14v.Beta()* n14v.Gamma()*300*abs(cos(n14v.Theta()) ))
        #print "N1 pt, beta*gamma: ",n1.pt(), n14v.Beta()* n14v.Gamma()
        isN1decay.set(jet, 1 )
    else:
        delayfactor.set(jet, 0 )
        Rdisplacementfactor.set(jet, 0 )
        Zdisplacementfactor.set(jet, 0 )
        isN1decay.set(jet, 0 )
    return True

#iterate over events, compute reco and L1 pt, flag delayed jets, compute N1 beta gamma
# - for a given signal time hypothesis, can weight events by the probability that the jets pass the timing cut
# -  jet time -> LLP time -> LLP beta gamma * proper time hypo
# - need to match for each jet the exact LLP it comes from, requires new signal where one decays to -1,-2,-3
#check if events passes:
# - nominal trigger cut (reco pT > 4xx gev)
# - nominal L1 + displaced HLT trigger cut A (jet pT > 2xx gev and jet time > 2 sigma) high pT delayed
# - nominal L1 + displaced HLT trigger cut B (jet pT > 2xx gev + jet pT > 40 and jet time > 2 sigma) ISR + delayed
# - pileup L1 + displaced HLT trigger cut (pass L1 and jet time > 2 sigma) pass L1 = weight event by the prob that it passes L1 adding to PU

def jetDelayProbability(lifetime, delayfactor):
    if delayfactor==0: return 0 #if not delay cannot pass
    #integral of a normalized exponential with tau*delayfactor, from cutpoint to cutpoint_high
    return exp(- time_cutpoint/(lifetime*delayfactor)) -exp(- time_cutpoint_high/(lifetime*delayfactor))

def jetDisplacementProbability(lifetime, jet):
    if Rdisplacementfactor(jet)==0: return 1 #if not displacement always pass
    #integral of a normalized exponential with tau*delayfactor, from cutpoint to cutpoint_high
    Rdecay = lifetime*Rdisplacementfactor(jet)
    Zdecay = lifetime*Zdisplacementfactor(jet)
    return 1-exp(- Rspace_cutpoint_high/Rdecay)*exp(- Zspace_cutpoint_high/Zdecay)

def passNominal(tree, lifetime):
    nonepass  = 1 #will do 1-probability that none pass
    for jet in tree.AntiKt4TruthDressedWZJets:
        if recopt(jet) < 450: continue
        nonepass *= 1- jetDisplacementProbability(lifetime, jet)
    return 1 - nonepass

def passDelayedHighPt(tree, lifetime):
    nonepass  = 1 #will do 1-probability that none pass
    for jet in tree.AntiKt4TruthDressedWZJets:
        if recopt(jet) < 220: continue
        nonepass *= 1- jetDelayProbability(lifetime, delayfactor(jet))*jetDisplacementProbability(lifetime, jet)
    return 1 - nonepass

def passDelayedISR(tree, lifetime):
    if truth2reco(tree.AntiKt4TruthDressedWZJets.at(0).pt()*gev) < 220: return 0
    maxprob = 0
    for i1,jet1 in enumerate(tree.AntiKt4TruthDressedWZJets):
      if truth2reco(jet1.pt()*gev) < 220: continue
      nonepass  = 1 #will do 1-probability that none pass
      for i2,jet2 in enumerate(tree.AntiKt4TruthDressedWZJets):
        if i1==i2: continue
        if truth2reco(jet2.pt()*gev) < 40: continue
        nonepass *= 1- jetDelayProbability(lifetime, delayfactor(jet) )*jetDisplacementProbability(lifetime, jet)
      maxprob = max(maxprob, 1-nonepass)
    return maxprob

def countL1J(tree, l1ptcut):
    count = 0
    for jet in tree.AntiKt4TruthDressedWZJets:
        count += (l1pt(jet) > l1ptcut)
    return count

#Solve[{Integrate[a Exp[-a x], {x, 190, Infinity}, Assumptions -> a > 0] == 2.3/28550, a > 0} , a, Reals]
#func_ht21 = TF1("ht21","expo(0)")
#ht21_lambda = 0.0496132
#func_ht21.SetParameters(ht21_lambda,-ht21_lambda)
##Currently using 3.1 eta
#func_ht31 = TF1("ht31","expo(0)")
#ht31_lambda = 0.0556166
#func_ht31.SetParameters(ht31_lambda,-ht31_lambda)
ht21cut = 190
ht31cut = 150
def getHTweight(tree, htcut):
    ht21, ht31 = 0,0
    for jet in tree.AntiKt4TruthDressedWZJets:
        if l1pt(jet) > 15 and abs( jet.eta() ) < 2.1: 
            ht21 += l1pt(jet)
        if l1pt(jet) > 20:
            ht31 += l1pt(jet)
    if ht21 >= ht21cut or ht31 >= ht31cut: return 1
    return 0

def passDelayedPU(tree, lifetime):
    nonepass  = 1 #will do 1-probability that none pass
    weight_J100_pu = 3.7/28550
    weight_4J15_pu = 4.8/28550
    weight_HT21_pu   = 2.3/28550
    weight_HT31_pu   = 6.8/28550
    weight_3J15_pu = 27.1/28550
    weight_2J15_pu = 170./28550 #guess
    weight_J15_pu  = 1101./28550
    countJ15  = countL1J(tree,15)
    countJ100 = countL1J(tree,100)

    weight_J100, weight_4J15 = 0,0
    if countJ100: weight_J100 = 1
    if countJ15>=4:   weight_4J15 = 1
    elif countJ15==3: weight_4J15 = weight_J15_pu
    elif countJ15==2: weight_4J15 = weight_2J15_pu
    elif countJ15==1: weight_4J15 = weight_3J15_pu

    weight_HT = getHTweight(tree, 150)

    #print (weight_HT,weight_4J15,weight_J100,weight_HT_pu+weight_4J15_pu+weight_J100_pu)
    l1weight = max(weight_HT,weight_4J15,weight_J100,weight_HT31_pu+weight_4J15_pu+weight_J100_pu)
    #l1weight = max(weight_4J15,weight_J100,weight_HT31_pu+weight_4J15_pu+weight_J100_pu)

    nonepass  = 1 #will do 1-probability that none pass
    for jet in tree.AntiKt4TruthDressedWZJets:
      if recopt(jet) < 40: continue
      nonepass *= 1- jetDelayProbability(lifetime, delayfactor(jet) )*jetDisplacementProbability(lifetime, jet)
    hltprob = 1-nonepass
    return hltprob*l1weight

def decorateEvent(tree):
    tree.pass_nominal = {}
    tree.pass_delayed_highpt = {}
    tree.pass_delayed_ISR = {}
    tree.pass_delayed_PU = {}

    htjetpt = sum([recopt(jet) for jet in tree.AntiKt4TruthDressedWZJets ])
    leadjetpt = max([recopt(jet) for jet in tree.AntiKt4TruthDressedWZJets ])
    leadN1jetpt = max([recopt(jet) for jet in tree.AntiKt4TruthDressedWZJets if isN1decay(jet) ])
    h_htjet.Fill(htjetpt)
    h_leadjet.Fill(leadjetpt)
    h_leadN1jet.Fill(leadN1jetpt)

    for lifetime in lifetimes:
        tree.pass_nominal[lifetime] = passNominal(tree, lifetime)

        tree.pass_delayed_highpt[lifetime] = passDelayedHighPt(tree, lifetime)
        tree.pass_delayed_ISR[lifetime] = passDelayedISR(tree, lifetime)
        tree.pass_delayed_PU[lifetime] = passDelayedPU(tree, lifetime)

        h_passleadjet['nominal'][lifetime].Fill(leadjetpt, tree.pass_nominal[lifetime])
        h_passleadjet['delayed_highpt'][lifetime].Fill(leadjetpt, tree.pass_delayed_highpt[lifetime])
        h_passleadjet['delayed_ISR'][lifetime].Fill(leadjetpt, tree.pass_delayed_ISR[lifetime])
        h_passleadjet['delayed_PU'][lifetime].Fill(leadjetpt, tree.pass_delayed_PU[lifetime])
        h_passleadN1jet['nominal'][lifetime].Fill(leadN1jetpt, tree.pass_nominal[lifetime])
        h_passleadN1jet['delayed_highpt'][lifetime].Fill(leadN1jetpt, tree.pass_delayed_highpt[lifetime])
        h_passleadN1jet['delayed_ISR'][lifetime].Fill(leadN1jetpt, tree.pass_delayed_ISR[lifetime])
        h_passleadN1jet['delayed_PU'][lifetime].Fill(leadN1jetpt, tree.pass_delayed_PU[lifetime])
        h_passhtjet['nominal'][lifetime].Fill(htjetpt, tree.pass_nominal[lifetime])
        h_passhtjet['delayed_highpt'][lifetime].Fill(htjetpt, tree.pass_delayed_highpt[lifetime])
        h_passhtjet['delayed_ISR'][lifetime].Fill(htjetpt, tree.pass_delayed_ISR[lifetime])
        h_passhtjet['delayed_PU'][lifetime].Fill(htjetpt, tree.pass_delayed_PU[lifetime])

        h_lifetime.Fill(lifetime)
        h_passlifetime['nominal'].Fill(lifetime, tree.pass_nominal[lifetime])
        h_passlifetime['delayed_highpt'].Fill(lifetime, tree.pass_delayed_highpt[lifetime])
        h_passlifetime['delayed_ISR'].Fill(lifetime, tree.pass_delayed_ISR[lifetime])
        h_passlifetime['delayed_PU'].Fill(lifetime, tree.pass_delayed_PU[lifetime])

h_lifetime = TH1F('h_lifetime','h_lifetime',len(lifetimes),array.array('f',[0]+[l*2 for l in lifetimes]))
h_htjet = TH1F('h_htjet','h_htjet',nbins,0,maxht)
h_leadjet = TH1F('h_leadjet','h_leadjet',nbins,0,maxjetpt)
h_leadN1jet = TH1F('h_leadN1jet','h_leadN1jet',nbins,0,maxjetpt)
h_passlifetime = {}
h_passhtjet = {}
h_passleadjet = {}
h_passleadN1jet = {}
for cut in cuts:
    h_passlifetime[cut] =  TH1F('h_lifetime'+cut,'h_lifetime'+cut,len(lifetimes),array.array('f',[0]+[l*2 for l in lifetimes]))
    h_passhtjet[cut] = {}
    h_passleadjet[cut] = {}
    h_passleadN1jet[cut] = {}
    for lifetime in lifetimes:
        h_passhtjet[cut][lifetime] = TH1F('h_htjet_%s_%f'%(cut,lifetime),'h_htjet',nbins,0,maxht)
        h_passleadjet[cut][lifetime] = TH1F('h_leadjet_%s_%f'%(cut,lifetime),'h_leadjet',nbins,0,maxjetpt)
        h_passleadN1jet[cut][lifetime] = TH1F('h_leadN1jet_%s_%f'%(cut,lifetime),'h_leadN1jet',nbins,0,maxjetpt)

for i in range(tree.GetEntries()):
    if i%1000==0: print i
    tree.GetEntry(i)
    dec1, dec2 = getDecays(tree)
    pdgs = [jet.auxdataConst['int']('TrueFlavor') for jet in tree.AntiKt4TruthDressedWZJets if jet.auxdataConst['int']('TrueFlavor')<4]
    for jet in tree.AntiKt4TruthDressedWZJets:
        ok = decorateJet(jet,dec1,dec2,pdgs)
        if not ok: break
    if not ok: continue #keep only events where partons are properly matched
    decorateEvent(tree) 

os.system('mkdir -p plots/%s'%sample)
canv = TCanvas()
line = TLine(0,1,maxjetpt,1)
effs = {}
for i,cut in enumerate(cuts):
    print cut
    effs[cut] = TEfficiency(h_passlifetime[cut], h_lifetime)
    effs[cut].SetMarkerStyle(20)
    effs[cut].SetLineColor(i+1)
    effs[cut].SetMarkerColor(i+1)
    effs[cut].Draw("ALP" if i==0 else "same LP")
    gPad.Update()
line.Draw()
canv.SetLogx()
gPad.Update()
effs[cuts[0]].GetPaintedGraph().GetYaxis().SetRangeUser(0,1.2)
canv.SaveAs("plots/%s/eff_lifetime.png"%sample)
canv.SetLogx(0)

norm = nbins/4.
for lifetime in lifetimes:
    effs_ht = {}
    effs_lead = {}
    effs_leadN1 = {}

    h = h_htjet.DrawNormalized("",norm)
    h.GetYaxis().SetRangeUser(0,1.2)
    for i,cut in enumerate(cuts):
        effs_ht[cut] = TEfficiency(h_passhtjet[cut][lifetime], h_htjet)
        effs_ht[cut].SetMarkerStyle(20)
        effs_ht[cut].SetLineColor(i+1)
        effs_ht[cut].SetMarkerColor(i+1)
        effs_ht[cut].Draw("same LP")
    line.Draw()
    canv.SaveAs("plots/%s/eff_ht_%f.png"%(sample,lifetime))

    h = h_leadjet.DrawNormalized("",norm)
    h.GetYaxis().SetRangeUser(0,1.2)
    for i,cut in enumerate(cuts):
        effs_lead[cut] = TEfficiency(h_passleadjet[cut][lifetime], h_leadjet)
        effs_lead[cut].SetMarkerStyle(20)
        effs_lead[cut].SetLineColor(i+1)
        effs_lead[cut].SetMarkerColor(i+1)
        effs_lead[cut].Draw("same LP")
    line.Draw()
    canv.SaveAs("plots/%s/eff_lead_%f.png"%(sample,lifetime))

    h = h_leadN1jet.DrawNormalized("",norm)
    h.GetYaxis().SetRangeUser(0,1.2)
    for i,cut in enumerate(cuts):
        effs_leadN1[cut] = TEfficiency(h_passleadN1jet[cut][lifetime], h_leadN1jet)
        effs_leadN1[cut].SetMarkerStyle(20)
        effs_leadN1[cut].SetLineColor(i+1)
        effs_leadN1[cut].SetMarkerColor(i+1)
        effs_leadN1[cut].Draw("same LP")
    line.Draw()
    canv.SaveAs("plots/%s/eff_leadN1_%f.png"%(sample,lifetime))
