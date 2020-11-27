from math import exp, sin, cos, log
from ROOT import *
gROOT.SetBatch(1)
import array
import os, sys

gROOT.Macro('$ROOTCOREDIR/scripts/load_packages.C')
# Initialize the xAOD infrastructure
xAOD.Init()

sample = 'higgsino_150'
samples = {
    'higgsino_200':"/afs/cern.ch/work/j/jmontejo/LLPtriggertiming/LLPtrigger_samples/n3n4_1TeVsquark_200gev.DAOD_TRUTH3.root",
    'higgsino_150':"/afs/cern.ch/work/j/jmontejo/LLPtriggertiming/LLPtrigger_samples/n3n4_1TeVsquark_150gev.DAOD_TRUTH3.root",
    'higgsino_150_pythiadecay':"/afs/cern.ch/work/j/jmontejo/LLPtriggertiming/LLPtrigger_samples/n3n4_1TeVsquark_150gev_pythiadecay.DAOD_TRUTH3.root",
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
#>>> ROOT.Math.normal_cdf_c(1)
#0.15865525393145707
#>>> ROOT.Math.normal_cdf_c(2)
#0.022750131948179216
#>>> ROOT.Math.normal_cdf_c(1)**2
#0.025171489600055125
#>>> ROOT.Math.normal_cdf_c(2)**2
#0.0005175685036595646
fiducial_ranges = {
    'contained_for_delayed' : ('randzup',(0,1000,300)),
    'contained_for_delayed' : ('randzup',(0,10000,3000)),
    'vtx_range1' : ('randzup',(4,40,300)),
    'vtx_range2' : ('randzup',(40,100,300)),
    'vtx_range3' : ('randzup',(100,200,300)),
    'vtx_range4' : ('randzup',(200,300,300)),
    'calratio_barrel_range1' : ('ronly',(2000,2200)),
    'calratio_barrel_range2' : ('ronly',(2200,3300)),
    'calratio_barrel_range3' : ('ronly',(3300,3700)),
    'calratio_endcap_range1' : ('zonly',(4200,4300)),
    'calratio_endcap_range2' : ('zonly',(4300,5400)),
    'calratio_endcap_range3' : ('zonly',(5400,5800)),
}
lifetimes = sorted([pow(10,i) for i in range(-2,3)]+[3*pow(10,i) for i in range(-2,3)])
cuts = ('nominal','delayed_highpt','delayed_ISR','delayed_PU','delayed2_PU','tracking_PU','CalRatio')

h_lifetime = TH1F('h_lifetime','h_lifetime',len(lifetimes),array.array('f',[0]+[l*2 for l in lifetimes]))
h_htjet = TH1F('h_htjet','h_htjet',nbins,0,maxht)
h_leadjet = TH1F('h_leadjet','h_leadjet',nbins,0,maxjetpt)
h_leadN1jet = TH1F('h_leadN1jet','h_leadN1jet',nbins,0,maxjetpt)
h_passlifetime = {}
h_passhtjet = {}
h_passleadjet = {}
h_passleadN1jet = {}

def calratioEndcapZEfficiency(llp, lifetime):
    eff = 0
    eff += 0.5*llpDisplacementProbability(lifetime, llp, fiducial_range='calratio_endcap_range1')
    eff += 1.0*llpDisplacementProbability(lifetime, llp, fiducial_range='calratio_endcap_range2')
    eff += 0.5*llpDisplacementProbability(lifetime, llp, fiducial_range='calratio_endcap_range3')
    assert eff<1, eff
    return eff

def calratioBarrelREfficiency(llp, lifetime):
    eff = 0
    eff += 0.5*llpDisplacementProbability(lifetime, llp, fiducial_range='calratio_barrel_range1')
    eff += 1.0*llpDisplacementProbability(lifetime, llp, fiducial_range='calratio_barrel_range2')
    eff += 0.5*llpDisplacementProbability(lifetime, llp, fiducial_range='calratio_barrel_range3')
    assert eff<1, eff
    return eff

def calratioPtEfficiency(LLPpt):
    if LLPpt<20: return 0
    if LLPpt>260: return 0.8
    return log((LLPpt-20)/10.)/4.

def vtxEfficiency(lifetime, llp):
    eff = 0
    eff += 0.9*llpDisplacementProbability(lifetime, llp, fiducial_range='vtx_range1')
    eff += 0.8*llpDisplacementProbability(lifetime, llp, fiducial_range='vtx_range2')
    eff += 0.7*llpDisplacementProbability(lifetime, llp, fiducial_range='vtx_range3')
    eff += 0.5*llpDisplacementProbability(lifetime, llp, fiducial_range='vtx_range4')
    assert eff<1, eff
    return eff

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
        bsm4v = bsm.p4()
        Rdisplacementfactor.set(bsm, bsm4v.Beta()* bsm4v.Gamma()*300*sin(bsm4v.Theta()) )
        Zdisplacementfactor.set(bsm, bsm4v.Beta()* bsm4v.Gamma()*300*abs(cos(bsm4v.Theta()) ))
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

def llpDelayProbability(lifetime, delayfactor):
    if delayfactor==0: return 0 #if not delay cannot pass
    #integral of a normalized exponential with tau*delayfactor, from cutpoint to cutpoint_high
    return exp(- time_cutpoint/(lifetime*delayfactor)) -exp(- time_cutpoint_high/(lifetime*delayfactor))

def llpDisplacementProbability(lifetime, jet, Rmin=0, Rmax=9e9, Zmin=0, Zmax=9e9, fiducial_range=None):
    if Rdisplacementfactor(jet)==0:
        if (Rmin or Zmin): return 0
        if not Rmin and not Zmin: return 1
        print "WTF",Rmin,Zmin
        sys.exit(1)
    if fiducial_range:
        if fiducial_ranges[fiducial_range][0] == 'randzup': 
            name, (Rmin, Rmax, Zmax) = fiducial_ranges[fiducial_range]
        elif fiducial_ranges[fiducial_range][0] == 'ronly': 
            name, (Rmin, Rmax) = fiducial_ranges[fiducial_range]
        elif fiducial_ranges[fiducial_range][0] == 'zonly': 
            name, (Zmin, Zmax) = fiducial_ranges[fiducial_range]
        else:
            print  "WTF",fiducial_ranges[fiducial_range]

    #integral of a normalized exponential with tau*delayfactor, from cutpoint to cutpoint_high
    Rdecay = lifetime*Rdisplacementfactor(jet)
    Zdecay = lifetime*Zdisplacementfactor(jet)
    if name == 'randzup':
        probr = exp(- Rmin/Rdecay)-exp(- Rmax/Rdecay)
        probz = 1.-exp(- Rmax/Rdecay)
    elif name == 'ronly':
        probr = exp(- Rmin/Rdecay)-exp(- Rmax/Rdecay)
        probz = 1.
    elif name == 'zonly':
        probr = 1.
        probz = exp(- Zmin/Zdecay)-exp(- Zmax/Zdecay)
    #print lifetime, Rdecay,Zdecay
    #print fiducial_ranges[fiducial_range]
    #print probr, probz
    return probr*probz

def passNominal(sortedjets, lifetime):
    nonepass  = 1 #will do 1-probability that none pass
    for jet in sortedjets:
        if recopt(jet) < 450: break
        nonepass *= 1- llpDisplacementProbability(lifetime, jet, fiducial_range='contained_for_delayed')
    return 1 - nonepass

def passDelayedHighPt(sortedjets, lifetime):
    return passDelayed(sortedjets, lifetime, 220)

def passCalRatio(llps, lifetime):
    nonepass  = 1 #will do 1-probability that none pass
    for llp in llps:
        pteff = calratioPtEfficiency(llp.pt())
        if abs(llp.eta())<1.4:
            poseff = calratioBarrelREfficiency(llp, lifetime)
        elif abs(llp.eta())<2.5:
            poseff = calratioEndcapZEfficiency(llp, lifetime)
        else:
            poseff = 0
        nonepass *= 1- pteff*poseff
    return 1 - nonepass

def passDelayedISR(sortedjets, lifetime):
    maxprob = 0
    for i1,jet1 in enumerate(sortedjets):
      if recopt(jet1) < 220: break
      maxprob = max(maxprob, passDelayed(sortedjets, lifetime) )
    return maxprob

def countL1J(sortedjets, l1ptcut):
    count = 0
    for jet in sortedjets:
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
def getHTweight(sortedjets):
    ht21, ht31 = 0,0
    for jet in sortedjets:
        if l1pt(jet) > 15 and abs( jet.eta() ) < 2.1: 
            ht21 += l1pt(jet)
        if l1pt(jet) > 20:
            ht31 += l1pt(jet)
    if ht21 >= ht21cut or ht31 >= ht31cut: return 1
    return 0

def passL1(sortedjets, lifetime):
    nonepass  = 1 #will do 1-probability that none pass
    weight_J100_pu = 3.7/28550
    weight_4J15_pu = 4.8/28550
    weight_HT21_pu   = 2.3/28550
    weight_HT31_pu   = 6.8/28550
    weight_3J15_pu = 27.1/28550
    weight_2J15_pu = 170./28550 #guess
    weight_J15_pu  = 1101./28550
    countJ15  = countL1J(sortedjets,15)
    countJ100 = countL1J(sortedjets,100)

    weight_J100, weight_4J15 = 0,0
    if countJ100: weight_J100 = 1
    if countJ15>=4:   weight_4J15 = 1
    elif countJ15==3: weight_4J15 = weight_J15_pu
    elif countJ15==2: weight_4J15 = weight_2J15_pu
    elif countJ15==1: weight_4J15 = weight_3J15_pu

    weight_HT = getHTweight(sortedjets)

    #print (weight_HT,weight_4J15,weight_J100,weight_HT_pu+weight_4J15_pu+weight_J100_pu)
    l1weight = max(weight_HT,weight_4J15,weight_J100,weight_HT31_pu+weight_4J15_pu+weight_J100_pu)
    #l1weight = max(weight_4J15,weight_J100,weight_HT31_pu+weight_4J15_pu+weight_J100_pu)
    return l1weight

def passTracking(llps, lifetime):
    nonepass  = 1 #will do 1-probability that none pass
    for llp in llps:
      nonepass *= (1-vtxEfficiency(lifetime, llp))
    hltprob = 1-nonepass
    return hltprob

def passDelayed2(sortedjets, lifetime, ptcut=40):
    nonepass  = 1 #will do 1-probability that none pass
    seenfactors = {}
    for i1,jet1 in enumerate(sortedjets):
      if recopt(jet1) < ptcut: break
      factor1 = delayfactor(jet1)
      for i2 in range(i1+1, len(sortedjets)):
        jet2 = sortedjets[i2]
        if recopt(jet2) < ptcut: break
        factor2 = delayfactor(jet2)
        prob1 =  llpDelayProbability(lifetime, delayfactor(jet1) )*llpDisplacementProbability(lifetime, jet1, fiducial_range='contained_for_delayed')
        prob2 =  llpDelayProbability(lifetime, delayfactor(jet2) )*llpDisplacementProbability(lifetime, jet2, fiducial_range='contained_for_delayed')
        prob = prob1*prob2
        factor = (factor1, factor2)
        if factor not in seenfactors or prob > seenfactors[factor]:
          seenfactors[factor] = prob
    for p in seenfactors.values():
      nonepass *= 1- p
    hltprob = 1-nonepass
    return hltprob


def passDelayed(sortedjets, lifetime, ptcut=40):
    nonepass  = 1 #will do 1-probability that none pass
    seenfactors = {}
    for jet in sortedjets:
      if recopt(jet) < ptcut: break
      factor = delayfactor(jet)
      prob =  llpDelayProbability(lifetime, delayfactor(jet) )*llpDisplacementProbability(lifetime, jet, fiducial_range='contained_for_delayed')
      if factor not in seenfactors or prob > seenfactors[factor]:
        seenfactors[factor] = prob
    for p in seenfactors.values():
      nonepass *= 1- p
    hltprob = 1-nonepass
    return hltprob

def decorateEvent(tree, llps):
    tree.pass_nominal = {}
    tree.pass_delayed_highpt = {}
    tree.pass_delayed_ISR = {}
    tree.pass_delayed_PU = {}
    tree.pass_delayed2_PU = {}
    tree.pass_tracking_PU = {}
    tree.pass_CalRatio = {}

    sortedjets = sorted(tree.AntiKt4TruthDressedWZJets, key=lambda x: recopt(x), reverse=True)

    htjetpt = sum([recopt(jet) for jet in tree.AntiKt4TruthDressedWZJets ])
    leadjetpt = max([recopt(jet) for jet in tree.AntiKt4TruthDressedWZJets ])
    leadN1jetpt = max([recopt(jet) for jet in tree.AntiKt4TruthDressedWZJets if isN1decay(jet) ])
    h_htjet.Fill(htjetpt)
    h_leadjet.Fill(leadjetpt)
    h_leadN1jet.Fill(leadN1jetpt)

    for lifetime in lifetimes:
        tree.pass_nominal[lifetime] = passNominal(sortedjets, lifetime)

        tree.pass_delayed_highpt[lifetime] = passDelayedHighPt(sortedjets, lifetime)
        tree.pass_delayed_ISR[lifetime] = passDelayedISR(sortedjets, lifetime)
        tracking = passTracking(llps, lifetime)
        delayed = passDelayed(sortedjets, lifetime)
        delayed2 = passDelayed2(sortedjets, lifetime)
        passl1  = passL1(sortedjets, lifetime)
        tree.pass_delayed2_PU[lifetime] = passl1*delayed2
        tree.pass_delayed_PU[lifetime] = passl1*delayed
        tree.pass_tracking_PU[lifetime] = passl1*tracking
        tree.pass_CalRatio[lifetime] = passCalRatio(llps,lifetime)
        #if leadN1jetpt > 220:
        #  if tree.pass_delayed_ISR[lifetime]> tree.pass_delayed_highpt[lifetime]:
        #    print  tree.pass_delayed_ISR[lifetime],tree.pass_delayed_highpt[lifetime], delayed
        #    for jet in sortedjets:
        #        print recopt(jet), delayfactor(jet), llpDelayProbability(lifetime, delayfactor(jet) )
        #    exit(1)

        h_passleadjet['nominal'][lifetime].Fill(leadjetpt, tree.pass_nominal[lifetime])
        h_passleadjet['delayed_highpt'][lifetime].Fill(leadjetpt, tree.pass_delayed_highpt[lifetime])
        h_passleadjet['delayed_ISR'][lifetime].Fill(leadjetpt, tree.pass_delayed_ISR[lifetime])
        h_passleadjet['CalRatio'][lifetime].Fill(leadjetpt, tree.pass_CalRatio[lifetime])
        h_passleadjet['delayed_PU'][lifetime].Fill(leadjetpt, tree.pass_delayed_PU[lifetime])
        h_passleadjet['delayed2_PU'][lifetime].Fill(leadjetpt, tree.pass_delayed2_PU[lifetime])
        h_passleadjet['tracking_PU'][lifetime].Fill(leadjetpt, tree.pass_tracking_PU[lifetime])
        h_passleadN1jet['nominal'][lifetime].Fill(leadN1jetpt, tree.pass_nominal[lifetime])
        h_passleadN1jet['delayed_highpt'][lifetime].Fill(leadN1jetpt, tree.pass_delayed_highpt[lifetime])
        h_passleadN1jet['delayed_ISR'][lifetime].Fill(leadN1jetpt, tree.pass_delayed_ISR[lifetime])
        h_passleadN1jet['CalRatio'][lifetime].Fill(leadN1jetpt, tree.pass_CalRatio[lifetime])
        h_passleadN1jet['delayed_PU'][lifetime].Fill(leadN1jetpt, tree.pass_delayed_PU[lifetime])
        h_passleadN1jet['delayed2_PU'][lifetime].Fill(leadN1jetpt, tree.pass_delayed2_PU[lifetime])
        h_passleadN1jet['tracking_PU'][lifetime].Fill(leadN1jetpt, tree.pass_tracking_PU[lifetime])
        h_passhtjet['nominal'][lifetime].Fill(htjetpt, tree.pass_nominal[lifetime])
        h_passhtjet['delayed_highpt'][lifetime].Fill(htjetpt, tree.pass_delayed_highpt[lifetime])
        h_passhtjet['delayed_ISR'][lifetime].Fill(htjetpt, tree.pass_delayed_ISR[lifetime])
        h_passhtjet['CalRatio'][lifetime].Fill(htjetpt, tree.pass_CalRatio[lifetime])
        h_passhtjet['delayed_PU'][lifetime].Fill(htjetpt, tree.pass_delayed_PU[lifetime])
        h_passhtjet['delayed2_PU'][lifetime].Fill(htjetpt, tree.pass_delayed2_PU[lifetime])
        h_passhtjet['tracking_PU'][lifetime].Fill(htjetpt, tree.pass_tracking_PU[lifetime])

        h_lifetime.Fill(lifetime)
        h_passlifetime['nominal'].Fill(lifetime, tree.pass_nominal[lifetime])
        h_passlifetime['delayed_highpt'].Fill(lifetime, tree.pass_delayed_highpt[lifetime])
        h_passlifetime['delayed_ISR'].Fill(lifetime, tree.pass_delayed_ISR[lifetime])
        h_passlifetime['CalRatio'].Fill(lifetime, tree.pass_CalRatio[lifetime])
        h_passlifetime['delayed_PU'].Fill(lifetime, tree.pass_delayed_PU[lifetime])
        h_passlifetime['delayed2_PU'].Fill(lifetime, tree.pass_delayed2_PU[lifetime])
        h_passlifetime['tracking_PU'].Fill(lifetime, tree.pass_tracking_PU[lifetime])

def main():
    for cut in cuts:
        h_passlifetime[cut] =  TH1F('h_lifetime'+cut,'h_lifetime'+cut,len(lifetimes),array.array('f',[0]+[l*2 for l in lifetimes]))
        h_passhtjet[cut] = {}
        h_passleadjet[cut] = {}
        h_passleadN1jet[cut] = {}
        for lifetime in lifetimes:
            h_passhtjet[cut][lifetime] = TH1F('h_htjet_%s_%f'%(cut,lifetime),'h_htjet',nbins,0,maxht)
            h_passleadjet[cut][lifetime] = TH1F('h_leadjet_%s_%f'%(cut,lifetime),'h_leadjet',nbins,0,maxjetpt)
            h_passleadN1jet[cut][lifetime] = TH1F('h_leadN1jet_%s_%f'%(cut,lifetime),'h_leadN1jet',nbins,0,maxjetpt)
    
    for i in range(4000):
    #for i in range(tree.GetEntries()):
        if i%1000==0: print i
        tree.GetEntry(i)
        dec1, dec2 = getDecays(tree)
        pdgs = [jet.auxdataConst['int']('TrueFlavor') for jet in tree.AntiKt4TruthDressedWZJets if jet.auxdataConst['int']('TrueFlavor')<4]
        for jet in tree.AntiKt4TruthDressedWZJets:
            ok = decorateJet(jet,dec1,dec2,pdgs)
            if not ok: break
            
        if not ok: continue #keep only events where partons are properly matched
        decorateEvent(tree, [dec1[0], dec2[0]]) 

        #if not any([delayfactor(jet) for jet in tree.AntiKt4TruthDressedWZJets]):
        #    print "WTF, not filled?"
        #    print dec1, dec2, pdgs
        #    for j in tree.AntiKt4TruthDressedWZJets:
        #        print j.pt(), recopt(j), j.auxdataConst['int']('TrueFlavor')
        #    exit(1)
    
    leg = TLegend(0.1,0.5,0.5,0.9)
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
        leg.AddEntry(effs[cut],cut,"LP")
        gPad.Update()
    line.Draw()
    leg.Draw()
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
        leg.Draw()
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
        leg.Draw()
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
        leg.Draw()
        canv.SaveAs("plots/%s/eff_leadN1_%f.png"%(sample,lifetime))

if __name__ == "__main__":  main()
