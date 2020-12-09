from math import exp, sin, cos, log, sqrt
import array
import os, sys
from collections import OrderedDict
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("cutset")
parser.add_argument("-n", "--nevents", default=9e9, type=int)
parser.add_argument("-s", "--sample", default="higgsino_200")
parser.add_argument("-d", "--delayed", action="store_true", help="Select only events with j40+delayed")
opts = parser.parse_args()

from ROOT import *
gROOT.SetBatch(1)
gROOT.Macro('$ROOTCOREDIR/scripts/load_packages.C')
# Initialize the xAOD infrastructure
xAOD.Init()
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

samples = {
    'higgsino_300':"/afs/cern.ch/work/j/jmontejo/LLPtriggertiming/LLPtrigger_samples/n3n4_1TeVsquark_300gev.DAOD_TRUTH3.root",
    'higgsino_200':"/afs/cern.ch/work/j/jmontejo/LLPtriggertiming/LLPtrigger_samples/n3n4_1TeVsquark_200gev.DAOD_TRUTH3.root",
    'higgsino_250':"/afs/cern.ch/work/j/jmontejo/LLPtriggertiming/LLPtrigger_samples/n3n4_1TeVsquark_250gev.DAOD_TRUTH3.root",
    'higgsino_150':"/afs/cern.ch/work/j/jmontejo/LLPtriggertiming/LLPtrigger_samples/n3n4_1TeVsquark_150gev.DAOD_TRUTH3.root",
    'higgsino_150_pythiadecay':"/afs/cern.ch/work/j/jmontejo/LLPtriggertiming/LLPtrigger_samples/n3n4_1TeVsquark_150gev_pythiadecay.DAOD_TRUTH3.root",
    'higgs_ss55': '/afs/cern.ch/work/j/jmontejo/LLPtriggertiming/LLPtrigger_samples/Higgs_to_scalars/mc16_13TeV.311314.MadGraphPythia8EvtGen_A14NNPDF31LO_HSS_LLP_mH125_mS55_ltlow.deriv.DAOD_TRUTH3.e7270_p3401.root', 
    'higgs_ss35': '/afs/cern.ch/work/j/jmontejo/LLPtriggertiming/LLPtrigger_samples/Higgs_to_scalars/mc15_13TeV.311312.MadGraphPythia8EvtGen_A14NNPDF31LO_HSS_LLP_mH125_mS35_ltlow.deriv.DAOD_TRUTH3.e7270_p3401.root',
    'my_hss35': '/eos/home-j/jmontejo/LLPtrigger_samples/Higgs_to_scalars/run/DAOD_TRUTH3.mc16_13TeV.311312.MadGraphPythia8EvtGen_A14NNPDF31LO_HSS_LLP_mH125_mS35_ltlow.root',
}

rfile = TFile.Open(samples[opts.sample])
if opts.delayed: opts.sample += "_selectdelayed"
tree  = xAOD.MakeTransientTree(rfile, "CollectionTree")
recopt = xAOD.Jet.Decorator('float')('recopt')
l1pt = xAOD.Jet.Decorator('float')('l1pt')
decaytime = xAOD.Jet.Decorator('float')('decaytime')
n1index = xAOD.Jet.Decorator('int')('n1index')
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
maxjet4pt = 200
#time_cutpoint = 2
time_cutpoint_high = 12
#>>> ROOT.Math.normal_cdf_c(1)
#0.15865525393145707
#>>> ROOT.Math.normal_cdf_c(2)
#0.022750131948179216 (*6kHz input = 136 Hz)
#>>> ROOT.Math.normal_cdf_c(1)**2
#0.025171489600055125
#>>> ROOT.Math.normal_cdf_c(2)**2
#0.0005175685036595646 (*6kHz input = 3 Hz)
#>>> ROOT.Math.normal_cdf_c(1.5)**2
#0.004463202141377714 (*6kHz input = 27 Hz)
fiducial_ranges = {
    'contained_for_delayedL1' : ('randzup',(0,2000,3000)),
    'contained_for_delayed' : ('randzup',(0,1500,3000)),
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
lifetimes = sorted([pow(10,i) for i in range(-2,3)]+[3*pow(10,i) for i in range(-3,3)])

h_lifetime = TH1F('h_lifetime','h_lifetime',len(lifetimes),array.array('f',[0]+[l*2 for l in lifetimes]))
h_htjet = TH1F('h_htjet','h_htjet',nbins,0,maxht)
h_leadjet = TH1F('h_leadjet','h_leadjet',nbins,0,maxjetpt)
h_jet4 = TH1F('h_jet4','h_jet4',nbins,0,maxjet4pt)
h_leadN1jet = TH1F('h_leadN1jet','h_leadN1jet',nbins,0,maxjetpt)
h_passlifetime = {}
h_passhtjet = {}
h_passleadjet = {}
h_passjet4 = {}
h_passleadN1jet = {}
h_passdelayedhtjet = {}
h_passdelayedleadjet = {}
h_passdelayedjet4 = {}
h_passdelayedleadN1jet = {}

def calratioEndcapZEfficiency(llp, lifetime):
    eff = 0
    eff += 0.5*llpDisplacementProbability(lifetime, llp, fiducial_range='calratio_endcap_range1')
    eff += 1.0*llpDisplacementProbability(lifetime, llp, fiducial_range='calratio_endcap_range2')
    eff += 0.5*llpDisplacementProbability(lifetime, llp, fiducial_range='calratio_endcap_range3')
    assert eff<=1, eff
    return eff

def calratioBarrelREfficiency(llp, lifetime):
    eff = 0
    eff += 0.5*llpDisplacementProbability(lifetime, llp, fiducial_range='calratio_barrel_range1')
    eff += 1.0*llpDisplacementProbability(lifetime, llp, fiducial_range='calratio_barrel_range2')
    eff += 0.5*llpDisplacementProbability(lifetime, llp, fiducial_range='calratio_barrel_range3')
    assert eff<=1, eff
    return eff

def calratioPtEfficiency(LLPpt):
    minpt = 20
    factor = 1
    if LLPpt<minpt: return 0
    if LLPpt>260: return 0.8/factor
    return log((LLPpt-minpt)/10.)/(4.*factor)

def vtxEfficiency(lifetime, llp):
    ff = 0.5
    eff = 0
    eff += 0.9*ff*llpDisplacementProbability(lifetime, llp, fiducial_range='vtx_range1')
    eff += 0.8*ff*llpDisplacementProbability(lifetime, llp, fiducial_range='vtx_range2')
    eff += 0.7*ff*llpDisplacementProbability(lifetime, llp, fiducial_range='vtx_range3')
    eff += 0.5*ff*llpDisplacementProbability(lifetime, llp, fiducial_range='vtx_range4')
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
      if 'higgsino' in opts.sample and bsm.status()==62:
        toret.append([bsm])
      elif 'hss' in opts.sample and bsm.status()==22:
        toret.append([bsm])
      else: continue
      toret[-1].extend([bsm.child(i) for i in range(bsm.nChildren())])
      bsm4v = bsm.p4()
      decaytime.set(bsm, rnd.Exp(1) )
      Rdisplacementfactor.set(bsm, decaytime(bsm)*bsm4v.Beta()* bsm4v.Gamma()*300*sin(bsm4v.Theta()) )
      Zdisplacementfactor.set(bsm, decaytime(bsm)*bsm4v.Beta()* bsm4v.Gamma()*300*abs(cos(bsm4v.Theta()) ))
      if len(toret)==2: break
    return toret

def decorateJet(jet, dec1, dec2, pdgs=None):
    recopt.set(jet, truth2reco(jet.pt()*gev) )
    l1pt.set(jet, reco2l1(recopt(jet)) )
    flav = jet.auxdataConst['int']('TrueFlavor') 
    dr1, dr2 = 9, 9
    for p in dec1:
        if p.pdgId() == flav: 
            dr1 = p.p4().DeltaR(jet.p4())
            break
    for p in dec2:
        if p.pdgId() == flav: 
            dr2 = p.p4().DeltaR(jet.p4())
            break

    n1 = None
    if dr1 < dr2 and dr1<0.6:
        n1 = dec1[0]
        n1index.set(jet, 0 )
    elif dr2 < dr1 and dr2<0.6:
        n1 = dec2[0]
        n1index.set(jet, 1 )
    elif min(dr1,dr2) < 9 and min(dr1,dr2>0.6):
        if 'higgsino' in opts.sample and pdgs.count(flav) <=2:
            #print "No good dR matching",dr1,dr2,flav,jet.pt(),jet.eta(),jet.phi(),[(p.pdgId(), p.pt(), p.eta(),p.phi()) for p in dec1+dec2 if p.pdgId()==flav]
            #print pdgs
            return False
        if 'hss' in opts.sample and pdgs.count(flav) <=4:
            return False
    if n1:
        #n14v = n1.p4()
        #delayfactor.set(jet, n14v.Gamma()*(1-n14v.Beta()*cos( n1.phi() -jet.phi() )))
        #deltal = getDeltaL(jet, n1)
        #time difference =  gamma*tau(1 - beta*cos) + sqrt( (beta*gamma*tau*sin)^2 + (d/c)^2 )  - d/c
        Rdisplacementfactor.set(jet, Rdisplacementfactor(n1) )
        Zdisplacementfactor.set(jet, Zdisplacementfactor(n1) )
        #print "N1 pt, beta*gamma: ",n1.pt(), n14v.Beta()* n14v.Gamma()
        isN1decay.set(jet, 1 )
    else:
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

def llpDelayProbability(delay, time_cutpoint):
    if delay==0: return 0 #if not delay cannot pass
    if delay < time_cutpoint or delay > time_cutpoint_high: return 0
    return 1

def llpDisplacementProbability(lifetime, jetorllp, Rmin=0, Rmax=9e9, Zmin=0, Zmax=9e9, fiducial_range=None):
    if Rdisplacementfactor(jetorllp)==0:
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

    Rdecay = lifetime*Rdisplacementfactor(jetorllp)
    Zdecay = lifetime*Zdisplacementfactor(jetorllp)
    if name == 'randzup':
        if Rdecay < Rmin or Rdecay > Rmax: return 0
        if Zdecay > Zmax: return 0
        return 1
    elif name == 'ronly':
        if Rdecay < Rmin or Rdecay > Rmax: return 0
        return 1
    elif name == 'zonly':
        if Zdecay < Zmin or Zdecay > Zmax: return 0
        return 1
    raise WTFError

def passPlateauJ100(sortedjets, lifetime, unused):
    return len(sortedjets)>=1 and recopt(sortedjets[0])>200
def passPlateau4J15(sortedjets, lifetime, unused):
    return len(sortedjets)>=4 and recopt(sortedjets[3])>60
def passPlateauHT(sortedjets, lifetime, unused):
    return len(sortedjets)>=1 and sum([recopt(j) for j in sortedjets])>550
def passStreamer(sortedjets, lifetime, unused):
    return 1

def passNominal(sortedjets, lifetime, unused):
    for jet in sortedjets:
        if recopt(jet) < 450: break
        if llpDisplacementProbability(lifetime, jet, fiducial_range='contained_for_delayed'):
            return 1
    count4j120 = 0
    for jet in sortedjets:
        if recopt(jet) < 120: break
        if llpDisplacementProbability(lifetime, jet, fiducial_range='contained_for_delayed'):
            count4j120 += 1
    if count4j120 >=4: return 1
    sumht = 0
    for jet in sortedjets:
        if llpDisplacementProbability(lifetime, jet, fiducial_range='contained_for_delayed'):
            sumht += recopt(jet)
    if sumht >=1000: return 1
    return 0

def passDelayedHighPt(sortedjets, lifetime, llps):
    return passDelayed(sortedjets, lifetime, llps, ptcut=200)

def passCalRatio(sortedjets, lifetime, llps):
    if len(sortedjets)==0 or recopt(sortedjets[0]) < 60: return 0
    nonepass  = 1 #will do 1-probability that none pass
    for llp in llps:
        pteff = calratioPtEfficiency(llp.pt()*gev)
        if abs(llp.eta())<1.4:
            poseff = calratioBarrelREfficiency(llp, lifetime)
        elif abs(llp.eta())<2.5:
            poseff = calratioEndcapZEfficiency(llp, lifetime)
        else:
            poseff = 0
        nonepass *= 1- pteff*poseff
    return 1 - nonepass

def passDelayedISR(sortedjets, lifetime, llps):
    if len(sortedjets)==0 or recopt(sortedjets[0]) < 200: return 0
    return passDelayed(sortedjets, lifetime, llps, count=1, timecut=1)

def countL1J(sortedjets, l1ptcut, lifetime):
    count = 0
    for jet in sortedjets:
        if not llpDisplacementProbability(lifetime, jet, fiducial_range='contained_for_delayedL1'): continue
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
def getHTweight(sortedjets, lifetime, only21=False, only31=False):
    ht21, ht31 = 0,0
    for jet in sortedjets:
        if not llpDisplacementProbability(lifetime, jet, fiducial_range='contained_for_delayed'): continue
        if l1pt(jet) > 15 and abs( jet.eta() ) < 2.1: 
            ht21 += l1pt(jet)
        if l1pt(jet) > 20:
            ht31 += l1pt(jet)
    if only21:
        return ht21 >= ht21cut
    if only31:
        return ht31 >= ht31cut
    if ht21 >= ht21cut or ht31 >= ht31cut: return 1
    return 0

def passJ100(sortedjets, lifetime):
    return countL1J(sortedjets,100,lifetime)>=1
def pass4J15(sortedjets, lifetime):
    countJ15  = countL1J(sortedjets,15,lifetime)
    return countJ15>=4
def passHT(sortedjets, lifetime,only21=False, only31=False):
    return getHTweight(sortedjets, lifetime,only21,only31)

def passL1(sortedjets, lifetime):
    weight_J100_pu = 3.7/28550
    weight_4J15_pu = 4.8/28550
    weight_HT21_pu   = 2.3/28550
    weight_HT31_pu   = 6.8/28550
    weight_3J15_pu = 27.1/28550
    weight_2J15_pu = 170./28550 #guess
    weight_J15_pu  = 1101./28550
    countJ15  = countL1J(sortedjets,15,lifetime)
    countJ100 = countL1J(sortedjets,100,lifetime)

    weight_J100, weight_4J15 = 0,0
    if countJ100: weight_J100 = 1
    if countJ15>=4:   weight_4J15 = 1
    elif countJ15==3: weight_4J15 = weight_J15_pu
    elif countJ15==2: weight_4J15 = weight_2J15_pu
    elif countJ15==1: weight_4J15 = weight_3J15_pu
    weight_HT = getHTweight(sortedjets, lifetime)
    l1weight = max(weight_HT,weight_4J15,weight_J100,weight_HT31_pu+weight_4J15_pu+weight_J100_pu)
    return l1weight

def passTracking4j60(sortedjets, lifetime, llps):
    if len(sortedjets)>=4 and recopt(sortedjets[3]) > 60:
        return passTracking(sortedjets, lifetime, llps)
    return 0
def passTrackingj300_4j60(sortedjets, lifetime, llps):
    if len(sortedjets)>=1 and recopt(sortedjets[0]) > 300:
        return passTracking(sortedjets, lifetime, llps)
    if len(sortedjets)>=4 and recopt(sortedjets[3]) > 60:
        return passTracking(sortedjets, lifetime, llps)
    return 0
def passTrackingAndDelayed(sortedjets, lifetime, llps):
    return passTracking(sortedjets, lifetime, llps)*passDelayed(sortedjets, lifetime, llps, count=1, timecut=1)

def passTracking(unused, lifetime, llps):
    nonepass  = 1 #will do 1-probability that none pass
    for llp in llps:
      nonepass *= (1-vtxEfficiency(lifetime, llp))
    hltprob = 1-nonepass
    return hltprob

def getDelay(jet, llps, lifetime):
    n1 = llps[ n1index(jet) ]
    n14v = n1.p4()
    length = 1500
    c = 300
    decayt = lifetime*decaytime(n1)
    timediff = n14v.Gamma()*decayt*(1-n14v.Beta()*cos( n1.phi() -jet.phi() )) + \
               sqrt( pow(n14v.Beta()*n14v.Gamma()*decayt*sin( n1.phi() -jet.phi() ), 2) + pow(length/c,2) ) - length/c
    return timediff
    
    #time difference =  gamma*tau(1 - beta*cos) + sqrt( (beta*gamma*tau*sin)^2 + (d/c)^2 )  - d/c


def passDelayed1_4J15plateau(sortedjets, lifetime, llps):
    return passDelayed(sortedjets, lifetime, llps, count=1, timecut=1)*passPlateau4J15(sortedjets, lifetime, llps)
def passDelayed1or2(sortedjets, lifetime, llps):
    return passDelayed1(sortedjets, lifetime, llps) or passDelayed2(sortedjets, lifetime, llps)
def passDelayed1(sortedjets, lifetime, llps):
    return passDelayed(sortedjets, lifetime, llps, count=1, timecut=2.65)
def passDelayed2(sortedjets, lifetime, llps):
    return passDelayed(sortedjets, lifetime, llps, count=2, timecut=1.5)
def passDelayed3(sortedjets, lifetime, llps):
    return passDelayed(sortedjets, lifetime, llps, count=3, timecut=0.99)

def passDelayed(sortedjets, lifetime, llps, ptcut=40, count=1, timecut=2):
    dopass = 0
    for jet in sortedjets:
      if recopt(jet) < ptcut: break
      dopass +=  llpDelayProbability(getDelay(jet, llps, lifetime), timecut )*llpDisplacementProbability(lifetime, jet, fiducial_range='contained_for_delayed')
      if dopass == count: return 1
    return 0

def decorateEvent(tree, llps, cuts):
    pass_cuts = {}

    sortedjets = sorted(tree.AntiKt4TruthDressedWZJets, key=lambda x: recopt(x), reverse=True)

    htjetpt = sum([recopt(jet) for jet in sortedjets])
    leadjetpt, leadN1jetpt = 0, 0
    if  len(sortedjets):
        leadjetpt = max([recopt(jet) for jet in sortedjets])
        n1jets = [recopt(jet) for jet in sortedjets if isN1decay(jet) ]
        if len(n1jets):
            leadN1jetpt = max(n1jets)
    jet4pt = 0
    if  len(sortedjets)>=4:
        jet4pt = recopt(sortedjets[3])
    h_htjet.Fill(htjetpt)
    h_leadjet.Fill(leadjetpt)
    h_leadN1jet.Fill(leadN1jetpt)
    h_jet4.Fill(jet4pt)

    for lifetime in lifetimes:
        h_lifetime.Fill(lifetime)
        pass_cuts[lifetime] = {}
        pL1   = passL1(sortedjets, lifetime)
        pJ100 = passJ100(sortedjets, lifetime)
        p4J15 = pass4J15(sortedjets, lifetime)
        pHT21 = passHT(sortedjets, lifetime, only21=True)
        pHT31 = passHT(sortedjets, lifetime, only31=True)
        for cutname, cutfunc in cuts.iteritems():
            if opts.delayed and not passDelayed1(sortedjets, lifetime, llps): continue

            pass_cuts[lifetime][cutname] = cutfunc(sortedjets, lifetime, llps)
            if "pass L1" in cutname:
                pass_cuts[lifetime][cutname] *= pL1
            if "pass J100" in cutname:
                pass_cuts[lifetime][cutname] *= pJ100
            if "pass 4J15" in cutname:
                pass_cuts[lifetime][cutname] *= p4J15
            if "pass HT21" in cutname:
                pass_cuts[lifetime][cutname] *= pHT21
            elif "pass HT31" in cutname:
                pass_cuts[lifetime][cutname] *= pHT31
            elif "pass HT" in cutname:
                pass_cuts[lifetime][cutname] *= max(pHT21,pHT31)

            h_passleadjet[cutname][lifetime].Fill(leadjetpt, pass_cuts[lifetime][cutname])
            h_passjet4[cutname][lifetime].Fill(jet4pt, pass_cuts[lifetime][cutname])
            h_passleadN1jet[cutname][lifetime].Fill(leadN1jetpt, pass_cuts[lifetime][cutname])
            h_passhtjet[cutname][lifetime].Fill(htjetpt, pass_cuts[lifetime][cutname])
            h_passlifetime[cutname].Fill(lifetime, pass_cuts[lifetime][cutname])

            h_passdelayedleadjet[cutname][lifetime].Fill(leadjetpt)
            h_passdelayedjet4[cutname][lifetime].Fill(jet4pt)
            h_passdelayedleadN1jet[cutname][lifetime].Fill(leadN1jetpt)
            h_passdelayedhtjet[cutname][lifetime].Fill(htjetpt)

def getColor(i):
    i += 1
    if i==3: return 4
    if i==4: return 8
    if i==5: return 28
    return i

def main():
    cutsets = {
    'l1' : OrderedDict([
        ('pass L1', passStreamer),
        ('pass 4J15', passStreamer),
        ('pass HT',  passStreamer),
        ('pass J100', passStreamer),
        #('streamer_pass HT21',  passStreamer),
        #('streamer_pass HT31',  passStreamer),
        #('plateau_pass J100', passPlateauJ100),
        #('plateau_pass 4J15', passPlateau4J15),
        #('plateau_pass HT',   passPlateauHT),
    ]),
    'hlt' : OrderedDict([
        ('Standard triggers',         passNominal),
        #('delayed_highpt',  passDelayedHighPt),
        #('j300 and delayed',     passDelayedISR),
        ('j300 or 4j60, and displaced vertex',   passTrackingj300_4j60),
        #('4j60 and displaced vertex',   passTracking4j60),
        #('delayed1_pass4J15', passDelayed1),
        #('delayed1_passHT', passDelayed1),
        #('pass L1 and 1 or 2 delayed', passDelayed1or2),
        ('pass L1 and 1 delayed j40 (2 sigma)', passDelayed1),
        ('pass L1 and 2 delayed j40 (1 sigma)', passDelayed2),
        #('4j60 and 1 delayed (1 sigma)', passDelayed1_4J15plateau),
        #('tracking_pass L1', passTracking),
        #('tracking_delayed_pass L1',   passTrackingAndDelayed),
        ('CalRatio trigger',        passCalRatio),
    ]),
    }
    cuts = cutsets[opts.cutset]
    for cut in cuts.iterkeys():
        h_passlifetime[cut] =  TH1F('h_lifetime'+cut,'h_lifetime'+cut,len(lifetimes),array.array('f',[0]+[l*2 for l in lifetimes]))
        h_passhtjet[cut] = {}
        h_passleadjet[cut] = {}
        h_passjet4[cut] = {}
        h_passleadN1jet[cut] = {}
        h_passdelayedhtjet[cut] = {}
        h_passdelayedleadjet[cut] = {}
        h_passdelayedjet4[cut] = {}
        h_passdelayedleadN1jet[cut] = {}
        for lifetime in lifetimes:
            h_passhtjet[cut][lifetime] = TH1F('h_htjet_%s_%f'%(cut,lifetime),'h_htjet',nbins,0,maxht)
            h_passleadjet[cut][lifetime] = TH1F('h_leadjet_%s_%f'%(cut,lifetime),'h_leadjet',nbins,0,maxjetpt)
            h_passjet4[cut][lifetime] = TH1F('h_jet4_%s_%f'%(cut,lifetime),'h_jet4',nbins,0,maxjet4pt)
            h_passleadN1jet[cut][lifetime] = TH1F('h_leadN1jet_%s_%f'%(cut,lifetime),'h_leadN1jet',nbins,0,maxjetpt)
            h_passdelayedhtjet[cut][lifetime] = TH1F('hh_htjet_%s_%f'%(cut,lifetime),'h_htjet',nbins,0,maxht)
            h_passdelayedleadjet[cut][lifetime] = TH1F('hh_leadjet_%s_%f'%(cut,lifetime),'h_leadjet',nbins,0,maxjetpt)
            h_passdelayedjet4[cut][lifetime] = TH1F('hh_jet4_%s_%f'%(cut,lifetime),'h_jet4',nbins,0,maxjet4pt)
            h_passdelayedleadN1jet[cut][lifetime] = TH1F('hh_leadN1jet_%s_%f'%(cut,lifetime),'h_leadN1jet',nbins,0,maxjetpt)
    
    for i in range(min(opts.nevents, tree.GetEntries())):
        if i%1000==0: print i
        tree.GetEntry(i)
        dec1, dec2 = getDecays(tree)
        if 'higgsino' in opts.sample:
            pdgs = [jet.auxdataConst['int']('TrueFlavor') for jet in tree.AntiKt4TruthDressedWZJets if jet.auxdataConst['int']('TrueFlavor')<4]
        elif 'hss' in opts.sample:
            pdgs = [jet.auxdataConst['int']('TrueFlavor') for jet in tree.AntiKt4TruthDressedWZJets if jet.auxdataConst['int']('TrueFlavor')==5]
        else:
            raise ImplementMe

        for jet in tree.AntiKt4TruthDressedWZJets:
            ok = decorateJet(jet,dec1,dec2,pdgs)
            if not ok: break
            
        if not ok: continue #keep only events where partons are properly matched
        decorateEvent(tree, [dec1[0], dec2[0]], cuts) 

        #if not any([delayfactor(jet) for jet in tree.AntiKt4TruthDressedWZJets]):
        #    print "WTF, not filled?"
        #    print dec1, dec2, pdgs
        #    for j in tree.AntiKt4TruthDressedWZJets:
        #        print j.pt(), recopt(j), j.auxdataConst['int']('TrueFlavor')
        #    exit(1)
    
    leg = TLegend(0.1,0.6,0.5,0.9)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    os.system('mkdir -p plots_sample/%s/%s'%(opts.sample,opts.cutset))
    canv = TCanvas()
    line = TLine(0,1,maxjetpt,1)
    effs = {}
    for i,cut in enumerate(cuts.keys()):
        print cut
        effs[cut] = TEfficiency(h_passlifetime[cut], h_lifetime)
        effs[cut].SetTitle("Efficiency ht;LLP proper lifetime [ns];Acceptance")
        effs[cut].SetMarkerStyle(20)
        effs[cut].SetLineColor(getColor(i))
        effs[cut].SetMarkerColor(getColor(i))
        effs[cut].Draw("ALP" if i==0 else "same LP")
        leg.AddEntry(effs[cut],cut,"LP")
        gPad.Update()
    leg.Draw()
    canv.SetLogx()
    gPad.Update()
    effs[cuts.keys()[0]].GetPaintedGraph().GetYaxis().SetRangeUser(0,1)
    canv.SaveAs("plots_sample/%s/%s/eff_lifetime.png"%(opts.sample,opts.cutset))
    canv.SetLogx(0)
    
    norm = nbins/4.
    for lifetime in lifetimes:
        effs_ht = {}
        effs_lead = {}
        effs_jet4 = {}
        effs_leadN1 = {}
    
        h = h_htjet.DrawNormalized("",norm)
        h.GetYaxis().SetRangeUser(0,1.0)
        h.SetTitle("Efficiency ht;H_{T} [GeV];Acceptance")
        for i,cut in enumerate(cuts):
            effs_ht[cut] = TEfficiency(h_passhtjet[cut][lifetime], h_htjet if not opts.delayed else h_passdelayedhtjet[cut][lifetime])
            effs_ht[cut].SetMarkerStyle(20)
            effs_ht[cut].SetLineColor(getColor(i))
            effs_ht[cut].SetMarkerColor(getColor(i))
            effs_ht[cut].Draw("same LP")
        line.Draw()
        leg.Draw()
        canv.SaveAs("plots_sample/%s/%s/eff_ht_%f.png"%(opts.sample,opts.cutset,lifetime))
    
        h = h_leadjet.DrawNormalized("",norm)
        h.GetYaxis().SetRangeUser(0,1.0)
        h.SetTitle("Efficiency lead;Leading jet p_{T} [GeV];Acceptance")
        for i,cut in enumerate(cuts):
            effs_lead[cut] = TEfficiency(h_passleadjet[cut][lifetime], h_leadjet if not opts.delayed else h_passdelayedleadjet[cut][lifetime])
            effs_lead[cut].SetMarkerStyle(20)
            effs_lead[cut].SetLineColor(getColor(i))
            effs_lead[cut].SetMarkerColor(getColor(i))
            effs_lead[cut].Draw("same LP")
        line.Draw()
        leg.Draw()
        canv.SaveAs("plots_sample/%s/%s/eff_lead_%f.png"%(opts.sample,opts.cutset,lifetime))
    
        h = h_jet4.DrawNormalized("",norm)
        h.GetYaxis().SetRangeUser(0,1.0)
        h.SetTitle("Efficiency jet4;Fourth leading jet p_{T} [GeV];Acceptance")
        for i,cut in enumerate(cuts):
            effs_jet4[cut] = TEfficiency(h_passjet4[cut][lifetime], h_jet4 if not opts.delayed else h_passdelayedjet4[cut][lifetime])
            effs_jet4[cut].SetMarkerStyle(20)
            effs_jet4[cut].SetLineColor(getColor(i))
            effs_jet4[cut].SetMarkerColor(getColor(i))
            effs_jet4[cut].Draw("same LP")
        line.Draw()
        leg.Draw()
        canv.SaveAs("plots_sample/%s/%s/eff_jet4_%f.png"%(opts.sample,opts.cutset,lifetime))
    
        h = h_leadN1jet.DrawNormalized("",norm)
        h.GetYaxis().SetRangeUser(0,1.0)
        h.SetTitle("Efficiency leadN1;Leading jet p_{T} from LLP decay [GeV];Acceptance")
        for i,cut in enumerate(cuts):
            effs_leadN1[cut] = TEfficiency(h_passleadN1jet[cut][lifetime], h_leadN1jet if not opts.delayed else h_passdelayedleadN1jet[cut][lifetime])
            effs_leadN1[cut].SetMarkerStyle(20)
            effs_leadN1[cut].SetLineColor(getColor(i))
            effs_leadN1[cut].SetMarkerColor(getColor(i))
            effs_leadN1[cut].Draw("same LP")
        line.Draw()
        leg.Draw()
        canv.SaveAs("plots_sample/%s/%s/eff_leadN1_%f.png"%(opts.sample,opts.cutset,lifetime))

if __name__ == "__main__":  main()
