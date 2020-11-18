from ROOT import *
from math import erfc
gROOT.SetBatch(1)

ran = TRandom3(123)

reco2hlt_res = TF1("reco2hlt_res","[0]/x+[1]",0,1000)
reco2hlt_res.SetParameters(2, 0.1)
reco2l1_res = TF1("reco2l1_res","[0]/x+[1]",0,1000)
minjetpt = 20
maxjetpt = 250

def calibrate(x, a, b):
    return (x-a)/(b)

def fitFcn(x, p):
    if x<10: return 0
    calibrated = max(1, calibrate(x[0],p[0],p[1]))
    reco2l1_res.SetParameters(p[2],p[3])
    res = calibrated*reco2l1_res.Eval(x[0])
    if x[0]<70: cutvalue = 15
    elif x[0]<130: cutvalue = 50
    else: cutvalue = 100
    prob = erfc((cutvalue - calibrated)/res)/2.
    return prob

fcn = TF1("fcn",fitFcn,minjetpt, maxjetpt ,4)
fcn.SetParameters(4,1.5,7.06285e+00, 2.74922e-02)
fcn.FixParameter(2)
fcn.FixParameter(3)

hdata = TH1F("hdata","hdata",50,0,maxjetpt)
hdata.SetBinContent( 1, 0 ) #0
hdata.SetBinContent( 2, 0 ) 
hdata.SetBinContent( 3, 0 ) 
hdata.SetBinContent( 4, 0 ) 
hdata.SetBinContent( 5, 0 )
hdata.SetBinContent( 6, 0.01 )
hdata.SetBinContent( 7, 0.06 )
hdata.SetBinContent( 8, 0.25 )
hdata.SetBinContent( 9, 0.52 )
hdata.SetBinContent(10, 0.75 )
hdata.SetBinContent(11, 0.88 ) #50
hdata.SetBinContent(12, 0.94 )
hdata.SetBinContent(13, 0.98 )
hdata.SetBinContent(14, 0.99 )
hdata.SetBinContent(15, 0.01 ) #70 - J50
hdata.SetBinContent(16, 0.07 )
hdata.SetBinContent(17, 0.18 )
hdata.SetBinContent(18, 0.35 )
hdata.SetBinContent(19, 0.55 )
hdata.SetBinContent(20, 0.71 )
hdata.SetBinContent(21, 0.83 ) #100
hdata.SetBinContent(22, 0.90 )
hdata.SetBinContent(23, 0.95 )
hdata.SetBinContent(24, 0.96 )
hdata.SetBinContent(25, 0.98 ) #120
hdata.SetBinContent(26, 0.99 )
hdata.SetBinContent(27, 0.02 )
hdata.SetBinContent(28, 0.05 )
hdata.SetBinContent(29, 0.12 )
hdata.SetBinContent(30, 0.21 )
hdata.SetBinContent(31, 0.33 ) #150
hdata.SetBinContent(32, 0.48 )
hdata.SetBinContent(33, 0.60 )
hdata.SetBinContent(34, 0.71 )
hdata.SetBinContent(35, 0.80 )
hdata.SetBinContent(36, 0.87 )
hdata.SetBinContent(37, 0.92 )
hdata.SetBinContent(38, 0.94 )
hdata.SetBinContent(39, 0.96 )
hdata.SetBinContent(40, 0.97 )
hdata.SetBinContent(41, 0.98 ) #200
hdata.SetBinContent(42, 0.99 )
hdata.SetBinContent(43, 0.99 )
hdata.SetBinContent(44, 1.0 )
hdata.SetBinContent(45, 1.0 )
hdata.SetBinContent(46, 1.0 )
hdata.SetBinContent(47, 1.0 )
hdata.SetBinContent(48, 1.0 )
hdata.SetBinContent(49, 1.0 )
hdata.SetBinContent(50, 1.0 )
can = TCanvas()
hdata.Fit(fcn,'W,R')
hdata.Draw()
can.SaveAs("plots/L1Jfit.png")
