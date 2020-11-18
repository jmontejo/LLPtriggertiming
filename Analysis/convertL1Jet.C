#include <vector>
#include <algorithm>
#include <iostream>

float toL1(float jet_pt, float jet_rnd){
    if (jet_pt < 25) return 0;
    float calib = (jet_pt - 2.15239e+01)/1.38092e+00;
    float smear = calib*(1+(0.042+17.93/jet_pt)*jet_rnd);
    return smear;
}

float convertL1Jet( int pos=0,
    float jet_pt0=0, float jet_eta0=0, float jet_rnd0=0,
    float jet_pt1=0, float jet_eta1=0, float jet_rnd1=0,
    float jet_pt2=0, float jet_eta2=0, float jet_rnd2=0,
    float jet_pt3=0, float jet_eta3=0, float jet_rnd3=0,
    float jet_pt4=0, float jet_eta4=0, float jet_rnd4=0,
    float jet_pt5=0, float jet_eta5=0, float jet_rnd5=0,
    float jet_pt6=0, float jet_eta6=0, float jet_rnd6=0,
    float jet_pt7=0, float jet_eta7=0, float jet_rnd7=0,
    float eta_threshold = 3.1
    ){

    vector<float> jets;
    if(jet_pt0>0 and fabs(jet_eta0)<eta_threshold) jets.push_back( toL1(jet_pt0,jet_rnd0) );
    if(jet_pt1>0 and fabs(jet_eta1)<eta_threshold) jets.push_back( toL1(jet_pt1,jet_rnd1) );
    if(jet_pt2>0 and fabs(jet_eta2)<eta_threshold) jets.push_back( toL1(jet_pt2,jet_rnd2) );
    if(jet_pt3>0 and fabs(jet_eta3)<eta_threshold) jets.push_back( toL1(jet_pt3,jet_rnd3) );
    if(jet_pt4>0 and fabs(jet_eta4)<eta_threshold) jets.push_back( toL1(jet_pt4,jet_rnd4) );
    if(jet_pt5>0 and fabs(jet_eta5)<eta_threshold) jets.push_back( toL1(jet_pt5,jet_rnd5) );
    if(jet_pt6>0 and fabs(jet_eta6)<eta_threshold) jets.push_back( toL1(jet_pt6,jet_rnd6) );
    if(jet_pt7>0 and fabs(jet_eta7)<eta_threshold) jets.push_back( toL1(jet_pt7,jet_rnd7) );
    std::sort (jets.begin(), jets.end());
    std::reverse (jets.begin(), jets.end());

    if(pos>=0 and pos<jets.size()) return jets.at(pos);
    if(pos>=0 and pos>=jets.size()) return 0.;
    //pos<0 return HT
    float sum=0;
    for (auto& n : jets){
        //if (n<10) continue;
        sum += n;
    }
    return sum;

    }
