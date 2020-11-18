#include "TLorentzVector.h"

float getResiduals( float gen_pt = 0., float gen_eta = 0., float gen_phi = 0., float gen_mass = 0., 
                    float j0_pt = 0, float j0_eta = 0, float j0_phi = 0, float j0_mass = 0, 
                    float j1_pt = 0, float j1_eta = 0, float j1_phi = 0, float j1_mass = 0, 
                    float j2_pt = 0, float j2_eta = 0, float j2_phi = 0, float j2_mass = 0, 
                    float j3_pt = 0, float j3_eta = 0, float j3_phi = 0, float j3_mass = 0, 
                    float j4_pt = 0, float j4_eta = 0, float j4_phi = 0, float j4_mass = 0, 
                    float j5_pt = 0, float j5_eta = 0, float j5_phi = 0, float j5_mass = 0, 
                    float j6_pt = 0, float j6_eta = 0, float j6_phi = 0, float j6_mass = 0, 
                    float j7_pt = 0, float j7_eta = 0, float j7_phi = 0, float j7_mass = 0, 
                    float j8_pt = 0, float j8_eta = 0, float j8_phi = 0, float j8_mass = 0, 
                    float j9_pt = 0, float j9_eta = 0, float j9_phi = 0, float j9_mass = 0 ){

	TLorentzVector gen;
	TLorentzVector j0, j1, j2, j3, j4, j5, j6, j7, j8, j9;
	gen.SetPtEtaPhiM(gen_pt, gen_eta, gen_phi, gen_mass);
	j0.SetPtEtaPhiM(j0_pt, j0_eta, j0_phi, j0_mass);
	j1.SetPtEtaPhiM(j1_pt, j1_eta, j1_phi, j1_mass);
	j2.SetPtEtaPhiM(j2_pt, j2_eta, j2_phi, j2_mass);
	j3.SetPtEtaPhiM(j3_pt, j3_eta, j3_phi, j3_mass);
	j4.SetPtEtaPhiM(j4_pt, j4_eta, j4_phi, j4_mass);
	j5.SetPtEtaPhiM(j5_pt, j5_eta, j5_phi, j5_mass);
	j6.SetPtEtaPhiM(j6_pt, j6_eta, j6_phi, j6_mass);
	j7.SetPtEtaPhiM(j7_pt, j7_eta, j7_phi, j7_mass);
	j8.SetPtEtaPhiM(j8_pt, j8_eta, j8_phi, j8_mass);
	j9.SetPtEtaPhiM(j9_pt, j9_eta, j9_phi, j9_mass);
    
    const int N=10;
    TLorentzVector jets[N] = {j0,j1,j2,j3,j4,j5,j6,j7,j8,j9};
    float mindr = 9;
    int mini = -1;
    float dr = 9;
    for(int i=0; i<N;i++){
        if(jets[i].E() <= 0) break;
        dr =gen.DeltaR(jets[i]);
        if (dr < mindr){
            mindr = dr;
            mini = i;
        }
    }
    if (mini<0) return 999999;
    
    return jets[mini].Pt() - gen.Pt();

}
