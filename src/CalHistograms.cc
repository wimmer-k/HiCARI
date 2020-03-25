#include <iostream>
#include <iomanip>
#include <string.h>
#include <sys/time.h>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2S.h"
#include "TH1S.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "TCutG.h"
#include "TEnv.h"
#include "TKey.h"
#include "TDirectory.h"

#include "CalHistograms.hh"
using namespace TMath;
using namespace std;

//not nice but works at least
static TCutG* TimeCut;
static bool foundTCut = false;

#ifdef SIMULATION
void CalHistograms::FillHistograms(GretinaCalc* gr, MiniballCalc* mb, ZeroDeg* zd, MINOS* mi){

  bool hasgret = gr->GetMult()>0;
  
  for (UShort_t g=0; g<gr->GetMult(); g++){
    HitCalc* hit = gr->GetHit(g);
    float energy = hit->GetEnergy();
    float energy_dc = hit->GetDCEnergy();
   Fill("egam",
	 8000,0,8000,energy);
    Fill("egam_tgam",
	 2000,0,4000,energy,
	 400,-200,200,hit->GetTime());
    Fill("egamdc",
	 8000,0,8000,energy_dc);
    Fill("egamdc_tgam",
	 2000,0,4000,energy_dc,
	 400,-200,200,hit->GetTime());
    Fill("egam_summary",
	 32,-0.5,31.5,4*fSett->Clu2Det(hit->GetCluster())+hit->GetCrystal(),
	 2000,0,4000,energy);
   }
  for(UShort_t g=0;g<gr->GetMultAB();g++){
    HitCalc* hit = gr->GetHitAB(g);
    float energy = hit->GetEnergy();
    float energy_dc = hit->GetDCEnergy();
    Fill("egamAB",
	 8000,0,8000,energy);
    Fill("egamAB_tgam",
	 2000,0,4000,energy,
	 400,-200,200,hit->GetTime());
    Fill("egamABdc",
	 8000,0,8000,energy_dc);
    Fill("egamABdc_tgam",
	 2000,0,4000,energy_dc,
	 400,-200,200,hit->GetTime());
    Fill("egamAB_summary",
	 32,-0.5,31.5,4*fSett->Clu2Det(hit->GetCluster())+hit->GetCrystal(),
	 2000,0,4000,energy);
   }

}
#else
void CalHistograms::FillHistograms(HiCARICalc* hi){
  Fill("h_mult",30,0,30,hi->GetMult());
  for(int i=0; i<hi->GetMult(); i++){
    int segs = 6;
    HiCARIHitCalc* hit = hi->GetHit(i);
    if(hit->IsTracking())
      segs = 40;
    if(hit->IsSuperClo())
      segs = 8;
    if(hit->IsBigRIPS())
      continue;

    Fill("h_cluster",12,0,12,hit->GetCluster());
    Fill("h_crystal",4,0,4,hit->GetCrystal());
    Fill("h_crystal_vs_cluster",12,0,12,hit->GetCluster(),4,0,4,hit->GetCrystal());
    Fill("h_en_summary",48,0,48,hit->GetCluster()*4+hit->GetCrystal(),4000,0,4000,hit->GetEnergy());
    Fill(Form("h_en_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),4000,0,4000,hit->GetEnergy());
    Fill("h_segsum_vs_en",2000,0,4000,hit->GetEnergy(),2000,0,4000,hit->GetSegSum());
    Fill(Form("h_segsum_vs_en_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),2000,0,4000,hit->GetEnergy(),2000,0,4000,hit->GetSegSum());
    Fill(Form("h_segmult_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),segs,0,segs,hit->GetSegmentNr().size());
    
    for(UShort_t j=0; j<hit->GetSegmentNr().size(); j++){
      Fill(Form("h_segen_vs_nr_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),segs,0,segs,hit->GetSegmentNr().at(j),4000,0,4000,hit->GetSegmentEn().at(j));
    }
    for(int j=i+1; j<hi->GetMult(); j++){
      HiCARIHitCalc* sec = hi->GetHit(j);
      Fill("h_engg_symm",4000,0,4000,hit->GetEnergy(),4000,0,4000,sec->GetEnergy());
      Fill("h_engg_symm",4000,0,4000,sec->GetEnergy(),4000,0,4000,hit->GetEnergy());
    }
  }//hit mult
  Fill("hAB_mult",30,0,30,hi->GetMultAB());
  for(int i=0; i<hi->GetMultAB(); i++){
    int segs = 6;
    HiCARIHitCalc* hit = hi->GetHitAB(i);
    if(hit->IsTracking())
      segs = 40;
    if(hit->IsSuperClo())
      segs = 8;
    if(hit->IsBigRIPS())
      continue;

    Fill("hAB_cluster",12,0,12,hit->GetCluster());
    Fill("hAB_crystal",4,0,4,hit->GetCrystal());
    Fill("hAB_crystal_vs_cluster",12,0,12,hit->GetCluster(),4,0,4,hit->GetCrystal());
    Fill("hAB_en_summary",48,0,48,hit->GetCluster()*4+hit->GetCrystal(),4000,0,4000,hit->GetEnergy());
    Fill(Form("hAB_en_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),4000,0,4000,hit->GetEnergy());
    Fill("hAB_segsum_vs_en",2000,0,4000,hit->GetEnergy(),2000,0,4000,hit->GetSegSum());
    Fill(Form("hAB_segsum_vs_en_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),2000,0,4000,hit->GetEnergy(),2000,0,4000,hit->GetSegSum());
    Fill(Form("hAB_segmult_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),segs,0,segs,hit->GetSegmentNr().size());
    
    for(UShort_t j=0; j<hit->GetSegmentNr().size(); j++){
      Fill(Form("hAB_segen_vs_nr_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),segs,0,segs,hit->GetSegmentNr().at(j),4000,0,4000,hit->GetSegmentEn().at(j));
    }
    for(int j=i+1; j<hi->GetMultAB(); j++){
      HiCARIHitCalc* sec = hi->GetHitAB(j);
      Fill("hAB_engg_symm",4000,0,4000,hit->GetEnergy(),4000,0,4000,sec->GetEnergy());
      Fill("hAB_engg_symm",4000,0,4000,sec->GetEnergy(),4000,0,4000,hit->GetEnergy());
    }
  }//hit mult

}
#endif
