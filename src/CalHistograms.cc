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

void CalHistograms::FillHistograms(HiCARICalc* hi){
  Fill("h_mult",30,0,30,hi->GetMult());
  
  long long int firstTS = -1;
  long long int BRTS = -1;
  
  for(int i=0; i<hi->GetMult(); i++){
    int segs = 6;
    HiCARIHitCalc* hit = hi->GetHit(i);
    if(hit->IsTracking())
      segs = 40;
    if(hit->IsSuperClo())
      segs = 4;
    if(hit->IsBigRIPS()){
      BRTS = hit->GetTS();
      continue;
    }
    if(i==0)
      firstTS = hit->GetTS();
    else
      Fill("h_tsdiff",500,0,1000,hit->GetTS()-firstTS);

    Fill("h_cluster",12,0,12,hit->GetCluster());
    Fill("h_crystal",4,0,4,hit->GetCrystal());
    Fill("h_crystal_vs_cluster",12,0,12,hit->GetCluster(),4,0,4,hit->GetCrystal());
    Fill("h_en_summary",48,0,48,hit->GetCluster()*4+hit->GetCrystal(),4000,0,4000,hit->GetEnergy());
    Fill("h_en_summary_fine",48,0,48,hit->GetCluster()*4+hit->GetCrystal(),16000,0,4000,hit->GetEnergy());
    Fill(Form("h_en_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),16000,0,4000,hit->GetEnergy());
    Fill("h_segsum_vs_en",2000,0,4000,hit->GetEnergy(),2000,0,4000,hit->GetSegSum());
    Fill(Form("h_segsum_vs_en_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),1000,0,4000,hit->GetEnergy(),1000,0,4000,hit->GetSegSum());
    Fill(Form("h_segmult_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),segs+1,0,segs+1,hit->GetSegmentNr().size());
    
    for(UShort_t j=0; j<hit->GetSegmentNr().size(); j++){
      Fill(Form("h_segen_vs_nr_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),segs,0,segs,hit->GetSegmentNr().at(j),4000,0,4000,hit->GetSegmentEn().at(j));
    }
    int clu = hit->GetCluster();
    int cry = hit->GetCrystal();
    
    for(int j=i+1; j<hi->GetMult(); j++){
      HiCARIHitCalc* sec = hi->GetHit(j);
      Fill("h_engg_symm",4000,0,4000,hit->GetEnergy(),4000,0,4000,sec->GetEnergy());
      Fill("h_engg_symm",4000,0,4000,sec->GetEnergy(),4000,0,4000,hit->GetEnergy());
      if(sec->GetEnergy() > hit->GetEnergy()){
    	Fill("h_hengg_tsdiff",500,-500,500,sec->GetTS()-hit->GetTS(),4000,0,4000,sec->GetEnergy());
    	Fill("h_lengg_tsdiff",500,-500,500,sec->GetTS()-hit->GetTS(),4000,0,4000,hit->GetEnergy());
      }
      else{
    	Fill("h_hengg_tsdiff",500,-500,500,hit->GetTS()-sec->GetTS(),4000,0,4000,hit->GetEnergy());
    	Fill("h_lengg_tsdiff",500,-500,500,hit->GetTS()-sec->GetTS(),4000,0,4000,sec->GetEnergy());
      }
      if(fabs(hit->GetTS()-sec->GetTS()) < fSett->CoincTimeDiff()){
    	Fill("h_engg_symm_TC",4000,0,4000,hit->GetEnergy(),4000,0,4000,sec->GetEnergy());
    	Fill("h_engg_symm_TC",4000,0,4000,sec->GetEnergy(),4000,0,4000,hit->GetEnergy());
      }
    }
    // for(int j=0; j<hi->GetMult(); j++){
    //   HiCARIHitCalc* sec = hi->GetHit(j);
    //   if(sec->GetCluster() == clu && sec->GetCrystal() == cry)
    // 	continue;
    //   //Fill(Form("h_engg_clus%02d_crys%02d",clu,cry),4000,0,4000,hit->GetEnergy(),4000,0,4000,sec->GetEnergy());
    //   if(fabs(hit->GetTS()-sec->GetTS()) < fSett->CoincTimeDiff()){
    // 	Fill(Form("h_engg_TC_clus%02d_crys%02d",clu,cry),4000,0,4000,hit->GetEnergy(),4000,0,4000,sec->GetEnergy());
    //   }      
    // }
  }//hit mult
  Fill("h_mult_multAB",30,0,30,hi->GetMult(),30,0,30,hi->GetMultAB());
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
    Fill("hAB_en_summary_fine",48,0,48,hit->GetCluster()*4+hit->GetCrystal(),16000,0,4000,hit->GetEnergy());
    Fill(Form("hAB_en_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),16000,0,4000,hit->GetEnergy());
    Fill("hAB_segsum_vs_en",2000,0,4000,hit->GetEnergy(),2000,0,4000,hit->GetSegSum());
    Fill(Form("hAB_segsum_vs_en_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),1000,0,4000,hit->GetEnergy(),1000,0,4000,hit->GetSegSum());
    Fill(Form("hAB_segmult_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),segs+1,0,segs+1,hit->GetSegmentNr().size());
    
    for(UShort_t j=0; j<hit->GetSegmentNr().size(); j++){
      Fill(Form("hAB_segen_vs_nr_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),segs,0,segs,hit->GetSegmentNr().at(j),4000,0,4000,hit->GetSegmentEn().at(j));
    }
    int clu = hit->GetCluster();
    int cry = hit->GetCrystal();
    for(int j=i+1; j<hi->GetMultAB(); j++){
      HiCARIHitCalc* sec = hi->GetHitAB(j);
      Fill("hAB_engg_symm",1000,0,4000,hit->GetEnergy(),1000,0,4000,sec->GetEnergy());
      Fill("hAB_engg_symm",1000,0,4000,sec->GetEnergy(),1000,0,4000,hit->GetEnergy());
      if(sec->GetEnergy() > hit->GetEnergy()){
    	Fill("hAB_hengg_tsdiff",500,-500,500,sec->GetTS()-hit->GetTS(),4000,0,4000,sec->GetEnergy());
    	Fill("hAB_lengg_tsdiff",500,-500,500,sec->GetTS()-hit->GetTS(),4000,0,4000,hit->GetEnergy());
      }
      else{
    	Fill("hAB_hengg_tsdiff",500,-500,500,hit->GetTS()-sec->GetTS(),4000,0,4000,hit->GetEnergy());
    	Fill("hAB_lengg_tsdiff",500,-500,500,hit->GetTS()-sec->GetTS(),4000,0,4000,sec->GetEnergy());
      }
      if(fabs(hit->GetTS()-sec->GetTS()) < fSett->CoincTimeDiff()){
    	Fill("hAB_engg_symm_TC",4000,0,4000,hit->GetEnergy(),4000,0,4000,sec->GetEnergy());
    	Fill("hAB_engg_symm_TC",4000,0,4000,sec->GetEnergy(),4000,0,4000,hit->GetEnergy());
      }
    }
    // for(int j=0; j<hi->GetMultAB(); j++){
    //   HiCARIHitCalc* sec = hi->GetHitAB(j);
    //   if(sec->GetCluster() == clu && sec->GetCrystal() == cry)
    // 	continue;
    //   //Fill(Form("hAB_engg_clus%02d_crys%02d",clu,cry),4000,0,4000,hit->GetEnergy(),4000,0,4000,sec->GetEnergy());
    //   if(fabs(hit->GetTS()-sec->GetTS()) < fSett->CoincTimeDiff()){
    // 	Fill(Form("hAB_engg_TC_clus%02d_crys%02d",clu,cry),4000,0,4000,hit->GetEnergy(),4000,0,4000,sec->GetEnergy());
    //   }      
    // }

  }//hit mult

  for(int i=0; i<hi->GetMult(); i++){
    HiCARIHitCalc* hit = hi->GetHit(i);
    if(hit->IsBigRIPS()){
      continue;
    }
    if(BRTS>0){
      Fill("h_TS_BRTS",2000,-1000,1000,hit->GetTS() - BRTS );
      Fill("h_TS_BRTS_egam",2000,-1000,1000,hit->GetTS() - BRTS ,2000,0,2000,hit->GetEnergy());
      Fill("h_TS_BRTS_egamDC",2000,-1000,1000,hit->GetTS() - BRTS ,2000,0,2000,hit->GetDCEnergy());
    }
  }


  
}
void CalHistograms::FillHistograms(GretinaCalc* gr){
  Fill("g_mult",30,0,30,gr->GetMult());
  long long int firstTS = -1;
  for(UShort_t g=0; g<gr->GetMult(); g++){
    HitCalc* hit = gr->GetHit(g);
    if(g==0)
      firstTS = hit->GetTS();
    else
      Fill("g_tsdiff",500,0,1000,hit->GetTS()-firstTS);

    Fill("g_cluster",20,0,20,hit->GetCluster());
    Fill("g_crystal",4,0,4,hit->GetCrystal());
    Fill("g_crystal_vs_cluster",20,0,20,hit->GetCluster(),4,0,4,hit->GetCrystal());
    Fill("g_en_summary",80,0,80,hit->GetCluster()*4+hit->GetCrystal(),4000,0,4000,hit->GetEnergy());
    Fill("g_en_summary_fine",80,0,80,hit->GetCluster()*4+hit->GetCrystal(),16000,0,4000,hit->GetEnergy());
    Fill(Form("g_en_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),16000,0,4000,hit->GetEnergy());
    Fill("g_maxhit_vs_en",2000,0,4000,hit->GetEnergy(),2000,0,4000,hit->GetMaxSingleHit());
    Fill(Form("g_maxhit_vs_en_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),1000,0,4000,hit->GetEnergy(),1000,0,4000,hit->GetMaxSingleHit());
    
  }//hits
  Fill("g_mult_multAB",30,0,30,gr->GetMult(),30,0,30,gr->GetMultAB());
  Fill("gAB_mult",30,0,30,gr->GetMultAB());
  for(UShort_t g=0; g<gr->GetMultAB(); g++){
    HitCalc* hit = gr->GetHit(g);
 
    Fill("gAB_cluster",20,0,20,hit->GetCluster());
    Fill("gAB_crystal",4,0,4,hit->GetCrystal());
    Fill("gAB_crystal_vs_cluster",20,0,20,hit->GetCluster(),4,0,4,hit->GetCrystal());
    Fill("gAB_en_summary",80,0,80,hit->GetCluster()*4+hit->GetCrystal(),4000,0,4000,hit->GetEnergy());
    Fill("gAB_en_summary_fine",80,0,80,hit->GetCluster()*4+hit->GetCrystal(),16000,0,4000,hit->GetEnergy());
    Fill(Form("gAB_en_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),16000,0,4000,hit->GetEnergy());
    Fill("gAB_maxhit_vs_en",2000,0,4000,hit->GetEnergy(),2000,0,4000,hit->GetMaxSingleHit());
    Fill(Form("gAB_maxhit_vs_en_clus%02d_crys%02d",hit->GetCluster(),hit->GetCrystal()),1000,0,4000,hit->GetEnergy(),1000,0,4000,hit->GetMaxSingleHit());
    
  }//hits

}
