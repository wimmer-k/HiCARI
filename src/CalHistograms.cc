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
}
#endif
