#include "MergeHistograms.hh"

#include <iostream>
#include <iomanip>
#include <string.h>
#include <sys/time.h>

#include <sstream>


#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2S.h"
#include "TH1S.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"

using namespace TMath;
using namespace std;

void MergeHistograms::Write(){
  fhlist->Sort();
  fhlist->Write();

}

void MergeHistograms::FillHistograms(int checkADC, HiCARICalc* hi, long long int brTS, long long int hiTS){
  fentry++;
  //Determine which of the systems are present in the data.
  bool hashicari = hi->GetMult()!=0;
  bool hasbigrips = true;//checkADC>0;

  if(hashicari && hasbigrips)
    FillCorrelationHistograms(checkADC, hi, brTS,hiTS);
}
void MergeHistograms::FillCorrelationHistograms(int checkADC, HiCARICalc* hi, long long int brTS, long long int hiTS){
  Fill("hTSdiff",1000,-500,500,brTS-hiTS);
  for(UShort_t i=0; i<hi->GetMult();i++){
    HiCARIHitCalc* hit = hi->GetHit(i);
    if(hit->GetCluster()==3 && hit->GetCrystal()==0)
      Fill("hEcorr",4096,0,4096,checkADC, 4000,0,4000,hit->GetEnergy());
  }
      
}


