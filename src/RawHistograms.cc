#include "RawHistograms.hh"

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

#include "Gretina.hh"

using namespace TMath;
using namespace std;

void RawHistograms::Write(){
  fhlist->Sort();
  fhlist->Write();

}

void RawHistograms::FillHistograms(Gretina* gr, Miniball* mb, ZeroDeg* zd, MINOS* mi){
  fentry++;
  //Determine which of the systems are present in the data.
  bool hasmode2 = gr->GetMult()!=0;


  if(hasmode2){
    FillMode2Histograms(gr);
  }
}

void RawHistograms::FillMode2Histograms(Gretina* gr){
  //Mode 2 histograms

  Fill("hmult",30,0,30,gr->GetMult());

  for(int i=0; i<gr->GetMult(); i++){
    Crystal* cr = gr->GetHit(i);
    Fill("hgamma",
	 10000,0,10000,cr->GetEnergy());
    Fill("hsegsum",
	 10000,0,10000,cr->GetSegmentSum());
    Fill("hgamma_segsum",
	 10000,0,10000,cr->GetEnergy(),
	 10000,0,10000,cr->GetSegmentSum());
    Fill("hgamma_segsum_diff",
	 20000,-10000,10000,cr->GetEnergy()-cr->GetSegmentSum());
    Fill("hgamma_IPsum",
	 10000,0,10000,cr->GetEnergy(),
	 10000,0,10000,cr->GetIPSum());
    Fill("hgamma_IPsum_diff",
	 20000,-10000,10000,cr->GetEnergy()-cr->GetIPSum());
    if (cr->GetEnergy()>300 && cr->GetEnergy()<10000){
      Fill("hgamma_segsum_diff_gated",
	   20000,-10000,10000,cr->GetEnergy()-cr->GetSegmentSum());
    }
    Fill("hm2_segmult",
	 41,-0.5,40.5,cr->GetMult());
    Fill("herror",
	 10,-0.5,9.5,cr->GetError());
    Fill("herror_ID",
	 150,-0.5,149.5,cr->GetID(),
	 10,-0.5,9.5,cr->GetError());
    Fill("hgamma_error",
	 10,-0.5,9.5,cr->GetError(),
	 500,0,10000,cr->GetEnergy());
    Fill(Form("hgamma_d%d_c%d",fSett->Clu2Det(cr->GetCluster()),cr->GetCrystal()),
	 3500,0,3500,cr->GetEnergy());
    Fill("hhitpattern",
	 150,-0.5,149.5,cr->GetID());
    Fill("hhitpattern_en",
	 150,-0.5,149.5,cr->GetID(),
	 2000,0,2000,cr->GetEnergy());

    for(int i=0; i<4; i++){
      Fill(Form("hcore_e_%d",i),
	   10000,0,3e4,cr->GetCoreE(i));
    }
    Fill("hprestep",
	 1000,0,1000,cr->GetPreStep());
    Fill("hpoststep",
	 1000,0,1000,cr->GetPostStep());

    for(int j=0; j<cr->GetMult(); j++){
      IPoint* IP = cr->GetIPoint(j);
      int segID = cr->GetID()*36 + IP->GetSeg();
      Fill("hseghitpattern",
	   30*4*36,-0.5,30*4*36-0.5,segID);
      Fill("hseghitpattern_en",
	   30*4*36,-0.5,30*4*36-0.5,segID,
	   2000,0,2000,IP->GetSegEnergy());
    }
  }
}
