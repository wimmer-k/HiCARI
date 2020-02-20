#ifndef SIM_HISTOGRAMS_HH__
#define SIM_HISTOGRAMS_HH__

#include <iostream>
#include <iomanip>

#include "TObject.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TEnv.h"
#include "math.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdexcept>

#include "Settings.hh"
#include "Gretina.hh"
#include "Miniball.hh"
#include "ZeroDeg.hh"
#include "GammaSim.hh"
#ifdef USELISA
#include "LISA.hh"
#endif

using namespace std;

class SimHistograms {
public:
  SimHistograms(Settings* set = NULL, int nentries=100000){
    Init(set, nentries);
  }
  ~SimHistograms(){
    delete fhlist;
  }
  void Init(Settings* set, int nentries){
    if (set!=NULL){
      fSett = set;
    } else {
      fSett = new Settings;
    }
    fhlist = new TList;
    fnentries = nentries;
    fentry =0;
  }
  void AddEntry(){
    fentry++;
  }
  TList* GetHList(){return fhlist;}
  void Write(){
    fhlist->Sort();
    fhlist->Write();
  }

#ifdef USELISA
  void FillHistograms(GretinaCalc* gr, MiniballCalc* mb, ZeroDeg* zd, GammaSim* gs, LISA* li);
#else
  void FillHistograms(GretinaCalc* gr, MiniballCalc* mb, ZeroDeg* zd, GammaSim* gs);
#endif

  void Fill(string name,int bins, double low, double high, double value){
    try{
      fhmap.at(name)->Fill(value);
    } catch(const out_of_range& e) {
      cout << "New 1-d histogram: " << name << endl;
      TH1F* newHist = new TH1F(name.c_str(),name.c_str(),bins,low,high);
      newHist->Fill(value);
      fhlist->Add(newHist);
      fhmap[name] = newHist;
    }
  }
  void Fill(string name,
	    int binsX, double lowX, double highX, double valueX,
	    int binsY, double lowY, double highY, double valueY){
    try{
      fhmap.at(name)->Fill(valueX,valueY);
    } catch(const out_of_range& e) {
      cout << "New 2-d histogram: " << name << endl;
      TH2F* newHist = new TH2F(name.c_str(),name.c_str(),
			       binsX,lowX,highX,
			       binsY,lowY,highY);
      newHist->Fill(valueX,valueY);
      fhlist->Add(newHist);
      fhmap[name] = newHist;
    }
  }

protected:
  TList* fhlist;
  map<string,TH1*> fhmap;
  int fnentries;
  int fentry;
  Settings* fSett;

};

#endif
