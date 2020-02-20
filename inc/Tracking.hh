#ifndef __TRACKING_HH
#define __TRACKING_HH

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "TEnv.h"
#include <algorithm>
#include "Gretina.hh"
#include "TrackSettings.hh"

using namespace std;

#define MAXNIPOINTS 8
#define NPERM MAXNIPOINTS*7*6*5*4*3*2*1

class Tracking {
public:
  Tracking(){
    fset = new TrackSettings();
    GeneratePermutations();
  };
  Tracking(TrackSettings* set){
    fset = set;
    GeneratePermutations();    
  };
  ~Tracking();
  void GeneratePermutations();
  void SetGretina(GretinaCalc *gr){fgr = gr;}
  void SortInClusters();
  vector<double> GetEsum(){return fesum;}
  GretinaEvent* GetEvent(){return fevt;}
  GretinaTrack* GetTrack(int n){return fevt->GetTrack(n);}
  void TrackCluster(GretinaTrack*);
private:
  vector<HitCalc*> SortInCluster(vector<HitCalc*> hits);
  double FigureOfMerit(vector<HitCalc*> hits, int nperm, double esum);
  
  TrackSettings* fset;
  GretinaCalc *fgr;
  GretinaEvent *fevt;
  vector<double> fesum;
  int fperm[MAXNIPOINTS+1][NPERM][MAXNIPOINTS+1];
  int fnperm[MAXNIPOINTS+1];
  int fcurrentsize;
};

#endif
