#ifndef __TRACKSETTINGS_HH
#define __TRACKSETTINGS_HH

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "Settings.hh"

#include "TSystem.h"
#include "TEnv.h"

using namespace std;

class TrackSettings : public Settings {
public:
  TrackSettings();
  TrackSettings(const char* filename);
  TrackSettings(vector<char*> files);
  ~TrackSettings();
  void ReadTrackSettings(TEnv* set);
  void PrintTrackSettings();
  int MaxBad(){
    return fmaxBad;
  }
  bool RedoMap(){
    return fredoMap;
  }
  double JumpFOM(){
    return fjumpFOM;
  }
  int Probabilities(){
    return fprob;
  }
  double RadiusGE(){
    return fradiusGE;
  }
protected:
  float fjumpFOM;				
  int fmaxBad;
  int fredoMap;
  int fprob;
  float fradiusGE;				

  ClassDef(TrackSettings, 1)
};

#endif
