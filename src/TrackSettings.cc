#include "TrackSettings.hh"

using namespace std;
TrackSettings::TrackSettings()
  : Settings("~/analysis/settings/nocal.dat"){
  TrackSettings("~/analysis/settings/nocal.dat");
}
TrackSettings::TrackSettings(const char* filename)
  : Settings(filename){
  TEnv set(filename);
  ReadTrackSettings(&set);
  if(fVerboseLevel>-1)
    PrintTrackSettings();
}
TrackSettings::TrackSettings(vector<char*> files): Settings(files){
  TEnv set;
  //Reverse order, because duplicate entries are ignored, instead of overwriting.
  for(vector<char*>::reverse_iterator it = files.rbegin(); it!=files.rend(); it++){
    set.ReadFile(*it,kEnvLocal);
    fInputFiles.push_back(*it);
  }
  ReadTrackSettings(&set);
  if(fVerboseLevel>1)
    PrintTrackSettings();
}
TrackSettings::~TrackSettings(){}
void TrackSettings::ReadTrackSettings(TEnv* set){
  fmaxBad= set->GetValue("Max.Bad.FOM",3);
  fjumpFOM = set->GetValue("Jump.FOM",0.8);
  fredoMap = set->GetValue("Redo.Map",0);
  fprob = set->GetValue("Probabilities",0);
  fradiusGE = set->GetValue("Gretina.Radius",185.0);
}
void TrackSettings::PrintTrackSettings(){
  cout << "Max.Bad.FOM\t" << fmaxBad << endl;
  cout << "Jump.FOM\t" << fjumpFOM << endl;
  cout << "Redo.Map\t" << fredoMap << endl;
  cout << "Probabilities\t" << fprob << endl;
  cout << "Gretina.Radius\t" << fradiusGE << endl;
  //cout << "\t" << f << endl;
  //cout << "\t" << f << endl;
  //cout << "\t" << f << endl;
 
}
