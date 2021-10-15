#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sys/time.h>
#include <signal.h>

#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"

#include "CommandLineInterface.hh"
#include "SimulatedEvent.hh"
#include "Settings.hh"
//#include "Calibration.hh"

using namespace TMath;
using namespace std;

bool signal_received = false;
void signalhandler(int sig){
  if (sig == SIGINT){
    signal_received = true;
  }
}

double get_time(){
    struct timeval t;
    gettimeofday(&t, NULL);
    double d = t.tv_sec + (double) t.tv_usec/1000000;
    return d;
}
int main(int argc, char* argv[]){
  double time_start = get_time();
  TStopwatch timer;
  timer.Start();
  signal(SIGINT,signalhandler);

  int LastBuffer =-1;
  char *InputFile = NULL;
  char *RootFile = NULL;
  vector<char*> SettingFile;
  int denom = 10000;
  bool wrawtree = false;
  bool wrawhist = false;
  bool wcaltree = true;
  bool wsimtree = true;
  bool wcalhist = false;

  //Read in the command line arguments
  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-lb", "last buffer to be read", &LastBuffer);
  interface->Add("-i", "input file", &InputFile);
  interface->Add("-o", "output file", &RootFile);
  interface->Add("-s", "settingsfile", &SettingFile);
  interface->Add("-rt", "write raw tree", &wrawtree);
  interface->Add("-st", "write sim tree", &wsimtree);
  interface->CheckFlags(argc, argv);

  //Complain about missing mandatory arguments
  if(InputFile == NULL){
    cout << "No input file given " << endl;
    return 1;
  }
  if(RootFile == NULL){
    cout << "No output ROOT file given " << endl;
    return 2;
  }
  if(SettingFile.size() == 0){
    cout << "No settings file given " << endl;
    return 3;
  }

  //Open the input and output files.
  TFile *ofile = new TFile(RootFile,"RECREATE");
  cout<<"input file: "<<InputFile<< endl;
  cout<<"writing to output file: "<<RootFile<< endl;
  cout << "------------------------------------" << endl;
  FILE *infile;
  string infilestr = string(InputFile);
  char extension[20];
  strcpy(extension, infilestr.substr(infilestr.find_last_of(".")+1).c_str());
  if(strcmp(extension,"dat")==0){
    infile = fopen(InputFile,"r");
  } 
  else{
    cout << "Unknown file type. Aborting" << endl;
    ofile->Close();
    return 7;
  }
  if(infile == NULL){
    cout << "Sorry I couldn't find the file: " << InputFile << ". Aborting ..." << endl;
    ofile->Close();
    return 3;
  }

  //Read the settings
  Settings* set = new Settings(SettingFile);
  ofile->cd();
  int vl = set->VLevel();
  if(vl>1)
    set->PrintSettings();
  set->Write("settings",TObject::kOverwrite);
  //Initialize the data structures for the event building.
  int buffers = 0;
  long long int bytes_read = 0;
  struct crys_ips_abcd5678 inbuf_abcd5678[1];
  sim_clust inbuf_hicari[1];
  G4SIM_EGS inbuf_g4sim_emitted_gammas[1];
  ZD_PHYSICSDATA inbuf_zero_physicsdata[1];
  SimulatedEvent *evt = new SimulatedEvent(set);
  evt->SetWrite(wrawtree, wrawhist, wcaltree, wcalhist, wsimtree);
  evt->SetVL(vl);
  Calibration *cal = new Calibration(set);
  evt->SetCalibration(cal);
  evt->Init();

  if(vl>0)
    cout << "start " << endl;


  int error = 0;
  //Loop over the entirety of the input file.
  while(!feof(infile) && !signal_received){

    //Finish reading if you have read as many buffers as the user has requested.
    if(LastBuffer > 0 && buffers >= LastBuffer)
      break;

    size_t bsize;
    //Read the header of the new buffer.
    unsigned int header[4];
    bsize = fread(header, sizeof(unsigned int), 4, infile);
    if(bsize==0){
      cout << "cannot read "<< InputFile << " (anymore)" << endl;
      break;
    }
    long long int ts = (long long int)header[3] << 32;
    ts += header[2];



    if(vl>1){
      for(int i=0;i<4;i++)
	cout <<(hex)<< header[i] <<"\t"<<(dec)<< header[i] <<endl;
    }
    bytes_read += 4*sizeof(unsigned int);

    //Decode the event differently depending on which type of data is identified in the header.
    //For each, pass the data from the file into the UnpackedEvent to be read.
    if(header[0]==GRETINA_ID){
      if(vl>0)
	cout << "gret timestamp:\t"<< ts << "\tlength: "<< header[1] << "\thex ts: " <<(hex) << ts << "\tlength: "<< header[1] <<(dec)<< endl;
      Crystal* crys;
      if (header[1]==sizeof(crys_ips_abcd5678)){
	bsize = fread(&inbuf_abcd5678[0], sizeof(crys_ips_abcd5678), 1, infile);
	bytes_read += sizeof(struct crys_ips_abcd5678);
	crys = new Crystal(inbuf_abcd5678[0]);
      } else {
	cout << "Unknown size for mode2 data" << endl;
	break;
      }
      if(bsize>0)
	error = evt->DecodeGretina(crys,ts);
      if(error){
	cout << "An error ("<<error<<") occured in DecodeGretina() while processing file: " << InputFile << ". Continuing ..." << endl;
	continue;
      }
    }
    else if(header[0]==HICARI_ID){
      if(vl>0)
	cout << "hicari timestamp:\t"<< ts << "\tlength: "<< header[1] << "\thex ts: " <<(hex) << ts << "\tlength: "<< header[1] <<(dec)<< endl;
      HiCARIHit* crys;
      if(header[1]==sizeof(sim_clust)){
	bsize = fread(&inbuf_hicari[0], sizeof(sim_clust), 1, infile);
	bytes_read += sizeof(struct sim_clust);
	crys = new HiCARIHit(inbuf_hicari[0],ts);
      }
      else{
	cout << "Unknown size for HiCARI data" << endl;
	break;
      }
      if(bsize>0)
	error = evt->DecodeHiCARI(crys,ts);
      if(error){
	cout << "An error ("<<error<<") occured in DecodeHiCARI() while processing file: " << InputFile << ". Continuing ..." << endl;
	continue;
      }
    }
    else if(header[0]==ZERO_PHYSDATA_ID){
      if(vl>0)
	cout << "zerodeg timestamp:\t"<< ts << "\tlength: "<< header[1] << "\thex: " <<(hex) << ts << "\tlength: "<< header[1] <<(dec) << "\tsizeof(ZERO_PHYSICSDATA): " << sizeof(ZD_PHYSICSDATA) << endl;
      if(header[1]==sizeof(ZD_PHYSICSDATA)){
	bsize = fread(&inbuf_zero_physicsdata[0], sizeof(ZD_PHYSICSDATA), 1, infile);
	bytes_read += sizeof(ZD_PHYSICSDATA);
      }
      else{
	cout << "Unknown size for ZD_PHYSICSDATA event" << endl;
	break;
      }
      if(bsize>0)
	error = evt->DecodeZeroDegPhysicsData(&inbuf_zero_physicsdata[0], ts);
      if(error){
	cout << "An error ("<<error<<") occured in DecodeZeroDegPhysicsData() while processing file: " << InputFile << ". Continuing ..." << endl;
	continue;
      }

    }
    else if(header[0]==GAMMA_G4SIM_ID){
      if(vl>0)
	cout << "gsim timestamp:\t"<< ts << "\tlength: "<< header[1] << "\thex: " <<(hex) << ts << "\tlength: "<< header[1] <<(dec) << "\tsizeof(G4SIM_EGS): " << sizeof(G4SIM_EGS) << endl;
      if (header[1]==sizeof(G4SIM_EGS)){
	bsize = fread(&inbuf_g4sim_emitted_gammas[0], sizeof(G4SIM_EGS), 1, infile);
	bytes_read += sizeof(G4SIM_EGS);
      } else {
	cout << "Unknown size for G4SIM_EGS event" << endl;
	break;
      }
      if(bsize>0)
	error = evt->DecodeGammaG4Sim(&inbuf_g4sim_emitted_gammas[0], ts);
      if(error){
	cout << "An error ("<<error<<") occured in DecodeGretinaG4Sim while processing file: " << InputFile << ". Continuing ..." << endl;
	continue;
      }
    }
    else{
      cout << "unidentified header " << header[0] << "\thex: " <<(hex) << header[0] <<(dec)<< endl;
      cout << "unidentified timestamp:\t"<< ts << "\tlength: "<< header[1] << "\thex: " <<(hex) << ts << "\tlength: "<< header[1] <<(dec)<< endl;
      cout << "at buffer nr " << buffers << endl;
      if(buffers==0)
	continue;
      break;
    }

    //Write the trees out to disk every denom events.
    buffers++;
    if(buffers % denom == 0){
      if(wcaltree)
	evt->GetCalTree()->AutoSave();
      if(wsimtree)
	evt->GetSimTree()->AutoSave();
      double time_end = get_time();
      cout << "\r" << buffers << " buffers read... "<<bytes_read/(1024*1024)<<" MB... "<<buffers/(time_end - time_start) << " buffers/s" << flush;
    }
  }

  //Finish reading the last event and close it out.
  evt->WriteLastEvent();

  cout << "Total of " << buffers << " data buffers ("<<bytes_read/(1024*1024)<<" MB) and" << endl;

  if(wcaltree){
    cout << "Total of " << evt->NrOfCalEvents() << " cal events ("<<evt->GetCalTree()->GetZipBytes()/(1024*1024)<<" MB) written."  << endl;
  }
  if(wsimtree){
    cout << "Total of " << evt->NrOfSimEvents() << " sim events ("<<evt->GetSimTree()->GetZipBytes()/(1024*1024)<<" MB) written."  << endl;
  }
  double time_end = get_time();
  cout << "Program Run time " << time_end - time_start << " s." << endl;
  cout << "Unpacked " << buffers/(time_end - time_start) << " buffers/s." << endl;
  timer.Stop();
  cout << "\n CPU time: " << timer.CpuTime() << "\tReal time: " << timer.RealTime() << endl;

  ofile->cd();
  //Final cleanup and writing of files.
  if(wcaltree){
    evt->GetCalTree()->Write("",TObject::kOverwrite);
  }
  if(wrawtree){
    evt->GetTree()->Write("",TObject::kOverwrite);
  }
  if(wsimtree){
    evt->GetSimTree()->Write("",TObject::kOverwrite);
  }
  //Cannot simply close ofile, even though it was the one opened at the beginning.
  //If the tree spills out into a new file after 1.8 GB, it closes the original file and opens a new one.
  //ofile then has a dead pointer.
  //However, X->GetCurrentFile() still points to the currently open file.
  //  X->GetCurrentFile()->Close();
  if(wcaltree)
    evt->GetCalTree()->GetCurrentFile()->Close();
  else if(wsimtree)
    evt->GetSimTree()->GetCurrentFile()->Close();

  return 0;
}
