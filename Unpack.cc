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
#include "UnpackedEvent.hh"
#include "RunInfo.hh"
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
  int denom = 100000;
  int wrawtree = 1;
  int wrawhist = 0;
  int wcaltree = 0;
  int wcalhist = 0;
  int makemode2 = 1;
  bool noHFC = false;
  
  //Read in the command line arguments
  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-lb", "last buffer to be read", &LastBuffer);
  interface->Add("-i", "input file", &InputFile);
  interface->Add("-o", "output file", &RootFile);
  interface->Add("-s", "settingsfile", &SettingFile);
  interface->Add("-rt", "write raw tree", &wrawtree);
  interface->Add("-ct", "write cal tree", &wcaltree);
  interface->Add("-rh", "write raw histos", &wrawhist);
  interface->Add("-ch", "write cal histos", &wcalhist);
  interface->Add("-m", "make mode2 data, no decomp", &makemode2);
  interface->Add("--no-HFC","do not use HFC as an intermediate step",&noHFC);
  interface->CheckFlags(argc, argv);

  if(wrawtree == 0 && wcaltree == 0)
    denom*=10;
  if(wcaltree || wcalhist)
    makemode2 = 1;
  
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

  //get the run number from the filename
  int run;
  TString ifname(InputFile);
  ifname.Remove(0,ifname.Length()-21); // Last 21 characters: RunXXXX/GlobalRaw.dat
  //cout << "ifname.Data() " << ifname.Data() << endl;
  sscanf(ifname.Data(),"run%04d/%*s",&run);

  //Open the input and output files.
  TFile *ofile = new TFile(RootFile,"RECREATE");
  cout<<"input file: "<<InputFile<< endl;
  cout<<"run number " << run << endl;
  cout<<"writing to output file: "<<RootFile<< endl;
  cout << "------------------------------------" << endl;
  FILE *infile;

  string infilestr = string(InputFile);
  char extension[20];
  strcpy(extension, infilestr.substr(infilestr.find_last_of(".")+1).c_str());
  if((strcmp(extension,"gz")==0 || strcmp(extension,"gzip")==0)
     && !noHFC){
    infile = popen(Form("zcat %s | GEB_HFC -p ", InputFile),"r");
  } else if (strcmp(extension,"dat")==0 && !noHFC){
    infile = popen(Form("GEB_HFC -p %s ", InputFile),"r");
  } else if ((strcmp(extension,"gz")==0 || strcmp(extension,"gzip")==0)
	     && noHFC){
    infile = popen(Form("zcat %s ", InputFile),"r");
  } else if (strcmp(extension,"dat")==0 && noHFC){
    infile = fopen(InputFile,"r");
  } else {
    cout << "Unknown file type.  Will try to open as a .dat file, using HFC" << endl;
    infile = popen(Form("GEB_HFC -p %s ", InputFile),"r");
  }

  if(infile == NULL){
    cout << "Sorry I couldn't find the file: " << InputFile << ". Aborting ..." << endl;
    ofile->Close();
    return 3;
  }
  //Read the settings
  Settings* set = new Settings(SettingFile);
  RunInfo* info = new RunInfo();
  ofile->cd();
  int vl = set->VLevel();
  info->SetHIRunNumber(run);
  set->Write("settings",TObject::kOverwrite);
  //Initialize the data structures for the event building.
  struct crys_ips_abcd6789 inbuf_abcd6789[1];
  int buffers = 0;
  long long int bytes_read = 0;
  UnpackedEvent *evt = new UnpackedEvent(set);
#ifdef SIMULATION
  evt->SetWrite(wrawtree, wrawhist, wcaltree, wcalhist,false);
#else
  evt->SetWrite(wrawtree, wrawhist, wcaltree, wcalhist);
#endif

  Calibration *cal = new Calibration(set);
  evt->SetCalibration(cal);
  evt->SetVL(vl);
  evt->Init();
  evt->SetMakeMode2(makemode2);
  
  //Loop over the entirety of the input file.
  while(!feof(infile) && !signal_received){

    //Finish reading if you have read as many buffers as the user has requested.
    if(LastBuffer > 0 && buffers >= LastBuffer)
      break;

    size_t bsize;
    unsigned short buffer[4096], evtlength[1], *pevent;

    //Read the header of the new buffer.
    unsigned int header[4];
    bsize = fread(header, sizeof(unsigned int), 4, infile);
    long long int ts = (long long int)header[3] << 32;
    ts += header[2];
    if(vl>1){
      cout << "header" << endl;
      for(int i=0;i<4;i++)
	cout << "0x"<<(hex)<< header[i] <<"\t"<<(dec)<< header[i] <<endl;
    }
    bytes_read += 4*sizeof(unsigned int);
    //Decode the event differently depending on which type of data is identified in the header.
    //For each, pass the data from the file into the UnpackedEvent to be read.

    if(header[0]==GRETINA_ID){
      if(vl>1){
	cout << "---------------------starting gretina------------------------------- " << endl;
	cout << "gret timestamp:\t"<< ts << "\tlength: "<< header[1] << "\thex: " <<(hex) << ts << "\tlength: "<< header[1] <<(dec)<< endl;
      }
      Crystal* crys;
      if(header[1]==sizeof(crys_ips_abcd6789)){
	bsize = fread(&inbuf_abcd6789[0], sizeof(crys_ips_abcd6789), 1, infile);
	bytes_read += sizeof(struct crys_ips_abcd6789);
	crys = new Crystal(inbuf_abcd6789[0]);
      } 
      else{
	cout << "Unknown size for mode2 data" << endl;
	break;
      }
      int error = evt->DecodeGretina(crys,ts);
      if(error){
	if(error==1)
	  cout << GREEN << "end of file reached " << DEFCOLOR << endl;
	else{
	  cout << "An error ("<<error<<") occured at buffer nr " << buffers << " in DecodeGretina() while processing file: " << InputFile << ". Continuing ..." << endl;
	  continue;
	}
      }
    }
    if(header[0]==TRACE_ID){
      if(vl>1){
	cout << "---------------------starting trace------------------------------- " << endl;
	cout << "trace timestamp:\t"<< ts << "\tlength: "<< header[1] << "\thex: " <<(hex) << ts << "\tlength: "<< header[1] <<(dec)<< endl;
      }
      char cBuf[16382*4];
      if(header[1]>4*16382)
	cout << "buffer nr " << buffers << ", size of cBuf not sufficient " << endl;
      bsize = fread(cBuf, sizeof(char), header[1], infile);
      bytes_read += header[1]*sizeof(char);
      int error = evt->DecodeMode3(cBuf, header[1], ts);
      if(error){
	if(error==1)
	  cout << GREEN << "end of file reached " << DEFCOLOR << endl;
	else{
	  cout << "An error ("<<error<<") occured at buffer nr " << buffers << " in DecodeMode3() while processing file: " << InputFile << ". Continuing ..." << endl;
	  continue;
	}
      }
    }
    //Write the trees out to disk every denom events.
    buffers++;
    if(buffers % denom == 0){
      // if(buffers % (denom*1000) == 0){
      // 	if(wrawtree)
      // 	  evt->GetTree()->AutoSave();
      // 	if(wcaltree)
      // 	  evt->GetCalTree()->AutoSave();
      // }
      double time_end = get_time();
      cout << "\r" << buffers << " buffers read... "<<bytes_read/(1024*1024)<<" MB... "<<buffers/(time_end - time_start) << " buffers/s" << flush;
    }
  }
  //Finish reading the last event and close it out.
  evt->WriteLastEvent();
  cout << endl;
  cout << "Total of " << BLUE << buffers << DEFCOLOR << " data buffers ("<< BLUE <<bytes_read/(1024*1024) << DEFCOLOR << " MB) read." << endl;
  info->SetHIBytes(bytes_read);
  info->SetHIRawEvents(evt->NrOfEvents());
  if(wrawtree||wrawhist){
    cout << "Total of " << BLUE << evt->NrOfEvents() << DEFCOLOR << " raw events";
    if(wrawtree)
      cout <<" ("<< BLUE << evt->GetTree()->GetZipBytes()/(1024*1024)<< DEFCOLOR << " MB)";
    cout <<" written."  << endl;
  }
  if(wcaltree){
    info->SetHICalEvents(evt->NrOfCalEvents());
    cal = evt->GetCalibration();
    info->SetBigRIPSCtr(cal->GetBigRIPSCtr());
    info->SetBigRIPSHitCtr(cal->GetBigRIPSHitCtr());
    info->SetHiCARICtr(cal->GetHiCARICtr());
    info->SetHiCARIHitCtr(cal->GetHiCARIHitCtr());
    cout << "Total of " <<BLUE << evt->NrOfCalEvents() << DEFCOLOR << " cal events ("<< BLUE <<evt->GetCalTree()->GetZipBytes()/(1024*1024) << DEFCOLOR << " MB) written."  << endl;
  }
  cout << endl;
  ofile->cd();
  //Final cleanup and writing of files.
  if(wrawtree){
    evt->GetTree()->Write("",TObject::kOverwrite);
  }
  if(wcaltree){
    evt->GetCalTree()->Write("",TObject::kOverwrite);
  }
  info->Write("info",TObject::kOverwrite);
  double time_end = get_time();
  cout << "Program Run time " << time_end - time_start << " s." << endl;
  cout << "Unpacked " << buffers/(time_end - time_start) << " buffers/s." << endl;
  timer.Stop();
  cout << "CPU time: " << timer.CpuTime() << "\tReal time: " << timer.RealTime() << endl;
  return 0;
}
