#include "Trace.hh"
//#include "Settings.hh"
#include <algorithm>

bool trace_comparison(Trace a, Trace b){
  return a.GetEnergy() > b.GetEnergy(); //Greater-than, for highest-first sort.
}
/*
crys_ips_abcd5678 Mode3Hit::MakeMode2(Settings* set){
  crys_ips_abcd5678 output;

  //Zero everything out that isn't used.
  output.tot_e = 0;
  output.trig_time = 0;
  output.t0 = 0;
  output.cfd = 0;
  output.chisq = 0;
  output.norm_chisq = 0;
  output.baseline = 0;
  output.prestep = 0;
  output.poststep = 0;
  output.pad = 0;
  for(int i=0; i<MAX_INTPTS; i++){
    output.ips[i].x = 0;
    output.ips[i].y = 0;
    output.ips[i].z = 0;
    output.ips[i].e = 0;
    output.ips[i].seg = 0;
    output.ips[i].seg_ener = 0;
  }

  //Sort traces by descending energy
  std::sort(ftrace.begin(), ftrace.end(), trace_comparison);

  //Fill parameters
  output.type = 0xabcd5678;
  output.num = 0;
  output.timestamp = fts;
  output.crystal_id = ftrace.at(0).GetID();


  for(int i=0; i<fmult; i++){
    Trace tr = ftrace.at(i);

    int segnum = tr.GetSegNum(set); //0-31 for real segnums, -4 to -1 for CC
    if(segnum>=0){
      if(output.num<MAX_INTPTS){
	output.ips[output.num].seg = segnum;
	output.ips[output.num].seg_ener = tr.GetEnergy();
	output.num++;
      }
    } else {
      int boardnum = -segnum - 1;
      output.core_e[boardnum] = tr.GetEnergy() >> 7;
    }
  }
  output.tot_e = output.core_e[0]/2;

  return output;
}

int Trace::GetSegNum(Settings* set){
  //Central contact is same for all crystals.
  if(fchn==9){
    return -(fboard-3)-1;
  }

  int standard_val = (fboard-3)*9+fchn;

  int detnum = set->Hole2Det(fhole);
  if (detnum==1 && (fcrystal==0 || fcrystal==2)) {
    //Q2, P0 and P2 have segments 16 and 22 swapped
    if(standard_val==16){
      return 22;
    } else if (standard_val==22){
      return 16;
    } else {
      return standard_val;
    }
  } else if (detnum==0 && (fcrystal==0 || fcrystal==2)){
    //Q1, P0 and P2 are just weird
    switch(standard_val){
    case 0:
      return 1;
    case 1:
      return 8;
    case 2:
      return 3;
    case 3:
      return 4;
    case 4:
      return 35;
    case 5:
      return 0;
    case 6:
      return 7;
    case 7:
      return 14;
    case 8:
      return 9;
    case 9:
      return 10;
    case 10:
      return 29;
    case 11:
      return 30;
    case 12:
      return 13;
    case 13:
      return 20;
    case 14:
      return 15;
    case 15:
      return 16;
    case 16:
      return 23;
    case 17:
      return 24;
    case 18:
      return 19;
    case 19:
      return 26;
    case 20:
      return 21;
    case 21:
      return 22;
    case 22:
      return 17;
    case 23:
      return 18;
    case 24:
      return 25;
    case 25:
      return 32;
    case 26:
      return 27;
    case 27:
      return 28;
    case 28:
      return 11;
    case 29:
      return 12;
    case 30:
      return 31;
    case 31:
      return 2;
    case 32:
      return 33;
    case 33:
      return 34;
    case 34:
      return 5;
    case 35:
      return 6;
    default:
      return -1;
    }
  } else {
    //Everything else is normal
    return standard_val;
  }

  cout << "This message in " << __PRETTY_FUNCTION__ << " should never be printed" << endl;
}
*/
