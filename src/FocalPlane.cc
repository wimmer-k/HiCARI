#include "FocalPlane.hh"
int fpID[6] = {3,5,7,8,9,11};
int fpNr(int id){
  if(id<3 || id>11)
    return -1;
  switch (id){
    case 3:
      return 0;
    case 5:
      return 1;
    case 7:
      return 2;
    case 8:
      return 3;
    case 9:
      return 4;
    case 11:
      return 5;
    default:
      return -1;
  }
}
int firstPPAC(int id){
  if(id<3 || id>11)
    return -1;
  switch (id){
    case 3:
      return 5;
    case 5:
      return 10;
    case 7:
      return 15;
    case 11:
      return 33;
    default:
      return -1;
  }
}
