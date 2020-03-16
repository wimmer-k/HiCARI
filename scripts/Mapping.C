//#include <bits/stdc++.h> 
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
using namespace std;
void Mapping(){
  // three dimensional key 
  map<int, map<int, map<int, int> > > m; 
  
  // for accessing first map 
  map<int, map<int, map<int, int> > >::iterator ftr; 
  
  // for accessing second map 
  map<int, map<int, int> >::iterator str; 
  
  // for accessing inner map 
  map<int, int>::iterator ptr; 
  m[11][0][0] = 01;
  m[11][0][1] = 02;
  m[11][0][2] = 03;
  m[11][0][3] = 11;
  m[11][1][0] = 12;
  m[11][1][1] = 13;
  m[11][1][2] = 21;
  m[11][1][3] = 22;
  m[11][2][0] = 23;
  m[11][2][1] = 31;
  m[11][2][2] = 32;
  m[11][2][3] = 33;
  m[11][3][0] = 41;
  m[11][3][1] = 42;
  m[11][3][2] = 43;
  m[11][3][3] = 51;
  
  m.insert(make_pair(12, map<int, map<int, int> >())); 
  m[12].insert(make_pair(0, map<int, int>())); 
  m[12][0].insert(make_pair(0,52));

  char* abc = "ABC";
  
  cout << m[11][3][1] << endl;
  cout << m[11][3][4] << endl;
  cout << m[12][3][1] << endl;

  cout << "MB mapping:" << endl;
  
  for(ftr = m.begin(); ftr != m.end(); ftr++) {   
    for(str = ftr->second.begin(); str != ftr->second.end(); str++) { 
      for(ptr = str->second.begin(); ptr != str->second.end(); ptr++) { 
	cout << "Hole " << ftr->first 
	     << ", Crystal " << str->first 
	     << ", Slot " << ptr->first 
	     << " is MB " << ptr->second/10 << abc[ptr->second%10] << endl; 
      } 
    } 
  }
}


/*
void Mapping(){
  map<vector<int>, int > hcs2mc; //hole/crystal/slot -> module/crystal
  vector<int> MB0A{ 11, 0, 0}; 
  vector<int> MB0B{ 11, 0, 1}; 
  vector<int> MB0C{ 11, 0, 2}; 
  vector<int> MB1A{ 11, 0, 3}; 
  vector<int> MB1B{ 11, 1, 0}; 
  vector<int> MB1C{ 11, 1, 1}; 
  vector<int> MB2A{ 11, 1, 2}; 
  vector<int> MB2B{ 11, 1, 3}; 
  vector<int> MB2C{ 11, 2, 0}; 
  vector<int> MB3A{ 11, 2, 1}; 
  vector<int> MB3B{ 11, 2, 2}; 
  vector<int> MB3C{ 11, 2, 3}; 
  vector<int> MB4A{ 11, 3, 0}; 
  vector<int> MB4B{ 11, 3, 1}; 
  vector<int> MB4C{ 11, 3, 2}; 
  vector<int> MB5A{ 11, 3, 3};

  // vector<int> mc_MB0A{0,0};
  // vector<int> mc_MB0B{0,1};
  int mc_MB0A = 00;
  int mc_MB0B = 01;

  hcs2mc[MB0A] = mc_MB0A;

  vector<int> test{11,0,0};
  hcs2mc.find(test);
}
*/
