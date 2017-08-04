/*
* @Author: sushil
* @Date:   2017-07-28 11:33:02
* @Last Modified by:   sushil
* @Last Modified time: 2017-08-02 12:46:07
* @Desc: CPU program to convert the local coordnate of
*  CPE to global coordinate for doblets
* @Input: Output of CPE which consists of hit in cm, GlobalPosition and 
*  Rotation matrix for each module
* @Output: Global coordinate of hits
*/
#include "LocalToGlobal.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;

void parseInput(std::vector<LocalHit> &v) {
  ifstream ifile("CPE_GPU.txt");
  if (!ifile) {
    cout<<"File not found: CPE_GPU.txt"<<endl;
  }
  string line;
  getline(ifile, line);
  long long unsigned hitid;
  int event, module;
  float xhit, yhit;
  LocalHit lhit;
  while(!ifile.eof()) {
    ifile>>event>>module>>hitid>>xhit>>yhit;
    lhit.module = module;
    lhit.HitId = hitid;
    lhit.x = xhit;
    lhit.y = yhit;
    v.push_back(lhit);
  }
}

RecHit toGlobal(GlobalPosition *gp, int module, float x, float y) {
  float xpos = gp[module].xpos;
  float ypos = gp[module].ypos;
  float zpos = gp[module].zpos;
  float r    = gp[module].r;
  float phi  = gp[module].phi;
  Rotation rot = gp[module].Rot;
  float R11  = rot.R11;
  float R12  = rot.R12;
  float R13  = rot.R13;
  float R21  = rot.R21;
  float R22  = rot.R22; 
  float R23  = rot.R23;
  float R31  = rot.R31;
  float R32  = rot.R32;
  float R33  = rot.R33;

  float z =0; // as there is no local z 2D module
  // local to global: Rota[]*local[] + pos[]
  float global_x = (R11*x + R21*y + R31*z) + xpos;
  float global_y = (R12*x + R22*y + R32*z) + ypos;
  float global_z = (R13*x + R23*y + R33*z) + zpos;
  
  RecHit hit;
  hit.x = global_x;
  hit.y = global_y;
  hit.z = global_z;
  // barrel: u=r, v=z, forward the opposite...
  if(module<1184) {
    hit.u = r;
    hit.v = global_z;
  }
  else {
    hit.u = global_z;
    hit.v = r;
  }
  hit.phi = phi;
  return hit;
}
void localToGlobal(GlobalPosition *gp, const vector<LocalHit> &lhitv,
  vector<RecHit> &hitv) {
  int module;
  float x,y;
  for(auto it = lhitv.begin(); it!=lhitv.end(); it++) {
    LocalHit lhit = *it;
    module = lhit.module;
    x = lhit.x;
    y = lhit.y;
    RecHit hit = toGlobal(gp, module, x,y);
    hit.HitId = lhit.HitId;
    hitv.push_back(hit);
  } 
}

void storeOutput(const std::vector<RecHit> &hitv, 
  const std::vector<LocalHit> lhitv)  {
  ofstream ofile("GlobalHit_CPU_standalone.txt");
  ofile<<"   HitId\t\t localx\t  localy\t  globalx\t   globaly \t  globalz"<<endl;
  ofile<<std::fixed;
  ofile<<setprecision(6);
  auto lit = lhitv.begin();
  for(auto git = hitv.begin(); git!=hitv.end(); git++) {
    LocalHit lhit = *lit;
    RecHit ghit = *git;
    ofile<<setw(14)<<lhit.HitId<<setw(14)<<lhit.x<<setw(14)<<lhit.y<<setw(14);
    ofile<<ghit.x<<setw(14)<<ghit.y<<setw(14)<<ghit.z<<endl;
    lit++;
  }
  ofile.close();
}
int main(){
  GlobalPosition *globalPosRot;
  globalPosRot = (GlobalPosition*)malloc(NMODULE*sizeof(GlobalPosition));
  uploadGlobal_Positon_Rotation_Matrix(globalPosRot);
  std::vector<LocalHit> lhit; //input
  std::vector<RecHit> hit;//output
  parseInput(lhit);
  localToGlobal(globalPosRot,lhit, hit);
  storeOutput(hit, lhit);
  return 0;
}

// upload the global position and rotation matrix for 
// local to global coordinate coversion
void uploadGlobal_Positon_Rotation_Matrix(GlobalPosition *gp) {
  // read the file and upload
  ifstream ifile("Global_Position_Rotation_forL2G.dat");
  if(!ifile) {
    cout<<"File not found: Global_Position_Rotation_forL2G.dat"<<endl;
  }
  string line;
  getline(ifile, line);
  for(int i=0;i<NMODULE;i++) {
    ifile>>gp[i].RawId>>gp[i].xpos>>gp[i].ypos>>gp[i].zpos>>gp[i].r>>gp[i].phi;
    ifile>>gp[i].Rot.R11>>gp[i].Rot.R12>>gp[i].Rot.R13;
    ifile>>gp[i].Rot.R21>>gp[i].Rot.R22>>gp[i].Rot.R23;
    ifile>>gp[i].Rot.R31>>gp[i].Rot.R32>>gp[i].Rot.R33;
  }
}
