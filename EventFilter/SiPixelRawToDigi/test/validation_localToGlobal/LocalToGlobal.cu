/*
* @Author: sushil
* @Date:   2017-07-28 11:33:02
* @Last Modified by:   sushil
* @Last Modified time: 2017-07-28 18:59:53
* @Desc: standalone GPU program to convert the local coordnate of
*  CPE to global coordinate for doblets
* @Input: Output of CPE which consists of hit in cm, GlobalPosition and 
*  Rotation matrix for each module
* @Output: Global coordinate of hits
*/
#include "LocalToGlobal.cuh"

#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cuda.h>
#include <cuda_runtime.h>
#include "PixelClusterUtil.h"

using namespace std;
// to be removed from cmssw
typedef long long unsigned uint64;
float *xhit_h, *yhit_h, *xhit_d, *yhit_d;
uint64 *hitId_h, *hitId_d;



int parseInput(float *xhit_h, float *yhit_h, uint64 *hitId_h, float *xhit_d,
  float *yhit_d,uint64 *hitId_d) {
  ifstream ifile("CPE_GPU.txt");
  if (!ifile) {
    cout<<"File not found: CPE_GPU.txt"<<endl;
  }
  string line;
  getline(ifile, line);
  int event, module;
  int i=0;
  uint64 hitid;
  float xhit, yhit;
  while(!ifile.eof()) {
    ifile>>event>>module>>hitid>>xhit>>yhit;
    hitId_h[i] = hitid;
    xhit_h[i] =xhit;
    yhit_h[i] =yhit;
    i++;
  }
  cudaMemcpy(xhit_d, xhit_h, i*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(yhit_d, yhit_h, i*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(hitId_d, hitId_h, i*sizeof(uint64), cudaMemcpyHostToDevice);
  return i;
}


void initDeviceMemory() {
  // allocate memory to hold the global hits and other parameter 
  const int NEVENT = 8;
  const int MAXCLUSTER = 20000;
  const int size = MAXCLUSTER*NEVENT*sizeof(RecHit);
  cudaMalloc((void**)&Hit, size);
  cudaMalloc((void**)&globalPosRot, NMODULE*sizeof(GlobalPosition));
  uploadGlobal_Positon_Rotation_Matrix(globalPosRot);
  //temporary to be removed from cmssw
  xhit_h = (float*)malloc(MAXCLUSTER*NEVENT*sizeof(float));
  yhit_h = (float*)malloc(MAXCLUSTER*NEVENT*sizeof(float));
  hitId_h  = (uint64*)malloc(MAXCLUSTER*NEVENT*sizeof(uint64));
  cudaMalloc((void**)&xhit_d, MAXCLUSTER*NEVENT*sizeof(float));
  cudaMalloc((void**)&yhit_d, MAXCLUSTER*NEVENT*sizeof(float));
  cudaMalloc((void**)&hitId_d, MAXCLUSTER*NEVENT*sizeof(uint64));
}
void freeDeviceMemory() {
  cudaFree(Hit);
  cudaFree(xhit_d);
  cudaFree(yhit_d);
  cudaFree(hitId_d);
  cudaFree(globalPosRot);
}

__host__ __device__ int getModuleId(uint64 clusterId) {
  //uint event = ((clusterId  >> EVENT_shift) & EVENT_mask);
  int module = ((clusterId >> MODULE_shift) & MODULE_mask);
  //uint xcor  = ((clusterId  >> XCOR_shift) & XCOR_mask);
  //uint ycor  = ((clusterId  >> YCOR_shift) & YCOR_mask);
  return module;
}

__device__ RecHit toGlobal(const GlobalPosition *gp, const int module,
  const float x,const float y) {
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

__global__ void localToGlobal_kernel(const int N, const GlobalPosition *globalPosRot,
  const float *lxhit, const float *lyhit, uint64 *hitId, RecHit *Hit) {
  int gIndex = threadIdx.x + blockIdx.x*blockDim.x;
  if(gIndex<N) {
    int module = getModuleId(hitId[gIndex]); // correct the first entry clusterId =0 bad hit
    RecHit hit = toGlobal(globalPosRot, module, lxhit[gIndex], lyhit[gIndex]);
    Hit[gIndex].HitId = hitId[gIndex];
    Hit[gIndex].x = hit.x;
    Hit[gIndex].y = hit.y;
    Hit[gIndex].z = hit.z;
    Hit[gIndex].u = hit.u;
    Hit[gIndex].v = hit.v;
    Hit[gIndex].phi = hit.phi;
  }
}

void storeOutput(const int N, float *lxhit, float *lyhit, const RecHit *Hit_d) {
  RecHit *Hit_h = (RecHit*)malloc(N*sizeof(RecHit));
  cudaMemcpy(Hit_h, Hit_d, N*sizeof(RecHit), cudaMemcpyDeviceToHost);

  ofstream ofile("GlobalHit_GPU_standalone.txt");
  ofile<<"   HitId\t\t localx\t  localy\t  globalx\t   globaly \t  globalz"<<endl;
  ofile<<std::fixed;
  ofile<<setprecision(6);
  for(int i=0;i<N;i++) {
    ofile<<setw(14)<<Hit_h[i].HitId<<setw(14)<<lxhit[i]<<setw(14)<<lyhit[i]<<setw(14);
    ofile<<Hit_h[i].x<<setw(14)<<Hit_h[i].y<<setw(14)<<Hit_h[i].z<<endl;
  }
  ofile.close();
  free(Hit_h);
}

int main(){

  initDeviceMemory();
  int N = parseInput(xhit_h, yhit_h, hitId_h, xhit_d, yhit_d, hitId_d);
  int threads = 512;
  int blocks  = N/threads +1; 
  localToGlobal_kernel<<<blocks, threads>>>(N, globalPosRot, xhit_d, yhit_d, hitId_d, Hit);
  storeOutput(N, xhit_h, yhit_h, Hit);
  freeDeviceMemory();
  cout<<"sizeof float: "<<sizeof(float)<<endl;
  return 0;
}

// upload the global position and rotation matrix for 
// local to global coordinate coversion
void uploadGlobal_Positon_Rotation_Matrix(GlobalPosition *globalPosRot) {
  GlobalPosition *gp;
  gp = (GlobalPosition*)malloc(NMODULE*sizeof(GlobalPosition));
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
  cudaMemcpy(globalPosRot, gp, NMODULE*sizeof(GlobalPosition), cudaMemcpyHostToDevice);
  free(gp);
}
