/*Copyright 2017 Sushil Dubey, Shashi Dugad, TIFR

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

// Cluster paramter Estimation
// Input : ClusterId[], xx[], yy[]
// Output: ClusterId  xhit(in cm)  yhit(in cm)

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <cuda.h>
#include <cuda_runtime.h>

#include "CudaError.h"
#include "PixelClusterUtil.h"
#include "EventInfoGPU.h"
//CPE specific
#include "CPEGPU.h"
// To intitalize memory
#include "CPEGPUMem.h"
using namespace std;

__host__ __device__ uint getModule(uint64 clusterId) {
  //uint event = ((clusterId  >> EVENT_shift) & EVENT_mask);
  uint module = ((clusterId >> MODULE_shift) & MODULE_mask);
  //uint xcor  = ((clusterId  >> XCOR_shift) & XCOR_mask);
  //uint ycor  = ((clusterId  >> YCOR_shift) & YCOR_mask);
  return module;
}
__host__ __device__ uint getEvent(uint64 clusterId) {
  uint event = ((clusterId  >> EVENT_shift) & EVENT_mask);
  return event;
}

// CPE kernel for a given cluster, it finds the xhit and yhit
// given in CMSSW CPE package(SiPixelRecHits)
// Input: clusterId, Index, xx, yy, adc
// output: xhit, yhit
__global__ void CPE_kernel(const CPE_cut_Param cpe_cut, const DetDB *detDB, 
                           const uint64 *ClusterId,const uint *Index,const uint *xx,const uint *yy,
                           const uint *adc, float *xhit, float *yhit ) 
{
   
  __shared__ float xmin, xmax;
  __shared__ float ymin, ymax;
  __shared__ float Q_l_X, Q_f_X, Q_l_Y, Q_f_Y;
  __shared__ float sizeX, sizeY;
  __shared__ LocalPoint lp_min, lp_max;
  __shared__ LorentzAngle cotAngle;
  __shared__ float lorentShiftX , lorentShiftY, shiftX, shiftY;
  float theThickness;

  uint blockId = blockIdx.x;
  //attention: index of last block 
  uint tid = threadIdx.x;
  uint64 clusterId = ClusterId[blockId];
  if(clusterId==0) {
    xhit[blockId] = 0;
    yhit[blockId] = 0;
    return ;
  }
  uint startIndex, size, moduleId;
  startIndex = Index[blockId];

  size  = Index[blockId+1] - startIndex;
  moduleId = getModule(clusterId); // to get the moduleId devide by 10^6
  theThickness = (moduleId<1184) ? thicknessBarrel: thicknessForward;

  // Jobs in kernel are 
  // Compute lorentzAngle (independent method) -------by 1st thread
  if(tid==0) {
    cotAngle = computeLorentzAngle(detDB,moduleId,startIndex,size, xx, yy,adc); 
    lorentShiftX = detDB->LorentzShiftX[moduleId];  // read from the database
    lorentShiftY = detDB->LorentzShiftY[moduleId];
    
    shiftX = 0.5f*lorentShiftX ;
    shiftY = 0.5f*lorentShiftY;
  
    lorentShiftY = lorentShiftY * widthLAFractionY;
    lorentShiftX = (moduleId<1184) ? lorentShiftX*widthLAFractionX_Barrel: lorentShiftX*widthLAFractionX_Forward;
  }

  // Find xmin, ymin, xmax, ymax (independent method) -------- by 2nd thread
  if(tid==1) {  
    min_max(startIndex, size, xx, xmin, xmax);
    min_max(startIndex, size, yy, ymin, ymax);
    sizeX = xmax - xmin + 1.0f;
    sizeY = ymax - ymin + 1.0f;
    //printf("xmin: %f, xmax: %f,  ymin: %f, ymax: %f\n",xmin, xmax, ymin, ymax); 
  }
  __syncthreads();

  // Find Q_f and Q_l which depend upon output of step 2.----- by 1st thread
  if(tid==0) {
    collectCharge (xx, yy, adc, startIndex, size, xmin, xmax,
                 ymin, ymax, Q_l_X, Q_f_X, Q_l_Y, Q_f_Y );
  }

  // Convert to localPosition in cm depends upon output of steps 2.-----by 2ns thread
  if(tid==1) {
    lp_min = localPositionInCm( xmin +1.0, ymin+1.0); // use the formula to convert
    lp_max = localPositionInCm( xmax, ymax);          // first pix and last pix
  }
  __syncthreads();

  // Compute x_hit using the formula depends upon output of step 1 to 4 ----- by 1st thread
  if(tid==0) {
    float x_hit=genericPixelHit(sizeX, lp_min.x(), lp_max.x(),
                  Q_f_X, Q_l_X,
                  cotAngle.cotAlpha,
                  pitchX,
                  theThickness,
                  lorentShiftX,
                  isItBigPixelInX((int)xmin),
                  isItBigPixelInX((int)xmax),
                  cpe_cut.the_eff_charge_cut_lowX,
                  cpe_cut.the_eff_charge_cut_highX,
                  cpe_cut.size_cutX
                  );
    x_hit = x_hit + shiftX;
    xhit[blockId] = x_hit;
  }
  // Compute y_hit using the formula depends upon output of step 1 to 4----- by 2nd thread
  if(tid==1) {
    float y_hit=genericPixelHit(sizeY, lp_min.y(), lp_max.y(),
                  Q_f_Y, Q_l_Y,
                  cotAngle.cotBeta,
                  pitchY,
                  theThickness,
                  lorentShiftY,
                  isItBigPixelInY((int)ymin),
                  isItBigPixelInY((int)ymax),
                  cpe_cut.the_eff_charge_cut_lowY,
                  cpe_cut.the_eff_charge_cut_highY,
                  cpe_cut.size_cutY
                  );
    y_hit = y_hit + shiftY;
    yhit[blockId] = y_hit;
  }
}
// end of CPE_kernel

// device function to calculate the actual pixel hit
// the function is taken from the CMSSW
__device__ float genericPixelHit(uint size, float first_pix, float last_pix,
                      float Q_f, float Q_l, float cot_angle, float pitch,
                      float theThickness, float lorentz_shift,
                      bool first_is_big, bool last_is_big,
                      float eff_charge_cut_low, 
                      float eff_charge_cut_high,
                      float size_cut)

{ 
  float geom_center = 0.5f*(first_pix + last_pix);
  // The case of only one pixel in this projection is separate.  Note that
  // here first_pix == last_pix, so the average of the two is still the
  // center of the pixel.
  if ( size == 1 ) {return geom_center;}

  // Width of the clusters minus the edge (first and last) pixels.
  // In the note, they are denoted x_F and x_L (and y_F and y_L)
  float W_inner = last_pix - first_pix;  // in cm

  // Predicted charge width from geometry
  float W_pred = theThickness * cot_angle - lorentz_shift;// geometric correction (in cm)

  //--- Total length of the two edge pixels (first+last)
  float sum_of_edge = 2.0f;

  if (first_is_big) sum_of_edge += 1.0f;
  if (last_is_big)  sum_of_edge += 1.0f;
  
  //--- The `effective' charge width -- particle's path in first and last pixels only
  if(W_pred<0) W_pred=0-W_pred;
  float W_eff = W_pred - W_inner;

  //--- If the observed charge width is inconsistent with the expectations
  //--- based on the track, do *not* use W_pred-W_innner.  Instead, replace
  //--- it with an *average* effective charge width, which is the average
  //--- length of the edge pixels.
  //
  //  bool usedEdgeAlgo = false;
  if ( (size >= size_cut) || (
       ( W_eff/pitch < eff_charge_cut_low ) |
       ( W_eff/pitch > eff_charge_cut_high ) ) ) {
      W_eff = pitch * 0.5f * sum_of_edge;  // ave. length of edge pixels (first+last) (cm)
  }
    
  //--- Finally, compute the position in this projection
  float Qdiff = Q_l - Q_f;
  float Qsum  = Q_l + Q_f;

  //--- Temporary fix for clusters with both first and last pixel with charge = 0
  if(Qsum==0) Qsum=1.0f;

  float hit_pos = geom_center + 0.5f*(Qdiff/Qsum) * W_eff;
  return hit_pos;
}


void CPE_wrapper(const uint total_cluster, const uint64 *ClusterId, const uint *Index, const uint *xx, const uint *yy,
                 const uint *adc ) 
{  
  CPE_cut_Param cpe_cut;
  int no_blocks = total_cluster;
  int no_threads = 2;
  // xhit_d, yhit_d, contains output
  CPE_kernel<<<no_blocks, no_threads>>>(cpe_cut,detDB,ClusterId, Index, xx, yy, adc, xhit_d, yhit_d); 
  cudaDeviceSynchronize();
  checkCUDAError("Error in CPE_kernel");
  cout<<"CPE execution finished!\n";
}

// compute cot alpha and beta for each cluster
// formula to calculate cot alpha and cot beta is taken from
// http://cmslxr.fnal.gov/source/RecoLocalTracker/SiPixelRecHits/src/PixelCPEBase.cc?v=CMSSW_8_1_0#0325
// https://cmssdt.cern.ch/lxr/source/DataFormats/SiPixelCluster/interface/SiPixelCluster.h#0104
// Input: as shown in fucntion
// Output: LorentzAngle cotAlpha and cotBeta
__device__ LorentzAngle computeLorentzAngle(const DetDB *detDB, const uint moduleId,
                        const uint startIndex, const uint size, const uint *xx,
                        const uint *yy, const uint *adc ) 
{

  float totalCharge = 0.0f;
  float xc=0.0f, yc=0.0f;
  uint i= startIndex;
  uint end = startIndex+size;
  for (; i<end; i++) {
    xc += (xx[i] + 0.5f)*adc[i];
    yc += (yy[i] + 0.5f)*adc[i];
    totalCharge += adc[i];
  }
  xc = xc/totalCharge;
  yc = yc/totalCharge;
  LocalPoint lp = localPositionInCm(xc, yc);

  float gvx = lp.x() - detDB->X0[moduleId];
  float gvy = lp.y() - detDB->Y0[moduleId];
  float gvz = -1.f/detDB->Z0[moduleId];
  float cot_aplha = gvx*gvz;
  float cot_beta  = gvy*gvz;
  LorentzAngle la = { cot_aplha, cot_beta }; 
  return la;
}

// device function to find min max
__device__ void min_max(uint startIndex,uint size,
                        const uint *xx, float &xmin,float &xmax) 
{
  
  xmin = xx[startIndex];
  xmax = xx[startIndex];
  uint i = startIndex+1;
  for(; i<(startIndex+ size); i++) {
    if(xmin>float(xx[i])) xmin = xx[i];
    if(xmax<float(xx[i])) xmax = xx[i];
  }
}

// device function to collect charges on
// the edge of the cluster
__device__ void 
collectCharge (const uint *xx, const uint *yy, const uint *adc,
               uint startIndex, uint size, float xmin, float xmax,
               float ymin, float ymax, float &Q_l_X, float &Q_f_X,
               float &Q_l_Y, float &Q_f_Y ) 
{
  Q_f_X = 0.0f;
  Q_l_X = 0.0f;
  Q_f_Y = 0.0f;
  Q_l_Y = 0.0f;
  float pix_adc = 0.0f;
  uint i=startIndex;
  for(; i<(startIndex+size); i++) {
    // upper cut is put on pixel charge but does not affect the result much 
    pix_adc = adc[i];//pixel.ADC > 13200.0f ? 13200.0f:pixel.ADC;
    if((float(xx[i])==xmin)) Q_f_X+= pix_adc;
    if((float(xx[i])==xmax)) Q_l_X+= pix_adc;
    if((float(yy[i])==ymin)) Q_f_Y+= pix_adc;
    if((float(yy[i])==ymax)) Q_l_Y+= pix_adc;
  }
} 

// this function converts pixel coordinates row and col in cm
// multiply the row and col by pitch size to convert in cm
// Input: x(0-159), y(0-145)
// Output: x(-0.81 cm to +0.81 cm), y(-3.24 cm to +3.24 cm) 
__device__ LocalPoint localPositionInCm(float x, float y) {
  //  m_xoffset = -(m_nrows + BIG_PIX_PER_ROC_X*m_nrows/ROWS_PER_ROC)/2. * 
  //  m_pitchx;
  //  m_yoffset = -(m_ncols + BIG_PIX_PER_ROC_Y*m_ncols/COLS_PER_ROC)/2. * 
  //  m_pitchy;
  // m_nrows = 160, BIG_PIX_PER_ROC_X=1, ROWS_PER_ROC=80,m_pitchx=0.01
  // m_ncols = 416, BIG_PIX_PER_ROC_Y=2, COLS_PER_ROC = 52, m_pitchy = 0.015
  // after calculating  
  float m_xoffset = 0.81f;
  float m_yoffset = 3.24f;
  // As big pixel issue is corrected in CMSSW_9_2_0
  int binoffx = int( x );        // truncate to int
  float fractionX = x - float(binoffx); // find the fraction 
  float local_pitchx = pitchX;   // default pitch
   
  if (binoffx>80) {            // ROC 1 - handles x on edge cluster
    binoffx=binoffx+2;
  } 
  else if (binoffx==80) {    // ROC 1
    binoffx=binoffx+1;
    local_pitchx *= 2;
  }
  else if (binoffx==79) {      // ROC 0
    binoffx=binoffx+0;
    local_pitchx *= 2;    
  } 
  // The final position in local coordinates 
  //float lpX = float( binoffx * m_pitchx ) + fractionX * local_pitchx + m_xoffset;
  float xcm = float(binoffx * pitchX) + fractionX * local_pitchx - m_xoffset;

  int binoffy = int( y );        // truncate to int
  float fractionY = y - float(binoffy); // find the fraction 
  float local_pitchy = pitchY;   // defaultpitch
  // 415 is last big pixel, 416 and above do not exists!
  //constexpr int bigYIndeces[]{0,51,52,103,104,155,156,207,208,259,260,311,312,363,364,415};
  //auto const j = std::lower_bound(std::begin(bigYIndeces),std::end(bigYIndeces),binoffy);
  //if (*j==binoffy) { local_pitchy  *= 2 ;}
  //binoffy += (j-bigYIndeces);
  if(binoffy>416) binoffy=432; //this is due to the bug in cmssw cpe, since origin is shifted by 1 remove this in cmssw
  else if(!(binoffy%52)) {
    binoffy += ((int)(binoffy/52))*2;
    local_pitchy  *= 2 ;
  }
  else {
    binoffy += ((int)(binoffy/52))*2 +1;
    if(!(binoffy+1)%52) local_pitchy  *= 2 ;
  }
  // The final position in local coordinates 
  float ycm = float(binoffy*pitchY) + fractionY*local_pitchy - m_yoffset;

  //float xcm = x*pitchX - m_xoffset;
  //float ycm = y*pitchY - m_yoffset;
  LocalPoint lp;
  lp.xcm=xcm;
  lp.ycm=ycm;
  return lp;
}

//-------------------------------------------------------------
// Return the BIG pixel information for a given pixel
// reference: http://cmslxr.fnal.gov/source/Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h?v=CMSSW_9_2_0#0119
__device__ bool isItBigPixelInX( const int ixbin ) {
  return (( ixbin == 79 ) || ( ixbin == 80 ));
}

__device__ bool isItBigPixelInY( const int iybin ) {
  int iybin0 = iybin%52;
  return(( iybin0 == 0 ) || ( iybin0 == 51 ));
}

void initDeviceMemCPE() {
  const int MAX_CLUSTER = 30000*NEVENT;
  cudaMalloc((void**)&xhit_d, MAX_CLUSTER*sizeof(float));
  cudaMalloc((void**)&yhit_d, MAX_CLUSTER*sizeof(float));
  cudaMallocManaged((void**)&detDB, sizeof(DetDB));
  uploadCPE_db(detDB);
}
void freeDeviceMemCPE() {
  cudaFree(xhit_d);
  cudaFree(yhit_d);
  cudaFree(detDB);
}

// upload the CPE database to the GPU memory
// they are constant for the module
void uploadCPE_db(DetDB *detDB) {
  uint moduleId,rawId,i=0;
  float X0, Y0, Z0, Rdet, Zdet, LShiftX, LShiftY;
  ifstream ifile("Pixel_CPE_Phase1_database_C_920.dat");
  if(!ifile) {
    cout<<"CPE database file does not exit !"<<endl;
    cout<<"File: Pixel_CPE_Phase1_database_C_920.dat"<<endl;
    exit(1);
  }
  string str;
  getline(ifile, str);
  while(!ifile.eof()) {
    ifile>>moduleId>>rawId>>X0>>Y0>>Z0>>Rdet>>Zdet>>LShiftX>>LShiftY;
    detDB->RawId[i] = rawId;
    detDB->X0[i] = X0;
    detDB->Y0[i] = Y0;
    detDB->Z0[i] = Z0;
    detDB->LorentzShiftX[i] = LShiftX;
    detDB->LorentzShiftY[i] = LShiftY;
    i++ ;
  }
  ifile.close();
  cout<<"CPE database uploaded successfully ! "<<endl;
}