/****** Filename: CPECP.cpp, cpe implementatio on cpu ******
 *  This program will calculate the pixel hit 
 *  i.e perform the cluster parameter estimation.
 *  Given cluster, it will find the xhit and yhit in 
 *  cm.
 *  It uses the input file extracted from the cmssw 
 *  package RecolTracker/SiRecHits/src/PixelGeneric.cc
 *  ps
 *  The output of this program will be validated 
 *  against the cmmsw output of the same.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <stdio.h>
#include "PixelClusterUtil.h"
#include "CPECPU.h"
typedef unsigned int uint;
typedef unsigned long long uint64;

using namespace std;

// to extract the moduleId from the clusterId
uint getModule(uint64 clusterId) {
  //uint event = ((clusterId  >> EVENT_shift) & EVENT_mask);
  uint module = ((clusterId >> MODULE_shift) & MODULE_mask);
  //uint xcor  = ((clusterId  >> XCOR_shift) & XCOR_mask);
  //uint ycor  = ((clusterId  >> YCOR_shift) & YCOR_mask);
  return module;
}
uint getEvent(uint64 clusterId) {
  uint event = ((clusterId  >> EVENT_shift) & EVENT_mask);
  return event;
}

LocalPoint pixelHit( const float (&detDB)[NMODULE][NTUPLE],const int moduleId, const std::vector<Pixel> &cluster) {
  // All these parameters to be moved in python conf
  // -----------------
  int size_cutX = 3;
  int size_cutY = 3; // to be moved into python conf
  float the_eff_charge_cut_lowX = 0.0 ;  // default
  float the_eff_charge_cut_lowY  = 0.0;
  float the_eff_charge_cut_highX = 1.0f ;// default 
  float the_eff_charge_cut_highY = 1.0f;
  // ----------------
  // temp variables
  float Q_l_X = 0.0f, Q_f_X=0.0f, Q_f_Y=0.0f, Q_l_Y =0.0f;
  float x_hit = 0.0f, y_hit=0.0f; // in cm
  float sizeX = 0.0f, sizeY =0.0f;
  float xmin, ymin, xmax, ymax;
  LorentzAngle cotAngle;
  float lorentShiftX = 0.0f, lorentShiftY = 0.0f;
  float theThickness = (moduleId<1184) ? thicknessBarrel: thicknessForward;
  ClusterParam clustParam;
  clustParam = getClusterParam(cluster);

  xmin = clustParam.xmin;
  ymin = clustParam.ymin;
  xmax = clustParam.xmax;
  ymax = clustParam.ymax;

  sizeX = xmax - xmin + 1.0f;
  sizeY = ymax - ymin + 1.0f; 
  
  collectCharge(cluster, clustParam, Q_l_X, Q_f_X, Q_l_Y, Q_f_Y);

  cotAngle =  comcputeLorentzAngle(moduleId, cluster,detDB ); // compute cot angle for each cluster
  
  lorentShiftX = detDB[moduleId][LorentzShiftX];//computeLorentShiftX(moduleId, detDB);  // read from the database
  lorentShiftY = detDB[moduleId][LorentzShiftY];//computeLorentShiftY(moduleId, detDB);
  //printf("cotAngle: %f  lorentShiftX: %f   lorentShiftY:  %f \n",cotAngle.cotAlpha,lorentShiftX, lorentShiftY );
  //printf("xmin: %f, xmax: %f,  ymin: %f, ymax: %f\n",xmin, xmax, ymin, ymax);
  //printf("Q_l_X: %f, Q_f_X: %f, Q_l_Y: %f, Q_f_Y: %f\n",Q_l_X, Q_f_X, Q_l_Y, Q_f_Y);
  float shiftX = 0.5f*lorentShiftX ;
  float shiftY = 0.5f*lorentShiftY;
  
  lorentShiftY = lorentShiftY * widthLAFractionY;
  lorentShiftX = (moduleId<1184) ? lorentShiftX*widthLAFractionX_Barrel: lorentShiftX*widthLAFractionX_Forward;
  
  LocalPoint lp_min = localPositionInCm( xmin +1.0, ymin+1.0); // use the formula to convert
  LocalPoint lp_max = localPositionInCm( xmax, ymax); // first pix and last pix

  x_hit=genericPixelHit(sizeX, lp_min.x(), lp_max.x(),
                  Q_f_X, Q_l_X,
                  cotAngle.cotAlpha,
                  pitchX,
                  theThickness,
                  lorentShiftX,
                  isItBigPixelInX(int(xmin)),
                  isItBigPixelInX(int(xmax)),
                  the_eff_charge_cut_lowX,
                  the_eff_charge_cut_highX,
                  size_cutX
                  );
  x_hit = x_hit + shiftX;

 // cout<<"lp_min().y: "<<lp_min.y()<<"  lp_max.y:"<<lp_max.y()<<endl;
  y_hit=genericPixelHit(sizeY, lp_min.y(), lp_max.y(),
                  Q_f_Y, Q_l_Y,
                  cotAngle.cotBeta,
                  pitchY,
                  theThickness,
                  lorentShiftY,
                  isItBigPixelInY(int(ymin)),
                  isItBigPixelInY(int(ymax)),
                  the_eff_charge_cut_lowY,
                  the_eff_charge_cut_highY,
                  size_cutY
                  );
  y_hit = y_hit + shiftY;
  LocalPoint lp;
  lp.xcm = x_hit;
  lp.ycm = y_hit;
  return lp;
}

float genericPixelHit(int size, float first_pix, float last_pix, float Q_f, float Q_l,
                      float cot_angle, float pitch, float theThickness, float lorentz_shift,
                      bool first_is_big, bool last_is_big,
                      float eff_charge_cut_low, 
                      float eff_charge_cut_high,
                      float size_cut) 
{ //charge_cut_high_x, charge_cut_low_x to be included
  float geom_center = 0.5f*(first_pix + last_pix);
  #ifdef DEBUG_XHIT
  cout<<"geom_center: "<<geom_center<<endl;
  #endif
  //cout<<"geom_center: "<<geom_center;
  // The case of only one pixel in this projection is separate.  Note that
  // here first_pix == last_pix, so the average of the two is still the
  // center of the pixel.
  if ( size == 1 ) {return geom_center;}

  // Width of the clusters minus the edge (first and last) pixels.
  // In the note, they are denoted x_F and x_L (and y_F and y_L)
  float W_inner = last_pix - first_pix;  // in cm
  //cout<<"   W_inner: "<<W_inner;
  // Predicted charge width from geometry
  float W_pred = theThickness * cot_angle - lorentz_shift;// geometric correction (in cm)
  //cout<<"  cot_angle: "<<cot_angle<<"   lorentShift:  "<<lorentz_shift;
  #ifdef DEBUG_XHIT
  cout<<"theThickness: "<<theThickness<<"  cot_angle: "<<cot_angle<<" lrentshift: "<<lorentz_shift<<"  size: "<<size<<endl;
  #endif 
  //--- Total length of the two edge pixels (first+last)
  float sum_of_edge = 2.0f;

  if (first_is_big) sum_of_edge += 1.0f;
  if (last_is_big)  sum_of_edge += 1.0f;
  
  //--- The `effective' charge width -- particle's path in first and last pixels only
  float W_eff = std::abs( W_pred ) - W_inner;
  //cout<<"W_inn: "<<W_inner<<" W_pre: "<<W_pred<<"  W_eff: "<<W_eff<<endl;
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
  //cout<<"W_eff: "<<W_eff<<" size: "<<size<<" size_cut: "<<size_cut<<endl;
  //cout<<"Q_f: "<<Q_f<<"  Q_l: "<<Q_l<<endl;


  //--- Temporary fix for clusters with both first and last pixel with charge = 0
  if(Qsum==0) Qsum=1.0f;

  float hit_pos = geom_center + 0.5f*(Qdiff/Qsum) * W_eff;
  //cout<<" hit_pos: "<<setprecision(9)<<hit_pos<<endl;
  //cout<<"\n";
  #ifdef DEBUG_XHIT
  cout<<"geom_center: "<<geom_center<<"  Q_diff: "<< Qdiff<<"  Qsum: "<<Qsum<<"  W_eff: "<<W_eff<<endl;
  #endif
  return hit_pos;
}

//new main it reads gpu friendly input file
// useful to compare reslt with GPU CPE
int main() {

  // upload the lorentzshift and x0 and y0 and R_det for all
  // module which is constant 
  uploadDetDB(detDB);

  pair<uint64, vector<Pixel> > mod_cl;
  vector<pair<uint, vector<Pixel> > > ClusterInfo;
  std::vector<uint64> ClusterId;
  std::vector<uint> Index; 
  ifstream input("CPE_Input_CPU_PartA.txt");
  string str;
  getline(input, str);
  uint index_temp ,i=0;
  uint64 clustId=0;
  
  while(!input.eof()) {
    input>>index_temp>>clustId;
    Index.push_back(index_temp);
    ClusterId.push_back(clustId);
    i++;
  }
  int total_cluster = i-1;
  i=0;
  cout<<"Total clusters: "<<total_cluster<<endl;
  //Index[total_cluster]= 78782;
  input.close();
  
  uint xx, yy, adc;
  uint64 pre_clustId = ClusterId[0];
  input.open("CPE_Input_CPU_PartB.txt");
  getline(input, str);
  std::vector<Pixel> v;
  Pixel p;
  while(!input.eof()) {
    input>>index_temp>>clustId>>xx>>yy>>adc;
    Pixel px = {xx,yy,adc};
    
    if(pre_clustId!=clustId) {
      mod_cl.first = pre_clustId;
      mod_cl.second = v;
      ClusterInfo.push_back(mod_cl);
      v.clear();
      pre_clustId = clustId;
    }
    v.push_back(px);
    i++;
  }
  // avoid reading last line twice
  v.pop_back();
  cout<<"Total pixels: "<<i<<endl;
  Index[total_cluster] = i-1;
  // to read last cluster
  mod_cl.first = pre_clustId;
  mod_cl.second = v;
  ClusterInfo.push_back(mod_cl);
  input.close();
  int total_hits = i;
  i = 0;
 // for each cluster call the cpe function
 using namespace chrono;
 ofstream cpeFile("CPE_CPU.txt");
 cpeFile<<"event\t moduleId\t clusterId\t xhit\t yhit"<<endl;
 high_resolution_clock::time_point t1 = high_resolution_clock::now();
 for(auto itm = ClusterInfo.begin(); itm!=ClusterInfo.end(); itm++) {  // for each module
   clustId = (*itm).first;
   if (clustId!=ClusterId[i++]) {
     cout<<"Cluster mismatched"<<clustId<<setw(14)<<ClusterId[i-1]<<endl;
     exit(EXIT_FAILURE);
   }
   uint moduleId = getModule(clustId); 
   uint event  = getEvent(clustId);
   std::vector<Pixel> pixVec = (*itm).second;
    
   LocalPoint lp = pixelHit(detDB,moduleId, pixVec);
   cpeFile<<setw(4)<<event<<setw(6)<<moduleId<<setw(16)<<clustId<<setw(20)<<lp.x()<<setw(20)<<lp.y()<<endl;
 }
 cpeFile.close();
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  double mic_sec = duration_cast<microseconds>(t2-t1).count();
  //cout<<"CPE CPU Time(mic. sec): "<<mic_sec<<endl;
}


ClusterParam getClusterParam(const std::vector<Pixel> &cluster) {
  int xmin = 0, ymin = 0;
  int xmax = 0, ymax = 0;
  int sizeX = 0, sizeY = 0;
  auto itc = cluster.begin();

  xmin = xmax = (*itc).row;
  ymin = ymax = (*itc).col;

  for (auto it=cluster.begin(); it!=cluster.end(); it++) {
    Pixel pixel = *it;
    if(xmin > pixel.row) xmin = pixel.row;
    if(xmax < pixel.row) xmax = pixel.row;

    if(ymin > pixel.col) ymin = pixel.col;
    if(ymax < pixel.col) ymax = pixel.col;
  }
  sizeX = xmax - xmin;
  sizeY = ymax - ymin;

  ClusterParam clustParam;

  clustParam.xmin = xmin;
  clustParam.xmax = xmax;
  clustParam.ymin = ymin;
  clustParam.ymax = ymax;
  clustParam.sizeX = sizeX;
  clustParam.sizeY = sizeY;
  return clustParam;
}

void 
collectCharge (const std::vector<Pixel> &cluster,
                   const ClusterParam  &clustParam,
                   float &Q_l_X, float &Q_f_X,
                   float &Q_l_Y, float &Q_f_Y ) 
{
  int xmin = clustParam.xmin, ymin = clustParam.ymin;
  int xmax = clustParam.xmax, ymax = clustParam.ymax;

  for(auto it=cluster.begin(); it!=cluster.end(); it++) {
    Pixel pixel = *it;
   // upper cut is put on pixel charge but does not affect the result much 
   float pix_adc = pixel.ADC;//pixel.ADC > 13200.0f ? 13200.0f:pixel.ADC;
    if((pixel.row==xmin)) Q_f_X+= pix_adc;
    if((pixel.row==xmax)) Q_l_X+= pix_adc;
    if((pixel.col==ymin)) Q_f_Y+= pix_adc;
    if((pixel.col==ymax)) Q_l_Y+= pix_adc;
  }
} 

LocalPoint localPositionInCm(float x, float y) {
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
  constexpr int bigYIndeces[]{0,51,52,103,104,155,156,207,208,259,260,311,312,363,364,415};
  auto const j = std::lower_bound(std::begin(bigYIndeces),std::end(bigYIndeces),binoffy);
  if (*j==binoffy) { local_pitchy  *= 2 ;}
  
  binoffy += (j-bigYIndeces);
  // The final position in local coordinates 
  float ycm = float(binoffy*pitchY) + fractionY*local_pitchy - m_yoffset;

  //float xcm = x*pitchX - m_xoffset;
  //float ycm = y*pitchY - m_yoffset;
  LocalPoint lp;
  lp.xcm=xcm;
  lp.ycm=ycm;
  return lp;
}   

LorentzAngle comcputeLorentzAngle(const int moduleId, const std::vector<Pixel>& clust, const float (&detDB)[NMODULE][NTUPLE]) {
  // compute cot alpha and beta for each cluster
  // formula to calculate cot alpha and cot beta as 
  // https://cmssdt.cern.ch/lxr/source/DataFormats/SiPixelCluster/interface/SiPixelCluster.h#0104
  float totalCharge = 0.0f;
  float xc=0.0f, yc=0.0f;
  for (auto it=clust.begin(); it!=clust.end(); it++) {
    Pixel pixel = *it;
    xc += (pixel.row + 0.5f)*pixel.ADC;
    yc += (pixel.col + 0.5f)*pixel.ADC;
    totalCharge += pixel.ADC;
  }
  xc = xc/totalCharge;
  yc = yc/totalCharge;
  //cout<<"xc: "<<xc<<"    yc: "<<yc<<endl;
  LocalPoint lp = localPositionInCm(xc, yc);
  //cout<<"lpx"<<lp.x()<<"     lpy: "<<lp.y()<<endl;
  auto gvx = lp.x() - detDB[moduleId][X0];
  auto gvy = lp.y() - detDB[moduleId][Y0];
  auto gvz = -1.f/detDB[moduleId][Z0];
  float cot_aplha = gvx*gvz;
  float cot_beta  = gvy*gvz;
  //printf("gvx: %f   x0: %f   gvy: %f   y0: %f   gvz: %f  z0: %f\n",
   // gvx,detDB[moduleId][X0],gvy,detDB[moduleId][Y0], gvz, detDB[moduleId][Z0] );
  LorentzAngle la = { cot_aplha, cot_beta }; 
  return la;
}
//-------------------------------------------------------------
// Return the BIG pixel information for a given pixel
//reference: http://cmslxr.fnal.gov/source/Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h?v=CMSSW_9_2_0#0119
bool isItBigPixelInX( const int ixbin ) {
  return (( ixbin == 79 ) || ( ixbin == 80 ));
}

bool isItBigPixelInY( const int iybin ) {
  int iybin0 = iybin%52;
  return(( iybin0 == 0 ) || ( iybin0 == 51 ));
     // constexpr int bigYIndeces[]{0,51,52,103,104,155,156,207,208,259,260,311,312,363,364,415,416,511};
     // return *std::lower_bound(std::begin(bigYIndeces),std::end(bigYIndeces),iybin) == iybin;
}
/*
float computeLorentShiftX(const int moduleId, const float (&detDB)[NMODULE][NTUPLE]) {
  return detDB[moduleId][LorentzShiftX]; 
}

float computeLorentShiftX(const int moduleId, const float (&detDB)[NMODULE][NTUPLE]) {
  return detDB[moduleId][LorentzShiftY]; 
}
*/
void  uploadDetDB(float (&det_DB)[NMODULE][NTUPLE]) {
  ifstream dbFile;
  dbFile.open("Pixel_CPE_Phase1_database_C_920.dat");
  if (!dbFile.is_open()) {
    cout<<"Error in opening file: detCoord_lorentzShift_database.txt"<<endl;
    exit(EXIT_FAILURE);
  }
  string str;
  getline(dbFile, str);
  int module, i=0;
  uint RawId;
  float x0, y0, z0, Rdet, Zdet, lshiftX, lshiftY;
  while(!dbFile.eof()) {
    dbFile>> module>> RawId>>x0>> y0>> z0>> Rdet >> Zdet >>lshiftX >>lshiftY;
    det_DB[i][X0] = x0;
    det_DB[i][Y0] = y0;
    det_DB[i][Z0] = z0;
    det_DB[i][LorentzShiftX] = lshiftX;
    det_DB[i][LorentzShiftY] = lshiftY;
    i++;
  }
  dbFile.close();
}
