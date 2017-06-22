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

#include "CPECPU.h"

using namespace std;

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
                      bool first_is_big, bool last_is_big, float eff_charge_cut_low, 
                      float eff_charge_cut_high, float size_cut) { //charge_cut_high_x, charge_cut_low_x to be included

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


int main(int argc, char const *argv[]) {
  
  // upload the lorentzshift and x0 and y0 and R_det for all
  // module which is constant 
  uploadDetDB(detDB);

  // to hold the input data
  typedef std::vector<std::vector<Pixel>> ClusterVec; 
  typedef std::pair<int, ClusterVec> ModuleSet;

  ClusterVec clusVec;
  ModuleSet  mod_ClusVec;  
  
  std::vector<ModuleSet> theModule;
  
  std::set<int> modSet;  
  std::set<int> clusterSet;

  std::set<int>::iterator im;
  std::set<int>::iterator ic;
  
  unsigned int module, old_module=0, clusterId, hit, x, y, adc;
  ifstream cpeinputFile;
  //cpeinputFile.open("../data/CPEInput_Forcpu.txt");
  cpeinputFile.open("CPE_Input_CPU.txt");
  string firstline;
  getline(cpeinputFile, firstline);
  //getline(cpeinputFile, firstline);
  #ifdef DEBUG
  cout<<"input file opened";
  int count =0;
  #endif
  std::vector<Pixel> v;
  modSet.insert(0);     // module id starts at 0
  clusterSet.insert(1); // cluster id starts at 1
  // read the cluster of all the module
  while(!cpeinputFile.eof()) {
    cpeinputFile>> module >> clusterId >> hit >> x >> y >>adc;
    //#ifdef DEBUG
    //count++;
    //#endif
    // store pixels of a cluster
    ic = clusterSet.find(clusterId);
    if(ic==clusterSet.end()) {
      clusterSet.insert(clusterId);
      clusVec.push_back(v);
      v.clear(); 
    }
    // store clusters of a module
    im = modSet.find(module);
    if(im==modSet.end()) {
      
      // insert the last cluster
      clusVec.push_back(v);
      v.clear(); 
      modSet.insert(module);
      mod_ClusVec.first = old_module;
      mod_ClusVec.second = clusVec;
      theModule.push_back(mod_ClusVec);
      clusVec.clear();
      clusterSet.clear();
      clusterSet.insert(1);
      old_module = module;
      //cout<<"module: "<<module<<endl;
    }
    Pixel px = {x,y,adc};
    v.push_back(px);
  }
  clusVec.push_back(v);
  mod_ClusVec.first = module;
  mod_ClusVec.second = clusVec;
  theModule.push_back(mod_ClusVec);

  #ifdef DEBUG
     cout<<"count: "<<count<<endl;
  #endif
 ofstream cpeOut("CPE_CPU_Standalone.txt");
 cpeOut<<"moduleId\t\t clusterId\t\t xhit \t yhit"<<endl;
 for (auto itm = theModule.begin(); itm!=theModule.end(); itm++) {  // for each module
     int moduleId       = (*itm).first;
     //if(moduleId==1) break;
     #ifdef DEBUG
     cout<<"moduleId: "<<moduleId<<endl;
     #endif
     ClusterVec clusVec = (*itm).second;
     int  clusterId = 1;
     for (auto itc=clusVec.begin(); itc!=clusVec.end();itc++) {  // for each cluster
      //auto pixVec = (*itc).begin();
      auto pixVec = *itc;
      LocalPoint lp = pixelHit(detDB,moduleId, pixVec);
      cpeOut<<setw(4)<<moduleId<<setw(14)<<clusterId++<<setw(14)<<lp.x()<<setw(14)<<lp.y()<<endl;
      //cout<<lp.x()<<"    "<<lp.y()<<endl; 
      //static int count =1;
      //count++;
      //if(count==5) break;

     }
  }
  cpeOut.close();

  return 0;
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
  LocalPoint lp = localPositionInCm(xc, yc);

  auto gvx = lp.x() - detDB[moduleId][X0];
  auto gvy = lp.y() - detDB[moduleId][Y0];
  auto gvz = -1.f/detDB[moduleId][Z0];
  float cot_aplha = gvx*gvz;
  float cot_beta  = gvy*gvz;
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
  dbFile.open("Pixel_CPE_Phase1_database.dat");
  if (!dbFile.is_open()) {
    cout<<"Error in opening file: detCoord_lorentzShift_database.txt"<<endl;
    exit(EXIT_FAILURE);
  }
  string str;
  getline(dbFile, str);
  int module, i=0;
  float x0, y0, z0, Rdet, Zdet, lshiftX, lshiftY;
  while(!dbFile.eof()) {
    dbFile>> module>> x0>> y0>> z0>> Rdet >> Zdet >>lshiftX >>lshiftY;
    det_DB[i][X0] = x0;
    det_DB[i][Y0] = y0;
    det_DB[i][Z0] = z0;
    det_DB[i][LorentzShiftX] = lshiftX;
    det_DB[i][LorentzShiftY] = lshiftY;
    i++;
  }
  dbFile.close();
}
