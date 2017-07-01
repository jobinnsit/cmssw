#ifndef CPECPU_H
#define CPECPU_H

#include <vector>

struct Pixel {
  unsigned int row;
  unsigned int col;
  unsigned int ADC;
};
struct LocalPoint {
  float xcm; //coordinate in cm
  float ycm; //coordinate in cm
  float x() {
  	return xcm;
  }
  float y() {
  	return ycm;
  }
};

struct ClusterParam {
  float xmin;
  float xmax;
  float ymin;
  float ymax;
  int sizeX;
  int sizeY;

};
struct LorentzAngle {
  float cotAlpha;
  float cotBeta;
};
const int NMODULE = 1856; // 1856 module
const int NTUPLE = 5;     // {X0, Y0, Z0, lorentshiftx, lorentShiftY}

LorentzAngle comcputeLorentzAngle(const int moduleId, const std::vector<Pixel>& clust, const float (&detDB)[NMODULE][NTUPLE]);

//float computeLorentShiftX(const int moduleId, const float (&detDB)[NMODULE][NTUPLE]);
//
//float computeLorentShiftY(const int moduleId, const float (&detDB)[NMODULE][NTUPLE]);

LocalPoint localPositionInCm(const float x, const float y);
LocalPoint pixelHit(const int module, const std::vector<Pixel> &v, const float (&detDB)[NMODULE][NTUPLE]);

float genericPixelHit(int size, float first_pix, float last_pix,
                      float Q_f, float Q_l, float cot_angle, float pitch,
                      float thickness, float lorentz_shift,
                      bool is_first_big, bool is_last_big,
                      float eff_charge_cut_low, 
                      float eff_charge_cut_high,float size_cut);

void 
collectCharge (const std::vector<Pixel> &cluster,
                   const ClusterParam  &clustParam,
                   float &Q_l_X, float &Q_f_X,
                   float &Q_l_Y, float &Q_f_Y );
ClusterParam getClusterParam(const std::vector<Pixel> &cluster);
void  uploadDetDB(float (&det_DB)[NMODULE][NTUPLE]);
bool isItBigPixelInX(int pixel);
bool isItBigPixelInY(int pixel);
void  uploadDetDB(float (&det_DB)[NMODULE][NTUPLE]);
// pitch size aligned with C_9_2_0 
const float pitchX = 0.01; // in cm instead of reading from database we have hard coded
const float pitchY = 0.015; // in cm 

// thickness for forward module has changed to 0.029
const float thicknessBarrel  = 0.0285;  // in cm instead of reading from database we have hard coded
const float thicknessForward = 0.0290;  // in cm
// thickness for forward module has changed to 0.029
// widthLAFractionX and Y has changed to 1 in CMSSW_9_2_0
const float widthLAFractionY = 1; // for both detector type
const float widthLAFractionX_Barrel = 1;  // barrel detector
const float widthLAFractionX_Forward = 1; // forward detector

float detDB[NMODULE][NTUPLE]; // contains database for X0, Y0, R_Det, lorentShiftX, lorentShiftY
enum DetParam{
  X0, 
  Y0, // x0, y0 are middle of the detector in the local frame in cm
  Z0, // R_Det in paper, radial detector's coordinate in cm
  LorentzShiftX,
  LorentzShiftY
};

#endif // CPECPU_H
