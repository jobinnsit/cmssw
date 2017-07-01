
#ifndef CPEGPU_H
#define CPEGPU_H

// common and constant parameters for all the DetId(modules)
// extracted from the CMSSW CPE
// all these parameters should be moved in python config
//const int NMODULE = 1856; // 1856 module
const int NTUPLE = 5;     // {X0, Y0, Z0, lorentshiftx, lorentShiftY}
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

typedef unsigned int uint;

struct Pixel {
  uint row;
  uint col;
  uint ADC;
};
struct LocalPoint {
  float xcm; //coordinate in cm
  float ycm; //coordinate in cm
  __device__ float x() {
    return xcm;
  }
  __device__ float y() {
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

struct DetDB {
  uint RawId[NMODULE+1];
  float X0[NMODULE+1];
  float Y0[NMODULE+1];
  float Z0[NMODULE+1];
  float LorentzShiftX[NMODULE+1];
  float LorentzShiftY[NMODULE+1];
};
DetDB *detDB; 

// global variable to store the 
float *xhit_d, *yhit_d;

struct CPE_cut_Param {
  // All these parameters to be intialzed from python conf
  // -----------------
  int size_cutX = 3;
  int size_cutY = 3; // to be moved into python conf
  float the_eff_charge_cut_lowX = 0.0 ;  // default
  float the_eff_charge_cut_lowY  = 0.0;
  float the_eff_charge_cut_highX = 1.0f ;// default 
  float the_eff_charge_cut_highY = 1.0f;
  // ----------------
};

void uploadCPE_db(DetDB *detDB);
__device__ LocalPoint localPositionInCm(const float x, const float y);
__device__ LorentzAngle computeLorentzAngle(const DetDB *detDB, const uint moduleId,
                        const uint startIndex, const uint size, const uint *xx,
                        const uint *yy, const uint *adc );
__device__ void min_max(uint startIndex,uint size,
                        const uint *xx, float &xmin,float &xmax);
__device__ void 
collectCharge (const uint *xx, const uint *yy, const uint *adc,
               uint startIndex, uint size, float xmin, float xmax,
               float ymin, float ymax, float &Q_l_X, float &Q_f_X,
               float &Q_l_Y, float &Q_f_Y );


__device__ float genericPixelHit(uint size, float first_pix, float last_pix,
                      float Q_f, float Q_l, float cot_angle, float pitch,
                      float theThickness, float lorentz_shift,
                      bool first_is_big, bool last_is_big,
                      float eff_charge_cut_low, 
                      float eff_charge_cut_high,
                      float size_cut);
__device__ bool isItBigPixelInX( const int ixbin );
__device__ bool isItBigPixelInY( const int iybin );
#endif // CPEGPU_H
