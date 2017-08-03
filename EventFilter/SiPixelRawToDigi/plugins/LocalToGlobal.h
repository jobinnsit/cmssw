#ifndef LOCALTOGLOBAL_H
#define LOCALTOGLOBAL_H

// this structure may change little bit in future
// the output formate is not finlized
// we may discard some of the members if not needed 
struct RecHit {
  long long unsigned HitId;// will store info about eventno, moduleno, layer etc.
  float lx;
  float ly;
  uint layer;
  float x;
  float y;
  float z;
  float u;  // barrel: u=r, v=z, forward the opposite...
  float v;
  float phi;
  float theta =0;
};

// since we don't know the definition of Rotatation matrix[3,3].
// we have hardcoded the elements of rotation matrix for each module
// which are sorted in ascending order.
struct Rotation {
  float R11;
  float R12;
	float R13;
	float R21;
	float R22;
	float R23;
	float R31;
	float R32;
	float R33;
};

// store the position of each module in global frame.
// store these constants in device memory for local to global converson 
struct GlobalPosition {
	unsigned int RawId; //unsigned 32 bit DetId or RawId
	float xpos;
	float ypos;
	float zpos;
	float r;
	float phi;
	Rotation Rot;                                                                  
};
GlobalPosition *globalPosRot; // pointer to store the array of rotatidon matrix for each module
RecHit *Hit; // pointer to store the array of global hits

const int NMODULE = 1856;

void uploadGlobal_Positon_Rotation_Matrix(GlobalPosition *globalPosRot);
#endif // CPEGPU_H
