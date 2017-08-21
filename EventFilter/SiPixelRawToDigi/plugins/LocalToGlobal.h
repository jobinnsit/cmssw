#ifndef LOCALTOGLOBAL_H
#define LOCALTOGLOBAL_H

// this structure may change little bit in future
// the output formate is not finlized
// we may discard some of the members if not needed 
// output structure
struct RecHit {
  bool  barrel =0; //true if barrel
  int   layer=0;
  int   disk=0;
  //HitID stores info about eventno and moduleno
  long long unsigned HitId;
  float lx;
  float ly;
  float x;
  float y;
  float z;
  float u;  //barrel: u=r, v=z, forward the opposite...
  float v;
  float phi;
  float theta;
};

//In barrel layers, each phi strip  has the same rotation matrix 
//which is hardocded 
//Each module in the forward pixel detector has its own rotation matrix
//which is hardcoded 
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

//Store the position of each module in the global frame.
//Store these constants in device memory needed for 
//local to global converson 
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

void uploadGlobal_Positon_Rotation_Matrix(GlobalPosition *globalPosRot);
#endif // CPEGPU_H
