/*
* @Author: sushildubey171
* @Date:   2017-07-29 17:04:19
* @Last Modified by:   sushil
* @Last Modified time: 2017-08-02 12:17:00
* @Desc: Validate the output of local to global cooridnate conversion of GPU program
* with standalone CPU program
* @Input: GlobalHit_CPU_standalone.txt, GlobalHit_GPU_standalone.txt
* @Output: if the difference between corresponding x, y and z hit is greater than 0.5 micron
* than it is displayed on the screen 
*/

#include <iostream>
#include <map>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <cmath>


using namespace std;
typedef long long unsigned uint64;
struct Hit {
  double x;
  double y;
  double z;
  Hit(double _x, double _y, double _z): x(_x), y(_y), z(_z){}
};

bool operator==(const Hit &hit1,const Hit &hit2) {
  if( fabs(hit1.x-hit2.x) > 0.00005 || 
      fabs(hit1.y-hit2.y) > 0.00005 ||
  	  fabs(hit1.z-hit2.z) > 0.00005
  	)
  	return false;
  else return true;
  }
const int NMODULE = 1856;

void readInput(const string filename, std::map<uint64 ,Hit> &map) {
  ifstream ifile(filename);
  if (!ifile) {
   	cout<<"\nfile not found: "<<filename<<endl;
	exit(-1);
  }
  string line;
  getline(ifile, line);
  uint64 hitId;
  double lx, ly, lz, gx, gy,gz;
  for(int i=0;i<NMODULE; i++) {
    ifile>>hitId>>lx>>ly>>gx>>gy>>gz;
    map.insert(make_pair(hitId, Hit(gx,gy,gz)));
  }
  ifile.close();
}

void compare(const std::map<uint64, Hit> &m1, const std::map<uint64, Hit> &m2) {
  cout<<fixed;
  cout<<setprecision(6);
  cout<<"\nMismatch hits (if any)"<<endl;
  cout<<"hitId     cpux     cpuy     cpuz    gpux     gpuy     gpuz";
  cout<<"     delta_x(in micron)  delta_y    delta_z"<<endl;
  int counter =0;
  for(auto it1=m1.begin(), it2=m2.begin(); it1!=m1.end(); it1++, it2++ ) {
  	Hit gpuhit = (*it1).second;
  	Hit cpuhit = (*it2).second;
  	if(gpuhit==cpuhit) {} // if difference is greater than 0.5 micron than print
  	else {
  		cout<<(*it1).first<<setw(15)<<gpuhit.x<<setw(11)<<gpuhit.y<<setw(11)<<gpuhit.z<<setw(15)
  		    <<cpuhit.x<<setw(11)<<cpuhit.y<<setw(11)<<cpuhit.z<<setw(15)<<fabs(gpuhit.x-cpuhit.x)*1e4
  		    <<setw(11)<<fabs(gpuhit.y-cpuhit.y)*1e4<<setw(11)<<fabs(gpuhit.z-cpuhit.z)*1e4<<endl;
  		counter++;    
  	}	
  }

  cout<<"\nTotal mismatch values: "<<counter<<endl;
  cout<<"\nValidation finished"<<endl;

} 

int main(){
  
  const string gpucmssw = "GlobalHit_GPU_CMSSW.txt";
  const string cpustandalone = "GlobalHit_CPU_standalone.txt";
  std::map<uint64, Hit> m1; // for cmssw
  std::map<uint64, Hit> m2; // for cpu standalone
  readInput(gpucmssw, m1);
  readInput(cpustandalone, m2);
  // display the output on console if the differnce is greater than 0.5 micron
  compare(m1, m2);
  return 0;
}
