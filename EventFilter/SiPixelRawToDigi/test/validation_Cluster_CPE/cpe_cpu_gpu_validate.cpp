#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <tuple>
#include <algorithm>

using namespace std;
typedef unsigned long long uint64;
typedef tuple<int,int, uint64, double, double> cpe_tuple;

int main(void){
  
  std::vector<cpe_tuple> cpu;
  std::vector<cpe_tuple> gpu;

  int  module,event;
  uint64 clusterId;
  float xhit, yhit;
  ifstream icpu("CPE_CPU.txt");
  ifstream igpu("CPE_GPU.txt");
  string str;
  getline(icpu, str);
  getline(igpu,str);
  while(!icpu.eof()) {
    icpu>>event>>module>>clusterId>>xhit>>yhit;
    cpu.push_back(make_tuple(event,module,clusterId,xhit,yhit));
  }
  cout<<"Total CPE hit: "<<cpu.size()-1<<endl;;
  icpu.close();
  while(!igpu.eof()) {
    igpu>>event>>module>>clusterId>>xhit>>yhit;
    gpu.push_back(make_tuple(event,module,clusterId,xhit,yhit));
  }
  igpu.close();

  for(auto it_cpu =cpu.begin(), it_gpu=gpu.begin(); it_cpu!=cpu.end()-1;it_cpu++, it_gpu++) {
    if(*it_cpu!=*it_gpu) {
      double xdiff = get<3>(*it_cpu) -get<3>(*it_gpu);
      double ydiff = get<4>(*it_gpu) -get<4>(*it_gpu);
      xdiff = abs(xdiff);
      ydiff = abs(ydiff);
      
      if(xdiff>0.00001 || ydiff>0.00001) {
      cout<<"Following entry did not match! line no:"<<it_cpu-cpu.begin()<<endl;
      cout<<"CPU: "<<get<0>(*it_cpu)<<"  "<<get<1>(*it_cpu)<<"  "<<get<2>(*it_cpu)<<"  "<<get<3>(*it_cpu)<<"  "<<get<4>(*it_cpu)<<endl;
      cout<<"GPU: "<<get<0>(*it_gpu)<<"  "<<get<1>(*it_gpu)<<"  "<<get<2>(*it_gpu)<<"  "<<get<3>(*it_gpu)<<"  "<<get<4>(*it_gpu)<<endl;
      }
    }
  }
  cout<<"CPE validation finished !"<<endl;
  return 0;
}
