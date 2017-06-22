#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <tuple>
#include <algorithm>

using namespace std;

typedef tuple<int, int, float, float> cpe_tuple;

bool mycompare(const cpe_tuple &lhs, const cpe_tuple &rhs){
 
  return get<0>(lhs) < get<0>(rhs);
}

int main(void){
  
  std::vector<cpe_tuple> cpu;
  std::vector<cpe_tuple> gpu;

  int clusterId, module;
  float xhit, yhit;
  ifstream icpu("CPE_CPU_Standalone.txt");
  ifstream igpu("CPE_CPU_CMSSW.txt");
  string str;
  getline(icpu, str);
  getline(igpu,str);
  while(!icpu.eof()) {
    icpu>>module>>clusterId>>xhit>>yhit;
    cpu.push_back(make_tuple(module,clusterId,xhit,yhit));
  }
  icpu.close();
  while(!igpu.eof()) {
    igpu>>module>>clusterId>>xhit>>yhit;
    gpu.push_back(make_tuple(module,clusterId,xhit,yhit));
  }
  igpu.close();

  for(auto it_cpu =cpu.begin(), it_gpu=gpu.begin(); it_cpu!=cpu.end();it_cpu++, it_gpu++) {
    if(*it_cpu!=*it_gpu) {
      float xdiff = get<2>(*it_cpu) -get<2>(*it_cpu);
      float ydiff = get<3>(*it_gpu) -get<3>(*it_gpu);
      xdiff = abs(xdiff);
      ydiff = abs(ydiff);
      if(xdiff>0.00001 || ydiff>0.00001) {
      cout<<"Following entry did not match! at line:"<<cpu.begin()-it_cpu<<endl;
      cout<<"CPU: "<<get<0>(*it_cpu)<<"  "<<get<1>(*it_cpu)<<"  "<<get<2>(*it_cpu)<<"  "<<get<3>(*it_cpu)<<endl;
      cout<<"GPU: "<<get<0>(*it_gpu)<<"  "<<get<1>(*it_gpu)<<"  "<<get<2>(*it_gpu)<<"  "<<get<3>(*it_gpu)<<endl;
      }
    }
  }
  cout<<"CPE validation finished !"<<endl;
  return 0;
}
