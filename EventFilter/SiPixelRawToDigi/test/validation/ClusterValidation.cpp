#include <cstdio>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <malloc.h>
#include <set>
#include <iomanip>

using namespace std;
#define N 200000
#define MAX_SIZE 200000
#define EN 1
#define MOD 1000000

typedef struct  {
  int xx[MAX_SIZE];
  int yy[MAX_SIZE];
  int clusterId[MAX_SIZE];
  int total_module;
  int module[1856];
  
} ClusterInfo;

void hashCluster(ClusterInfo *CPUcluster,int counter) {
  static int call=0;
  call++;
  std::set<int> setOfcluster;
  for(int i=0;i<1856;i++) {
    CPUcluster->module[i] = 0;

  }
  int moduleId; 
  //printf("Inside hashCluster: %d\n",counter );
  setOfcluster.insert(CPUcluster->clusterId[0]);
  moduleId = CPUcluster->clusterId[0]/MOD;
  CPUcluster->module[moduleId] += 1;
 // printf("Before for loop\n");
  for(int i=1;i<counter;i++) {
    // Search for element in set using find member function
    std::set<int>::iterator it = setOfcluster.find(CPUcluster->clusterId[i]);
    if(it == setOfcluster.end())  {
      moduleId = CPUcluster->clusterId[i]/MOD;
      CPUcluster->module[moduleId] += 1;
      //pcluId = CPUcluster->clusterId[i];
      //if(moduleId==68) printf("%d\n",CPUcluster->clusterId[i]);
      setOfcluster.insert(CPUcluster->clusterId[i]);
    }
  }
  CPUcluster->total_module = moduleId;

}

int main(int argc, char* argv[]) {
  //printf("Inside main\n");
  ClusterInfo CPUcluster;
  ClusterInfo GPUcluster;
  
  int xx[N],yy[N],clusterId[N],IndEvent[EN], Event_size[EN];
  
  ifstream file ;
  ifstream file1 ;
  file.open("Cluster_CPU.txt");
  string line;
  getline(file, line);
  int i=0,e=0,p,x,y,cl;
  while(!file.eof()) {
    file >> cl  >> x >>y;
    xx[i] = x;
    yy[i] = y;
    clusterId[i] = cl;
    i++;
  }
  file.close();
 
  //printf("aftercpu i: %d\n",i);
    for(int j=0;j<i;j++) {
        //printf("index: %d\n",Index);
      CPUcluster.xx[j] = xx[j];
      CPUcluster.yy[j] = yy[j];
      CPUcluster.clusterId[j] = clusterId[j];
    }
   
    //assign_cluster(&CPUcluster, Event_size[ev]);
    hashCluster(&CPUcluster, i);
  
  // gpu 

  //printf("gpu start\n");
  i=0;
  file1.open("Cluster_GPU.txt");
  getline(file1, line);
  while(!file1.eof()) {
    file1 >>cl>>x>>y;
    xx[i] = x;
    yy[i] = y;
    clusterId[i] = cl;
    i++;
  }
  file1.close();
 
  // printf("after gpui: %d\n",i);

  for(int j=0;j<i;j++) {
     //printf("index: %d\n",Index);
    GPUcluster.xx[j] = xx[j];
    GPUcluster.yy[j] = yy[j];
    GPUcluster.clusterId[j] = clusterId[j];
  }
   
   // printf("before hash\n");
    hashCluster(&GPUcluster, i);
   // printf("after hash\n");

  ofstream outputfile;
  outputfile.open("GPUCPUClustercomparison.txt");
  outputfile<<"module_no\t No_cluster(gpu)\tNo_cluster(cpu)\tDifference"<<endl;
  for(int j=0;j<1856;j++)
    outputfile<<setw(6)<<j<<setw(6)<<GPUcluster.module[j]<<setw(6)<<CPUcluster.module[j]<<setw(6)
	<<GPUcluster.module[j]-CPUcluster.module[j]<<endl;
  
  outputfile.close();
 cout<<"Finished successfully!"<<endl;
 return 0;
}
