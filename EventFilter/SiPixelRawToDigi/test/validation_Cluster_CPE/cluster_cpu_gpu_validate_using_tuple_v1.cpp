#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <tuple>
#include <algorithm>

// The code is successfully compiled on lxplus by 
// > g++ -O2 cluster_cpu_gpu_validate_using_tuple.cpp -std=gnu++11 -o out
// > ./out


using namespace std;
typedef unsigned long long uint64;
typedef tuple<uint64, int, int> cluster_tuple;

bool mycompare(const cluster_tuple &lhs, const cluster_tuple &rhs){
  //return get<1>(lhs) < get<1>(rhs);
  return get<0>(lhs) < get<0>(rhs);
}

int main(void){
  int icnt, icpu_mismatch, igpu_mismatch;
  int tot_cpu_clusters, tot_gpu_clusters;
  uint64 iprev_id, i0, j0;
  int i1, i2, j1, j2;
  uint64 ClusterID;
  int XX, YY;
  string str; // Have a buffer string
  vector<cluster_tuple> cpu_data;
  vector<cluster_tuple> gpu_data;
  cluster_tuple cpu_row, gpu_row;

  //Fill up Cluster_CPU tuple
  icnt=0;
  ifstream inCluster_CPU;  inCluster_CPU.open("Cluster_CPU.txt");
  getline(inCluster_CPU, str); //Header Line
  while(!inCluster_CPU.eof()) {
    inCluster_CPU >> ClusterID >> XX >> YY;
    icnt++;
    if(ClusterID==1) continue; // ignore the false cluster
    cpu_data.push_back(make_tuple(ClusterID, XX, YY));
  }
  inCluster_CPU.close();
  sort(cpu_data.begin(),cpu_data.end(), mycompare);


  //Fill up Cluster_GPU tuple
  icnt=0;
  ifstream inCluster_GPU;  inCluster_GPU.open("Cluster_GPU.txt");
  getline(inCluster_GPU, str); //Header Line
  while(!inCluster_GPU.eof()) {
    icnt++;
    inCluster_GPU >> ClusterID >> XX >> YY;
    if(ClusterID==0)continue; // ignore the bad cluster
    gpu_data.push_back(make_tuple(ClusterID, XX, YY));
  }
  inCluster_GPU.close();
  sort(gpu_data.begin(),gpu_data.end(), mycompare);

  cout << "Please Wait (~5 minutes): CPU Cluster data matching with GPU Cluster data is under progress ..." << endl;
  icnt=0; icpu_mismatch=0; tot_cpu_clusters=0; iprev_id=0;
  for(vector<cluster_tuple>::iterator iter_cpu=cpu_data.begin(); iter_cpu!= cpu_data.end(); iter_cpu++){
    cpu_row=*iter_cpu;
    i0=get<0>(*iter_cpu); i1=get<1>(*iter_cpu); i2=get<2>(*iter_cpu);
    if(iprev_id!=i0){iprev_id=i0; tot_cpu_clusters++;}
    icnt++;
    //cout << "CPU Data:  " << icnt << "  " << i0 << "\t" << i1 << "\t" << i2 << endl;
    vector<cluster_tuple>::iterator iter_gpu = std::find_if(gpu_data.begin(), gpu_data.end(), [&cpu_row](const cluster_tuple &e) {return e == cpu_row;});
    if(iter_gpu != gpu_data.end()) {
      j0=get<0>(*iter_gpu);  j1=get<1>(*iter_gpu); j2=get<2>(*iter_gpu);
    }else{
      icpu_mismatch++;
      if(icpu_mismatch==1){
        cout << endl << "Following CPU Entries Not Matching with GPU" << endl;
        cout << "      Row#      ClusterID       XX       YY" << endl;
      }
      cout << "CPU:  " << icnt << "\t" << i0 << "\t" << i1 << "\t" << i2 << endl;
      //exit(0);
    }
  }//for(vector<cluster_tuple>::iterator iter_cpu=cpu_data.begin(); ...
  cout << "Total Number of CPU Clusters=" << tot_cpu_clusters << endl;
  cout << endl << endl;
  if(icpu_mismatch==0){
    cout << "Matching CPU Cluster data with GPU Cluster data completed succesfully!!!";
    cout << endl << endl << endl;
  }


  cout << "Please Wait (~5 minutes): GPU Cluster data matching with CPU Cluster data is under progress ..." << endl;
  icnt=0; igpu_mismatch=0; igpu_mismatch=0;tot_gpu_clusters=0, iprev_id=0;
  for(vector<cluster_tuple>::iterator iter_gpu=gpu_data.begin(); iter_gpu!= gpu_data.end(); iter_gpu++){
    gpu_row=*iter_gpu;
    i0 = get<0>(*iter_gpu); i1 = get<1>(*iter_gpu); i2 = get<2>(*iter_gpu);
    icnt++;
    if(iprev_id!=i0){iprev_id=i0; tot_gpu_clusters++;}
    vector<cluster_tuple>::iterator iter_cpu = std::find_if(cpu_data.begin(), cpu_data.end(), [&gpu_row](const cluster_tuple &e) {return e == gpu_row;});
    if(iter_cpu != cpu_data.end()) {
      j0=get<0>(*iter_cpu); j1=get<1>(*iter_cpu); j2=get<2>(*iter_cpu);
    }else{
      igpu_mismatch++;
      if(igpu_mismatch==1){
        cout << endl << "Following GPU Entries Not Matching with CPU" << endl;
        cout << "      Row#      ClusterID       XX       YY" << endl;
      }
      cout << "GPU:  " << icnt << "\t" << i0 << "\t" << i1 << "\t" << i2 << endl;
      //exit(0);
    }
  }//for(vector<cluster_tuple>::iterator iter_gpu=cpu_data.begin(); ...
  cout << "Total Number of GPU Clusters=" << tot_gpu_clusters << endl;
  cout << endl << endl;
  if(igpu_mismatch==0)cout << "Matching GPU Cluster data with CPU Cluster data complted succesfully!!!" << endl << endl;
  if(icpu_mismatch==0 && igpu_mismatch==0)cout << "Congratulations: Validation of CPU and GPU Cluster data completed succesfully" << endl;
  cout<<"Total clusters: "<<max(tot_cpu_clusters, tot_gpu_clusters)<<endl;
  cout<<"Number of cluster matching: "<<min(tot_cpu_clusters, tot_gpu_clusters)<<endl;

}//int main(void){
