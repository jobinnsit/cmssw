#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <tuple>
#include <algorithm>

// The code is successfully compiled on lxplus by 
// > g++ r2d_cpu_gpu_validate_using_tuple.cpp -std=gnu++11 -o out
// > ./out

//This is my first tuple code for comparing data obtained
//from Raw2Digi_CPU and Raw2Digi_GPU codes !!!

using namespace std;

typedef tuple<int, int, int, int, int> r2dtuple;

bool mycompare(const r2dtuple &lhs, const r2dtuple &rhs){
  return get<1>(lhs) < get<1>(rhs);
}

int main(void){
  int icnt,FedID, RawID, XX, YY, ADC;
  string str; // Have a buffer string
  int i0, i1, i2, i3, i4;
  int j0, j1, j2, j3, j4;
  vector<r2dtuple> cpu_data;
  vector<r2dtuple> gpu_data;
  r2dtuple cpu_row, gpu_row;


  icnt=0;
  ifstream inR2D_CPU;  inR2D_CPU.open("R2D_CPU.txt");
  getline(inR2D_CPU, str); //Header Line
  while(!inR2D_CPU.eof()) {
    inR2D_CPU >> FedID >> RawID >> XX >> YY >> ADC;
    cpu_data.push_back(make_tuple(FedID, RawID, XX, YY, ADC));
  }
  inR2D_CPU.close();
  sort(cpu_data.begin(),cpu_data.end(), mycompare);


  icnt=0;
  ifstream inR2D_GPU;  inR2D_GPU.open("R2D_GPU.txt");
  getline(inR2D_GPU, str); //Header Line
  while(!inR2D_GPU.eof()) {
    inR2D_GPU >> FedID >> RawID >> XX >> YY >> ADC;
    gpu_data.push_back(make_tuple(FedID, RawID, XX, YY, ADC));
  }
  inR2D_GPU.close();
  sort(gpu_data.begin(),gpu_data.end(), mycompare);

  //vector<r2dtuple>::iterator iter_gpu = std::find_if(gpu_data.begin(), gpu_data.end(), [](const r2dtuple &e) 
  //{return std::get<0>(e) == 1200 && std::get<1>(e) ==  303050776;});

  cout << "Please Wait (5 minutes): CPU data matching with GPU data is under progress ..." << endl;
  icnt=0;
  for(vector<r2dtuple>::iterator iter_cpu=cpu_data.begin(); iter_cpu!= cpu_data.end(); iter_cpu++){
    static int row=0;
	row++;
	cout<<"row: "<<row<<endl;
	cpu_row=*iter_cpu;
    i0 = get<0>(*iter_cpu);    i1 = get<1>(*iter_cpu);
    i2 = get<2>(*iter_cpu);    i3 = get<3>(*iter_cpu);
    i4 = get<4>(*iter_cpu);
    icnt++;
    if(i1<300000000 || i1>400000000)continue;
    //cout << "CPU Data:  " << icnt << "  " << i0 << "\t" << i1 << "\t" << i2;
    //cout <<  "\t" << i3 << "\t" << i4 << endl;
    vector<r2dtuple>::iterator iter_gpu = std::find_if(gpu_data.begin(), gpu_data.end(), [&cpu_row](const r2dtuple &e) {return e == cpu_row;});
    if(iter_gpu != gpu_data.end()) {
      //cout << "Found" << endl;
      j0 = get<0>(*iter_gpu);    j1 = get<1>(*iter_gpu);
      j2 = get<2>(*iter_gpu);    j3 = get<3>(*iter_gpu);
      j4 = get<4>(*iter_gpu);
    }else{
      cout << "Following CPU Entry Not Found, Therefore quitting immediately" << endl;
      cout << "CPU:  " << i0 << "\t" << i1 << "\t" << i2 <<  "\t" << i3 << "\t" << i4 << endl;
      exit(0);
    }
  }//for(vector<r2dtuple>::iterator iter_cpu=cpu_data.begin(); ...
  cout << "CPU data matching with GPU data!!!" << endl << endl;


  cout << "Please Wait (5 minutes): GPU data matching with CPU data is under progress ..." << endl;
  icnt=0;
  for(vector<r2dtuple>::iterator iter_gpu=gpu_data.begin(); iter_gpu!= gpu_data.end(); iter_gpu++){
    gpu_row=*iter_gpu;
    i0 = get<0>(*iter_gpu);    i1 = get<1>(*iter_gpu);
    i2 = get<2>(*iter_gpu);    i3 = get<3>(*iter_gpu);
    i4 = get<4>(*iter_gpu);
    icnt++;
    if(i1<300000000 || i1>400000000)continue;
    //cout << "CPU Data:  " << icnt << "  " << i0 << "\t" << i1 << "\t" << i2;
    //cout <<  "\t" << i3 << "\t" << i4 << endl;
    vector<r2dtuple>::iterator iter_cpu = std::find_if(cpu_data.begin(), cpu_data.end(), [&gpu_row](const r2dtuple &e) {return e == gpu_row;});
    if(iter_cpu != cpu_data.end()) {
      //cout << "Found" << endl;
      j0 = get<0>(*iter_cpu);    j1 = get<1>(*iter_cpu);
      j2 = get<2>(*iter_cpu);    j3 = get<3>(*iter_cpu);
      j4 = get<4>(*iter_cpu);
    }else{
      cout << "Following GPU Entry Not Found, Therefore quitting immediately" << endl;
      cout << "GPU:  " << i0 << "\t" << i1 << "\t" << i2 <<  "\t" << i3 << "\t" << i4 << endl;
      exit(0);
    }
  }//for(vector<r2dtuple>::iterator iter_gpu=cpu_data.begin(); ...
  cout << "GPU data matching with CPU data!!!" << endl << endl;
  cout << "Congratulations: Validation of CPU and GPU data completed succesfully" << endl;

return 0;

}//int main(void){


//Some more useful access to the tuple is given below
//Example of accessing tuple by index:     cpu_row = cpu_data[25];
//Example of assigning iterator to tuple:  cpu_row = *iter_cpu;
