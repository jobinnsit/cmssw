/** File Name: GPUCluster_cc5.cu
* Date: 01/03/2017
* Clusterizer Algorithm which takes input from RawToDigi
* and find clusters and subclusters.
* The output of RawToDigi servers input to this program.
* Input: XX[], YY[], mIndexStart[], mIndexEnd[]
* XX: 0 - 159 rows
* YY: 0- 415  cols (source CMSSW ReadPixelCluster.cc)
* origin is shifted by 1 along both the axes
* Output: Cluster number associated with each pixels.
*--------Change Log ---------
* Date 06/03/2017
* produces clusters and sub cluster each part of code is tested by 
* printing the output each after cluster, after sorting and after
* unique operation and sub-cluster-kernel tested and working fine
* The illegal memory problem was due to the last block index in sub-cluster
* kernel which was garbage value, temp solution skip last block
*--------Change Log-----------
* Date 11/04/2017
* ADC array is added, now the cluster will also contain the adc of each pixel
* 
*-------change log--------------
* Date 27/04/2017
* Fixed bug in sub-clustering, if X- coordinate differs by more
* than blockDim.x then extra clusters were formed. Now it is fixed 
* by iterating upto xmax limit.
*/ 
// System includes
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <assert.h>
#include <iomanip>
#include <string>
#include <iostream>
#include <fstream>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/execution_policy.h>
// CUDA runtime
#include <cuda.h>
#include <cuda_runtime.h>
#include "CudaError.h"
#include "PixelCluster.h"

using namespace std;
using namespace thrust;

/*
  This functions checks for cuda error
  Input: debug message
  Output: returns cuda error message
*/
void checkCUDAError(const char *msg) {
  cudaError_t err = cudaGetLastError();
  if( cudaSuccess != err) {
    fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
    exit(-1);
  }
}

/*
  The origin of (x,y) was shifted by (1,1) for clustering.
  Reshift the origin to (0,0)
  Input: x[], y[], size
  Output: x[]=x[]-1, y[]=y[]-1; 
*/
__global__ void shift_origin_kernel(uint wordCounter,uint *xx, uint *yy) {
  
  uint gIndex = threadIdx.x + blockIdx.x*blockDim.x;
  if(gIndex<wordCounter) {
    //since bad pixel has x=0,y=0
    if(xx[gIndex]>0) { // either both are 0 or none
      xx[gIndex] = xx[gIndex] -1;
      yy[gIndex] = yy[gIndex] -1;
    }
  }
}
/*
  This kernel sorts the xx[] and yy[] as per the sorted Index[]
  Input: Index[], xx[], yy[]
  Output: xx[], yy[]
*/
__global__ void copy_kernel(const uint *Index,const uint *xx, const uint *yy,
                            const uint *ADC, const uint size, uint *xx1,
                            uint *yy1, uint *ADC1) {
  uint tid = blockDim.x * blockIdx.x + threadIdx.x;
    if(tid < size) {
    xx1[tid] = xx[Index[tid]];
    yy1[tid] = yy[Index[tid]];
    ADC1[tid] = ADC[Index[tid]];
    }
}

/*
  This kernel will check for the subdivision and forms the new cluster
  Input: Index[], xx[], yy[], number of pixels
  Output: gClusterId[]

*/
__global__ void sub_cluster_kernel80(const uint *Index, const uint *xx, 
                                     const uint *yy,const uint gpCounter,
                                     uint *gClusterId) 
{
  uint tid     = threadIdx.x;
  uint blockid = blockIdx.x;
  uint start   = Index[blockid];
  uint end     = Index[blockid+1];
  __shared__ uint old_clusterId, moduleId;
  __shared__ int nstripx,nstripy;

  //if(tid==0 && blockid<100) {
    //printf("sc-kernel blockId: %d  totalBlock:  %d\n", blockid, gridDim.x);
  //}
  // skip the empty clusters
  if(gClusterId[start] == 0) return;
  if(end==(start+1) ) return;

  // sub-cluster kecd te  
  
  //kernel to handle cluster with size <= 80
  // found after analysing around 300 events
  if(end-start <= 80) {

    // assuming that cluster size is less than 80
    __shared__ uint xp[162],yp[417];
    __shared__ int tid0x,tid0y;
    __shared__ uint xmin, ymin, xmax, ymax;
    nstripx=0;  nstripy=0;
    tid0x=-1;   tid0y=-1;
    
    #pragma unroll    // tells the compiler to unroll the loop 
    for(uint i=0;i<6;i++) {
      uint gtid = 6*tid + i;  // gtid upto 416 from 80 threads.
      if(gtid<MAX_Y)
        yp[gtid] = 0;
    }
        
    // intialize xp[0] to xp[160] from 80 threads.
    xp[2*tid+ 0]  = 0;
    xp[2*tid + 1] = 0; 
    // find xmin, ymin and subtract from all the pixels
    if(tid==0) {
      xmin = xx[start];
      xmax = xx[start];
      for(int i=1;i<(end-start);i++) {  
        if(xmin > xx[start+i]) {
          xmin = xx[start+i];
        }
        if(xmax<xx[start+i]) {xmax=xx[start+i];}
      }
    }
    if(tid==1) {
      ymin = yy[start];
      for(int i=1;i<(end-start);i++) {  
        if(ymin >yy[start+i]) {
          ymin = yy[start+i];
        }
        if(ymax<yy[start+i]) {ymax=yy[start+i];}
      }
    }
    if(tid==2) {
      xp[160] =0;
      xp[161] =0;
      yp[416] =0;
      old_clusterId = gClusterId[start];
      moduleId = old_clusterId / 1000000; 
    }
    __syncthreads();  
  
    // to find the projection
    if((start+tid) < end)  {
      uint xc = xx[start + tid] - (xmin-1);
      uint yc = yy[start + tid] - (ymin-1);
      xp[xc] = xc;
      yp[yc] = yc;                                    
    }
    __syncthreads();
  
    if(tid && (nstripx==0 || nstripy==0)){
      if(xp[tid] && !xp[tid - 1]){
        nstripx=1; tid0x=tid;
      }
      if(yp[tid] && !yp[tid - 1]){
        nstripy=1; tid0y=tid;
      }  
    }  // end of if(tid)
    __syncthreads();

    if(tid != tid0x) {
      if(xp[tid] && !xp[tid - 1]){
        nstripx=2;
      }    
    }  // end of if(tid)

    if(tid != tid0y) {
      if(yp[tid] && !yp[tid - 1]) {
        nstripy=2;
      }
    }       
    __syncthreads();
  
    if(nstripx==2 || nstripy==2){ // form cluster only if it is divisible
      if(tid) {
        if(xp[tid] && !xp[tid - 1]) {
          uint i = tid;
          while(xp[i]) {
            xp[i] = tid;
            i++;
            if(i == 162) break;
          }
        }
        if(yp[tid] && !yp[tid - 1]) {
          uint i = tid;
          while(yp[i]) {
            yp[i] = tid;
            i++;
            if(i == 416) break;
          }
        }
      } //end of if(tid)
      // if the difference is grater that 80 then they are not 
      // covered by the thread (rarest case)
      if(xmax-xmin>=blockDim.x) { // if the difference is greater than blockDim.x
        uint ext_tid = blockDim.x + tid;
        if(xp[ext_tid] && !xp[ext_tid-1]) {
          uint i = ext_tid;
          while(xp[i]) {
            xp[i] = ext_tid;
            i++;
            if (i==162) break;
          }
        }
      }

      __syncthreads();

      //assign the cluster id to each pixel
      if((start + tid )< end) {
        uint px=0, py=0,new_ClusterId=0;
        px = xp[xx[start + tid]-(xmin-1)] + xmin-1;            
        py = yp[yy[start + tid]-(ymin-1)] + ymin-1;         
        new_ClusterId = 1000000*(moduleId) + 1000 * py + px;
        if(old_clusterId!=new_ClusterId && (px!=0 && py!=0))
          gClusterId[start+tid] = new_ClusterId;  
      }
    }//end of if(nstripx==2 || nstripy==2) 
  } // end of if(end-start<=80)
  else if(blockid!=gridDim.x-1){  // sub-cluster kernel to handle cluster with size > 80
    //f(tid==0) printf("inside else blockId: %d\n",blockid );
    __shared__ uint xp[162],yp[417];
    __shared__ uint cluster_size,itrn;
    nstripx=0;  nstripy=0;
    #pragma unroll    // tells the compiler to unroll the loop 
    for(uint i=0;i<6;i++) {
      uint gtid = 6*tid + i;  // gtid upto 416 from 80 threads.
      if(gtid<MAX_Y)
        yp[gtid] = 0;
    }
        
    // intialize yp[0] to yp[160] from 80 threads.
    xp[2*tid+ 0]  = 0;
    xp[2*tid + 1] = 0;  
   
    if(tid==0) {
      xp[160] = 0;
      cluster_size = end-start;
      itrn = cluster_size/blockDim.x + 1;
      old_clusterId = gClusterId[start];
      moduleId = old_clusterId / 1000000;
    }
    //if(tid==0) printf("before intializatio in else block\n");
    __syncthreads();  
    //if(tid==0) printf("after intializatio in else block\n");
    //  to find the projection
    #pragma unroll 
    for(uint i=0;i<itrn;i++) {
      uint gtid  = 80*i + tid;
      if(gtid < cluster_size)  {
        uint xc = xx[start+gtid];
        uint yc = yy[start+gtid];
        xp[xc] = xc;
        yp[yc] = yc;                                                                            
      }
      __syncthreads();
    }
  //  if(tid==0) printf("after projection in else block\n");
    #pragma unroll
    for(uint j =0;j<6; j++) {  // generates gtid form 0.....416 from 80 threads
      uint gtid = 6*tid + j;
      if(gtid && gtid<MAX_Y) {
        if(yp[gtid] && !yp[gtid - 1]) {
          uint i = gtid;
          uint old_nstripy = atomicAdd(&nstripy,1);  // to store the no of xstrip
          while(yp[i]) {
            yp[i] = gtid;
            i++;
            if(i >= MAX_Y) break;
          }
        }
        if(gtid < MAX_X) {
          if(xp[gtid] && !xp[gtid - 1]) {
            uint k = gtid;
            uint old_nstripx = atomicAdd(&nstripx,1);
            while(xp[k]) {
              xp[k] = gtid;
              k++;
              if(k >= MAX_X ) break;
            }
          }
        }//end of if(gtid<161)
      }
      __syncthreads();
    } // end of for loop  
  // if(tid==0) printf("after strip in else block\n");
    //assign the cluster id to each pixel
    // only if it divisible 
    if(nstripx>1 || nstripy>1) {
      for(uint i=0;i<itrn;i++) {
        uint gtid = blockDim.x*i + tid;
        if(gtid < cluster_size)  {
          uint px=0, py=0,new_ClusterId=0;
          px = xp[xx[start+gtid]];                    // find location of pixel on x strip
          py = yp[yy[start+gtid]];                    // find location of strip on y strip
          new_ClusterId = 1000000*(moduleId) + 1000 * py + px;
          if(old_clusterId!=new_ClusterId && (px!=0 && py!=0))
            gClusterId[start+gtid] = new_ClusterId;    
        }
      }
    }
  // if(tid==0) printf("after gClusterId in else block\n");
  } // end of else

} // end of sub_cluster

// fills the Index[] with 0, 1,2,3.. upto the size of event
__global__ void createIndex_kernel(const uint wordCounter, uint *Index) {
  uint tid = blockDim.x * blockIdx.x + threadIdx.x;
    if(tid < wordCounter) {
      Index[tid] = tid;
    }
}
void createIndex(const uint wordCounter, uint *Index) {
  // to fill the index array Index[] 
  // which will be used in sort by key 
  int nthreads = 1024;
  int nblocks = wordCounter/nthreads +1;
  createIndex_kernel<<<nblocks, nthreads>>>(wordCounter, Index);
  cudaDeviceSynchronize();
  checkCUDAError("Error in createIndex_kernel");
}

/*
  This function sorts the cluster id, removes the duplicate
  and finds the sub-cluster within each cluster
  Input: d_Index[], d_gClusterId[], d_xx[], d_yy[], wordCounter, h_Index[]
  Output:d_xx1[], d_yy1[]
*/
void sub_cluster(uint *d_xx, uint *d_yy,const uint *d_ADC, uint *d_Index, uint *d_gClusterId, 
                  const uint wordCounter, uint *d_gClusterId1, uint *d_xx1, uint *d_yy1, uint *d_ADC1 ) 
{
  
  //cout<<"Inside sub_cluster function: wordCounter"<<endl;
  createIndex(wordCounter, d_Index);
  
  // get device_ptr needed for thrust operations
  thrust::device_ptr<uint> Index(d_Index); //Index is the index array
  thrust::device_ptr<uint> gClusterId(d_gClusterId);

  // sort the cluster id by key
  thrust::sort_by_key(gClusterId , gClusterId + wordCounter, Index);
 
  cudaMemcpy(d_gClusterId1, d_gClusterId, wordCounter * sizeof(uint),  cudaMemcpyDeviceToDevice );

  // launch kernel for sorting xx[] and yy[]
  uint N_threads = 1024;
  uint N_blocks  = wordCounter / N_threads +1;
  copy_kernel<<<N_blocks,N_threads>>>(d_Index, d_xx, d_yy,d_ADC, wordCounter,d_xx1,d_yy1,d_ADC1); 

  checkCUDAError("Error in copy kernel");
  
  // removes the consecutive duplicate
  // new_end.first gives size of gClusterId with no ducplicate
  uint total_cluster=0;
  thrust::pair<thrust::device_ptr<uint>, thrust::device_ptr<uint> > new_end;
  
  // Fill the index again which will be used for uniuqe 
  createIndex(wordCounter, d_Index);
  
  new_end = thrust::unique_by_key(gClusterId , gClusterId + wordCounter, Index );
  total_cluster = new_end.first - gClusterId;
  checkCUDAError(" Failed after unique operation");
 
  //launch the kernel for subdivision
  dim3 no_threads =  80; // maximum size of cluster found after analysis
  dim3 no_blocks  =  total_cluster;
  //cout<<"Total_clusters: "<<total_cluster<<endl;
  // Ignore first few cluster as they might contain 0s
  sub_cluster_kernel80<<<no_blocks,no_threads>>>(d_Index, d_xx1, d_yy1,wordCounter, d_gClusterId1);
  cudaDeviceSynchronize();
  checkCUDAError(" Failed after sub-cluster-kernel");

} // End of sort_cluster

/* 
  This is the main kernel for clustarisation
  Inputs:  xx[],yy[],module[], Imodule[]
  Outputs: gClusterId[]
*/
__global__ void cluster_kernel(uint *xx, uint *yy, const int *mIndexStart,
                               const int *mIndexEnd, uint *gClusterId) 
{

  __shared__ uint xp[MAX_X+1], yp[MAX_Y+1];   // Array to store x and y projection
  uint moduleId = blockIdx.x;                 // to get block id
  uint tid = threadIdx.x;                     // to get thread id
  uint moduleBegin = mIndexStart[moduleId];
  uint moduleEnd   = mIndexEnd[moduleId];
  if(tid==0) {
    //printf("moduleId: %d  mstart: %d  mend: %d\n",moduleId, moduleBegin, moduleEnd);
  }
  // Module empty 
  if((moduleBegin==0 && moduleEnd==0 ) ) {
    return;
  }
  //module contains only one pixel
  if(moduleBegin==moduleEnd) {
    int px = xx[moduleBegin];
    int py = yy[moduleBegin];
    gClusterId[moduleBegin]= 1000000 * (moduleId) + 1000 * py + px;
    return;
  }
  
  uint module_size = (moduleEnd - moduleBegin) + 1;

  yp[tid] = 0;                // Initialize X projection to false
  if(tid < MAX_X+1) {
    xp[tid] = 0;              // Initialize Y projection to false
  }
  
  __syncthreads();            // let all threads finish intialization
  // To get projection on x and y axis
  // For loop is used to deal with module with hits more
  // than blockDim.x
  uint noItr = module_size/blockDim.x+1;
  for(uint i=0; i<noItr; i++) {
    uint idx = moduleBegin + tid + i*blockDim.x;
    if(idx <= moduleEnd) {
      uint tx = xx[idx];
      uint ty = yy[idx];
      xp[tx] = tx;
      yp[ty] = ty;
      gClusterId[idx] = 0;
    } // if(tid < module_size)
  } // End of for(uint i=0; i<noItr...
 //printf("After P:\n"); 
  __syncthreads();

  //  Store the unique strip# for all element in one strip.
  //  Divide the projections array into multiple strips.
  //  Distribute each strip to one thread.
  //  Each thread will process the strip independently
  //  and will store the start location of strip at all location of strip
  //  This procedure will be repeated for both x and y strip.
  if(tid) {
    if(yp[tid] && !yp[tid - 1]) {
      uint i = tid;
      while(yp[i]) {
        yp[i] = tid;
        i++;
        if(i > MAX_Y) break;
      }
    }

    if(tid < MAX_X) {
      if(xp[tid] && !xp[tid - 1]) {
        uint i = tid;
        while(xp[i]) {
          xp[i] = tid;
          i++;
          if(i > MAX_X) break;
        }
      }
    }//end of if(tid<MAX_X){
  } //end of if(tid)
 __syncthreads();
  //assign the cluster id to each pixel
  for(uint i=0; i<noItr; i++) {
    uint idx = moduleBegin + tid + i*blockDim.x;
    if(idx <= moduleEnd) {
      uint px = xp[xx[idx]];                    // find location of pixel on x strip
      uint py = yp[yy[idx]];                    // find location of strip on y strip
      if(px!=0 && py!=0)
        gClusterId[idx] = 1000000 * (moduleId) + 1000 * py + px;     // gclusterId(idx) = f(px, py);
    } // if(idx <= moduleEnd)
  } // End of for(uint i=0; i<noItr...  

} //End of cluster_kernel

__global__ void init_kernel(const uint wordCounter, const uint total_cluster, uint *Index) {
  Index[total_cluster] = wordCounter;
}

/************* origin of call to the kernel***************/

void PixelCluster_Wrapper(uint *d_xx, uint *d_yy, uint *d_ADC,const uint wordCounter,
                         const int *mIndexStart,const int *mIndexEnd) 
{
    checkCUDAError("Error in RawToDigi, didn't enter in Cluster");    
    
    cout<<"Clustering started on GPU!"<<endl;
    //to measure the time
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);

    cudaMemset(d_gClusterId, 0, wordCounter*sizeof(uint));
    checkCUDAError("Error in setting memory to 0");

    // launch clustering kernel
    dim3 blockSize  = MAX_MODULE_SIZE;   // no of blocks
    dim3 threadsize = NO_THREADS;        // no of threads
    
    cluster_kernel <<< blockSize, threadsize>>>(d_xx, d_yy, mIndexStart, mIndexEnd, d_gClusterId);
    cudaDeviceSynchronize();
    checkCUDAError(" Failed after main kernel call");
    
    //cout<<"Finding Sub-clusters..."<<endl;
    sub_cluster(d_xx, d_yy, d_ADC, Index, d_gClusterId, wordCounter, d_gClusterId1, d_xx1, d_yy1, d_ADC1);
    
    // FOR CPE
    // formate the output of cluster before giving to CPE
    // sort the clusterIds and corrseponding attributes
    // to get the start and end index of cluster
    createIndex(wordCounter, Index);
    thrust::device_ptr<uint> ClusterId_ptr(d_gClusterId1);
    thrust::device_ptr<uint> Index_ptr(Index);
    
    thrust::sort_by_key(ClusterId_ptr , ClusterId_ptr + wordCounter, Index_ptr);
    cudaMemcpy(d_gClusterId, d_gClusterId1, wordCounter * sizeof(uint),  cudaMemcpyDeviceToDevice );
    // now sort the xx yy and ADC
    uint N_threads = 1024;
    uint N_blocks = wordCounter/N_threads +1;
    copy_kernel<<<N_blocks,N_threads>>>(Index, d_xx1, d_yy1,d_ADC1, wordCounter,d_xx,d_yy,d_ADC); 
    cudaDeviceSynchronize();
    checkCUDAError("Error in sorting ");

    uint total_cluster=0;
    thrust::pair<thrust::device_ptr<uint>, thrust::device_ptr<uint> > new_end;
  
    // Fill the index again which will be used for uniuqe 
    createIndex(wordCounter, Index);
    new_end = thrust::unique_by_key(ClusterId_ptr , ClusterId_ptr + wordCounter, Index_ptr );
    total_cluster = new_end.first - ClusterId_ptr;
    checkCUDAError(" Failed at Clustering");
    //cout<<"Total_clusters after sub-clustering: "<<total_cluster<<endl;
    //cout<<"Clustering finished on GPU!"<<endl;
    // call CPE function
    init_kernel<<<1,1>>>(wordCounter, total_cluster, Index);
    //Index[total_cluster] = wordCounter;
    //since origin is shifted by (1,1) move it back to (0,0) before giving it CPE
    shift_origin_kernel<<<N_blocks, N_threads>>>(wordCounter,d_xx,d_yy); 
    cudaDeviceSynchronize();
    
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    float microseconds = milliseconds*1000;

    ofstream GPUTime("Cluster_GPU_Time.txt",ios::out | ios::app);
    static int eventNo = 1;
    cout<<"Event "<<eventNo<<endl;
    if (eventNo==1) {
      GPUTime<<"Event#\t  "<<"Total pixel\t  "<<"Time(us)"<<endl;
    }
    GPUTime<<setw(4)<<eventNo<<setw(12)<<wordCounter<<setw(14)<<microseconds<<endl;
    
    GPUTime.close();
    uint *xx, *yy, *adc_h, *gClusterId, *index_h;
    xx = (uint*)malloc(wordCounter*sizeof(uint));
    yy = (uint*)malloc(wordCounter*sizeof(uint));
    adc_h = (uint*)malloc(wordCounter*sizeof(uint));
    gClusterId = (uint*)malloc(wordCounter*sizeof(uint));
    index_h = (uint*)malloc(wordCounter*sizeof(uint));
    int count = wordCounter*sizeof(uint);
    cudaMemcpy(xx, d_xx,count , cudaMemcpyDeviceToHost);
    cudaMemcpy(yy, d_yy,count , cudaMemcpyDeviceToHost);
    cudaMemcpy(gClusterId, d_gClusterId1, count , cudaMemcpyDeviceToHost);
    cudaMemcpy(adc_h, d_ADC, count, cudaMemcpyDeviceToHost);
    cudaMemcpy(index_h, Index, count, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    ofstream cfile;
    cfile.open("CPE_Input_CPU_PartA"+to_string(eventNo)+".txt");
    cfile<<"Index\tClusterId"<<endl;
    for (int i = 0; i < total_cluster; i++) {
      //cfile<<setw(10)<<gClusterId[i]<<setw(6)<<xx[i]<<setw(6)<<yy[i]<<setw(10)<<adc_h[i]<<endl;
      cfile<<setw(6)<<index_h[i]<<setw(14)<<gClusterId[i]<<endl;
    }
    cfile.close();
    cudaMemcpy(gClusterId, d_gClusterId, count , cudaMemcpyDeviceToHost);
    cfile.open("CPE_Input_CPU_PartB"+to_string(eventNo)+".txt");
    cfile<<"ClusterId\t\txx\tyy\tADC"<<endl;
    for (int i = 0; i < wordCounter; i++) {
      //cfile<<setw(10)<<gClusterId[i]<<setw(6)<<xx[i]<<setw(6)<<yy[i]<<setw(10)<<adc_h[i]<<endl;
      cfile<<setw(6)<<i<<setw(14)<<gClusterId[i]<<setw(6)<<xx[i]<<setw(6)<<yy[i]<<setw(10)<<adc_h[i]<<endl;
    }
    cfile.close();
    eventNo++;
    free(xx);
    free(yy);
    free(adc_h);
    free(gClusterId);
    free(index_h);
    CPE_wrapper(total_cluster,d_gClusterId1, Index, d_xx, d_yy, d_ADC);

  
    // move it into destructor
    
} //end of pixel clusterizer

void initDeviceMemCluster() {
    const int MAX_FED = 150; // not all are present typically 108
    const int MAX_WORD = 2000; // don't know the exact max word, for PU70 max was 2900
    const int size = MAX_FED*MAX_WORD*sizeof(uint);
    cudaMalloc((void**)&Index , size*sizeof(uint));
    checkCUDAError("Error in unified memory allocation");
    cudaMalloc((void**)&d_xx1, size*sizeof(uint));
    cudaMalloc((void**)&d_yy1, size*sizeof(uint));
    cudaMalloc((void**)&d_gClusterId, size*sizeof(uint));
    cudaMalloc((void**)&d_gClusterId1, size*sizeof(uint));
    cudaMalloc((void**)&d_ADC1,  size*sizeof(uint));
    checkCUDAError("Error in memory allocation for clustering");
}

void freeDeviceMemCluster() {
    cudaFree(Index);
    cudaFree(d_xx1);
    cudaFree(d_yy1);
    cudaFree(d_gClusterId);
    cudaFree(d_gClusterId1);
    cudaFree(d_ADC1);
}