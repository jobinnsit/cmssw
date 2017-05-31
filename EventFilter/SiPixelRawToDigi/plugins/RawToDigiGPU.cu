/**2017-03-02  Sushil Dubey  <sdubey@felk40.cern.ch>
 *
 * File Name: RawToDigiGPU.cu
 * Description: It converts Raw data into Digi data using GPU
 * then it applies the adc threshold to drop the dead pixels
 * The Output of RawToDigi data is given to pixelClusterizer
 *
**/ 
// System includes
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <assert.h>
#include <iomanip>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/execution_policy.h>
// CUDA runtime
#include <cuda.h>
#include <cuda_runtime.h>
#include "RawToDigiGPU.h"
#include "RawToDigiCPUGPU.h"
using namespace std;

// forward declaration to be moved in header file
void PixelCluster_Wrapper(uint *xx_adc, uint *yy_adc, uint *adc_d,const uint wordCounter, 
                          const int *mIndexStart, const int *mIndexEnd, uint *xx, uint *yy);

/*
void initCablingMap() {

  ifstream mapFile;
  mapFile.open("RawId_ModuleId_CablingMap_ArrayFile.txt");
  string str;
  getline(mapFile, str);
  uint rawId, moduleId, rocInDU;
  int i =1;  // cabling map index starts at 1
  while(!mapFile.eof()) {
    mapFile >> rawId >> rocInDU >> moduleId;
    Map->RawId[i]    = rawId;
    Map->rocInDet[i] = rocInDU;
    Map->moduleId[i] = moduleId;
    i++;
  }
  mapFile.close();
}
*/
// New cabling Map
void initCablingMap() {

  ifstream mapFile;
  mapFile.open("Pixel_Phase1_Raw2Digi_GPU_Cabling_Map.dat");
  string str;
  getline(mapFile, str);
  uint Index, FedId, Link, idinLNK, B_F, RawID, idinDU, ModuleID;
  int i =1;  // cabling map index starts at 1
  while(!mapFile.eof()) {
    mapFile >> Index>>FedId>>Link>>idinLNK>>B_F>>RawID>>idinDU>>ModuleID;
    Map->RawId[i] = RawID;
    Map->rocInDet[i] = idinDU;
    Map->moduleId[i] = ModuleID;
    i++;
  }
  mapFile.close();
  cout<<"Cabling Map uploaded successfully!"<<endl;
}

void initDeviceMemory() {
  int sizeByte = MAX_FED * MAX_LINK * MAX_ROC * sizeof(uint)+sizeof(uint);
  // Unified memory for cabling map
  cudaMallocManaged((void**)&Map,  sizeof(CablingMap));
  cudaMallocManaged((void**)&Map->RawId,    sizeByte);
  cudaMallocManaged((void**)&Map->rocInDet, sizeByte);
  cudaMallocManaged((void**)&Map->moduleId, sizeByte);
    // Number of words for all the feds 
  const uint MAX_WORD_SIZE = MAX_FED*MAX_WORD*sizeof(uint); 
  
  // all the below host memory are only for testing purpose
  // host memory only be used in CPE
  adc_h  = (uint*)malloc(MAX_WORD_SIZE);
  word_h = (uint*)malloc(MAX_WORD_SIZE);
  fedIndex_h = (uint*)malloc((MAX_FED+1)*sizeof(uint)); // +1 for last fed index
  xx     =   (uint*)malloc(MAX_WORD_SIZE);
  yy     =   (uint*)malloc(MAX_WORD_SIZE);
  RawId  =   (uint*)malloc(MAX_WORD_SIZE);
  moduleId = (uint*)malloc(MAX_WORD_SIZE);
  int mSize = totalModule*sizeof(int);
  mIndexStart = (int*)malloc(mSize); 
  mIndexEnd = (int*)malloc(mSize);
 
  cudaMalloc((void**)&word_d,       MAX_WORD_SIZE);
  cudaMalloc((void**)&fedIndex_d,   (MAX_FED+1)*sizeof(uint));
  cudaMalloc((void**)&xx_d,         MAX_WORD_SIZE); // to store the x and y coordinate
  cudaMalloc((void**)&yy_d,         MAX_WORD_SIZE);
  cudaMalloc((void**)&xx_adc,         MAX_WORD_SIZE); // to store the x and y coordinate
  cudaMalloc((void**)&yy_adc,         MAX_WORD_SIZE);
  cudaMalloc((void**)&adc_d,        MAX_WORD_SIZE);
  cudaMalloc((void**)&layer_d ,     MAX_WORD_SIZE);
  cudaMalloc((void**)&RawId_d,      MAX_WORD_SIZE);
  cudaMalloc((void**)&moduleId_d,   MAX_WORD_SIZE);
  cudaMalloc((void**)&mIndexStart_d, mSize);
  cudaMalloc((void**)&mIndexEnd_d, mSize);
  
  cout<<"Memory Allocated successfully !\n";
  // Upload the cabling Map
  initCablingMap();
  
}

void freeMemory() {

  //GPU specific
  // memory used for testing purpose will be released during
  // deployment.
  free(word_h);
  free(fedIndex_h);
  free(adc_h);
  free(xx);
  free(yy);
  free(RawId);
  free(moduleId);
  free(mIndexStart);
  free(mIndexEnd);
  cudaFree(word_d);
  cudaFree(fedIndex_d);
  cudaFree(adc_d);
  cudaFree(layer_d);
  cudaFree(xx_d);
  cudaFree(yy_d);
  cudaFree(xx_adc);
  cudaFree(yy_adc);
  cudaFree(RawId_d);
  cudaFree(moduleId_d);
  cudaFree(mIndexStart_d);
  cudaFree(mIndexEnd_d);
  cudaFree(Map->RawId);
  cudaFree(Map->rocInDet); 
  cudaFree(Map->moduleId);
  cudaFree(Map);
  cout<<"Memory Released !\n";

}

__device__ uint getLink(uint ww)  {
  //printf("Link_shift: %d  LINK_mask: %d\n", LINK_shift, LINK_mask);
  return ((ww >> LINK_shift) & LINK_mask);
}

__device__ uint getRoc(uint ww) {
  return ((ww >> ROC_shift ) & ROC_mask);
}
__device__ uint getADC(uint ww) {
  return ((ww >> ADC_shift) & ADC_mask);
}

__device__ bool isBarrel(uint rawId) {
  return (1==((rawId>>25)&0x7));
}
//__device__ uint FED_START = 1200;

__device__ DetIdGPU getRawId(const CablingMap *Map, uint fed, uint link, uint roc) {
  uint index = fed * MAX_LINK* MAX_ROC + (link-1)* MAX_ROC + roc;
  DetIdGPU detId = {Map->RawId[index], Map->rocInDet[index], Map->moduleId[index]};
  return detId;  
}
// Convert local pixel to global pixel
__device__ Pixel frameConversion(bool bpix, int side, uint rocIdInDetUnit, Pixel local) {
  
  int slopeRow  = 0,  slopeCol = 0;
  int rowOffset = 0, colOffset = 0;

  if(bpix) {
    if(side==1) {
      if(rocIdInDetUnit <8) {
        slopeRow  = -1;
        slopeCol  =  1;
        rowOffset = 2*numRowsInRoc-1;
        colOffset = rocIdInDetUnit * numColsInRoc;
      }
      else {
        slopeRow  = 1;
        slopeCol  = -1;
        rowOffset = 0;
        colOffset = (16-rocIdInDetUnit)*numColsInRoc-1;
      }
    }
    else {
      if (rocIdInDetUnit <8) {
        slopeRow = 1;     slopeCol = -1;
        rowOffset = 0;
        colOffset = (8-rocIdInDetUnit)*numColsInRoc-1;
      }
      else {
        slopeRow  = -1;
        slopeCol  = 1;
        rowOffset = 2*numRowsInRoc-1;
        colOffset = (rocIdInDetUnit-8)*numColsInRoc;
      } // if roc
    }
  }
  else { // fpix
    if(side==-1) { // pannel 1
      if (rocIdInDetUnit < 8) {
        slopeRow = 1;
        slopeCol = -1;
        rowOffset = 0;
        colOffset = (8-rocIdInDetUnit)*numColsInRoc-1;
      }
      else {
        slopeRow = -1;
        slopeCol = 1;
        rowOffset = 2*numRowsInRoc-1;
        colOffset = (rocIdInDetUnit-8)*numColsInRoc;
      }
    }
    else { // pannel 2
      if (rocIdInDetUnit < 8) {
        slopeRow = 1;
        slopeCol = -1;
        rowOffset = 0;
        colOffset = (8-rocIdInDetUnit)*numColsInRoc-1;
      }
      else {
        slopeRow = -1;
        slopeCol = 1;
        rowOffset = 2*numRowsInRoc-1;
        colOffset = (rocIdInDetUnit-8)*numColsInRoc;
      }

    } // side

  }

  uint gRow = rowOffset+slopeRow*local.row;
  uint gCol = colOffset+slopeCol*local.col;
  //printf("Inside frameConversion gRow: %u  gCol: %u\n",gRow, gCol);
  Pixel global = {gRow, gCol};
  return global;
}


/*----------
* Name: applyADCthreshold_kernel()
* Desc: converts adc count to electrons and then applies the 
* threshold on each channel. 
* make pixel to 0 if it is below the threshold
* Input: xx_d[], yy_d[], layer_d[], wordCounter, adc[], ADCThreshold
*-----------
* Output: xx_adc[], yy_adc[] with pixel threshold applied 
*/
// kernel to apply adc threshold on the channels	
__global__ void applyADCthreshold_kernel
(const uint *xx_d, const uint *yy_d, const uint *layer_d, uint *adc, const uint wordCounter,
 const ADCThreshold adcThreshold, uint *xx_adc, uint *yy_adc ) {
  int tid = threadIdx.x;
  int gIndex = blockDim.x*blockIdx.x+tid;
  if(gIndex<wordCounter) {
    //int i=0;
    //for(DigiIterator di = begin; di != end; ++di) {
      uint adcOld = adc[gIndex];
      const float gain = adcThreshold.theElectronPerADCGain_; // default: 1 ADC = 135 electrons
      const float pedestal = 0; //
      int adcNew = int(adcOld*gain+pedestal);
      // rare chance of entering into the if()
      if (layer_d[gIndex]>=adcThreshold.theFirstStack_) {
        if (adcThreshold.theStackADC_==1 && adcOld==1) {
          adcNew = int(255*135); // Arbitrarily use overflow value.
        }
        if (adcThreshold.theStackADC_ >1 && adcThreshold.theStackADC_!=255 && adcOld>=1){
          adcNew = int((adcOld-1) * gain * 255/float(adcThreshold.theStackADC_-1));
        }
      }
  
    if(adcNew >adcThreshold.thePixelThreshold ) {
      xx_adc[gIndex]=xx_d[gIndex];
      yy_adc[gIndex]=yy_d[gIndex];
    }
    else {
      xx_adc[gIndex]=0; // 0: dead pixel
      yy_adc[gIndex]=0;
    }
    adc[gIndex] = adcNew;
  }
}  


// Kernel to perform Raw to Digi conversion
__global__ void RawToDigi_kernel(const CablingMap *Map,const uint *Word,const uint *fedIndex, 
                                 uint *XX, uint *YY, uint *RawId, uint *moduleId, 
                                 int *mIndexStart, int *mIndexEnd, uint *ADC, uint *layerArr ) 
{
  //printf("Inside GPU: \n");
  int fedId    = blockIdx.x;
  int threadId = threadIdx.x;

  int begin  = fedIndex[fedId];
  int end    = fedIndex[fedId+1];
 
  int no_itr = (end - begin)/ blockDim.x + 1; // to deal with number of hits greater than blockDim.x 
  #pragma unroll
  for(int i =0; i<no_itr; i++) { // use a static number to optimize this loop
    int gIndex = begin + threadId + i*blockDim.x;  // *optimize this
    if(gIndex <end) {
      uint ww    = Word[gIndex]; // Array containing 32 bit raw data
      if(ww == 0 ) {
        //noise and dead channels are ignored
        XX[gIndex] = 0;  // 0 is an indicator of a noise/dead channel
        YY[gIndex]  = 0; // skip these pixels during clusterization
        RawId[gIndex] = 0; 
        ADC[gIndex]   = 0; 
        moduleId[gIndex] = 9999; //9999 is the indication of bad module, taken care later  
        layerArr[gIndex] = 0;
        //fedIdArr[gIndex] = fedId; // used for testing
        continue ;         // 0: bad word, 
      } 
      uint link  = getLink(ww);            // Extract link
      uint roc   = getRoc(ww);             // Extract Roc in link
      DetIdGPU detId = getRawId(Map, fedId, link, roc);
      uint rawId  = detId.RawId;
      uint rocIdInDetUnit = detId.rocInDet;

      bool barrel = isBarrel(rawId);
  
      //printf("ww: %u    link:  %u  roc: %u   rawId: %u\n", ww, link, roc, rawId);
      //printf("from CablingMap  rocInDU: %u  moduleId: %u", rocIdInDetUnit, detId.moduleId);
      //printf("barrel: %d\n", barrel);
      uint layer =0, ladder =0;
      int side =0, panel =0, disk =0, blade =0, module=0;
    
      if(barrel) {
        layer  = (rawId >> layerStartBit_)  & layerMask_;
        ladder = (rawId >> ladderStartBit_) & ladderMask_;
        module = (rawId >> moduleStartBit_) & moduleMask_;
        side   = (module<5)? -1:1;
     
      }
      else {
        // endcap ids
        layer = 0;
        panel = (rawId >> panelStartBit_) & panelMask_;
        disk  = (rawId >> diskStartBit_)  & diskMask_ ;
        side  = (panel==1)? -1:1;
        //blade = (rawId>>bladeStartBit_) & bladeMask_;
      }
      // ***special case of layer to 1 be handled here
      Pixel localPix;
      if(layer==1) {
        uint col = (ww >> COL_shift) & COL_mask;
        uint row = (ww >> ROW_shift) & ROW_mask;
        localPix.row = row;
        localPix.col = col;
        //if(event==0 && fedId==0)
         //printf("col: %u  row: %u\n",col, row);
      }
      else {
        // ***conversion rules for dcol and pxid
        uint dcol = (ww >> DCOL_shift) & DCOL_mask;
        uint pxid = (ww >> PXID_shift) & PXID_mask;
        uint row  = numRowsInRoc - pxid/2;
        uint col  = dcol*2 + pxid%2;
        localPix.row = row;
        localPix.col = col;
      }

      Pixel globalPix = frameConversion(barrel, side, rocIdInDetUnit, localPix);
      XX[gIndex]    = globalPix.row +1 ; // origin shifting by 1 0-159
      YY[gIndex]    = globalPix.col +1 ; // origin shifting by 1 0-415
      ADC[gIndex]   = getADC(ww);
      RawId[gIndex] = detId.RawId; // only for testing
      layerArr[gIndex] = layer;
      //fedIdArr[gIndex] = fedId;     // used for testing
      // only for testing purpose: White box testing
      //XX[gIndex] = fedId; // fedId for pattern
      //YY[gIndex] = gIndex; // wwIndex for pattern
      //RawId[gIndex] = ww;
      moduleId[gIndex] = detId.moduleId;
    } // end of if(gIndex < end)
  } // end of for(int i =0;i<no_itr...)
  
  __syncthreads();
  // three cases possible
  // case 1: 21 21 21 22 21 22 22
  // pos   : 0  1  2  3  4  5  6
  // solution swap 21 with 22 : 21 21 21 21 22 22 22
  // atomicExch(address, value), set the variable at address to value.
  // do the swapping for above case and replace the 9999 with 
  // valid moduleId
  for(int i =0; i<no_itr; i++) { 
    int gIndex = begin + threadId + i*blockDim.x;  
    if(gIndex <end) {
      //rare condition 
      if(moduleId[gIndex]==moduleId[gIndex+2] && moduleId[gIndex]<moduleId[gIndex+1]) {
        atomicExch(&moduleId[gIndex+2], atomicExch(&moduleId[gIndex+1], moduleId[gIndex+2]));
        //*swap all the digi id
        atomicExch(&XX[gIndex+2], atomicExch(&XX[gIndex+1], XX[gIndex+2]));
        atomicExch(&YY[gIndex+2], atomicExch(&YY[gIndex+1], YY[gIndex+2]));
        atomicExch(&RawId[gIndex+2], atomicExch(&RawId[gIndex+1], RawId[gIndex+2])); 
        //atomicExch(&fedIdArr[gIndex+2], atomicExch(&fedIdArr[gIndex+1], fedIdArr[gIndex+2]));
        atomicExch(&ADC[gIndex+2], atomicExch(&ADC[gIndex+1], ADC[gIndex+2]));
        atomicExch(&layerArr[gIndex+2], atomicExch(&layerArr[gIndex+1], layerArr[gIndex+2]));
      }
      __syncthreads();
      //rarest condition
      // above condition fails at 361 361 361 363 362 363 363
      // here we need to swap 362 with previous 363
      if(moduleId[gIndex]==moduleId[gIndex+2] && moduleId[gIndex]>moduleId[gIndex+1]) {
        atomicExch(&moduleId[gIndex+1], atomicExch(&moduleId[gIndex], moduleId[gIndex+1]));
        //*swap all the digi id
        atomicExch(&XX[gIndex+1], atomicExch(&XX[gIndex], XX[gIndex+1]));
        atomicExch(&YY[gIndex+1], atomicExch(&YY[gIndex], YY[gIndex+1]));
        atomicExch(&RawId[gIndex+1], atomicExch(&RawId[gIndex], RawId[gIndex+1])); 
        //atomicExch(&fedIdArr[gIndex+2], atomicExch(&fedIdArr[gIndex+1], fedIdArr[gIndex+2]));
        atomicExch(&ADC[gIndex+1], atomicExch(&ADC[gIndex], ADC[gIndex+1]));
        atomicExch(&layerArr[gIndex+1], atomicExch(&layerArr[gIndex], layerArr[gIndex+1]));
      }
      // moduleId== 9999 then pixel is bad with x=y=layer=adc=0
      // this bad pixel will not affect the cluster, since for cluster
      // the origin is shifted at (1,1) so x=y=0 will be ignored
      // assign the previous valid moduleId to this pixel to remove 9999
      // so that we can get the start & end index of module easily.
      __syncthreads(); // let the swapping finish first
      if(moduleId[gIndex]==9999) {
        int m=gIndex;
        while(moduleId[--m]==9999) {} //skip till you get the valid module
        moduleId[gIndex]=moduleId[m];
      } 
    } // end of if(gIndex<end)
  } //  end of for(int i=0;i<no_itr;...)
  __syncthreads();

  // mIndexStart stores staring index of module 
  // mIndexEnd stores end index of module 
  // both indexes are inclusive 
  // check consecutive module numbers
  // for start of fed
  for(int i =0; i<no_itr; i++) { 
    int gIndex = begin + threadId + i*blockDim.x;  
    if(gIndex <end) {
      if(gIndex == begin) {
        mIndexStart[moduleId[gIndex]] = gIndex;
      }
      // for end of the fed
      if(gIndex == (end-1)) {  
        mIndexEnd[moduleId[gIndex]] = gIndex;
      }   
      // point to the gIndex where two consecutive moduleId varies
      if(gIndex!= begin && (gIndex<(end-1)) && moduleId[gIndex]!=9999) {
        if(moduleId[gIndex]<moduleId[gIndex+1] ) {
          mIndexEnd[moduleId[gIndex]] = gIndex;
        }
        if(moduleId[gIndex] > moduleId[gIndex-1] ) {
          mIndexStart[moduleId[gIndex]] = gIndex;
        } 
      } //end of if(gIndex!= begin && (gIndex<(end-1)) ...  
    } //end of if(gIndex <end) 
  }
} // end of Raw to Digi kernel

// kernel wrapper called from runRawToDigi_kernel
void RawToDigi_kernel_wrapper(const uint wordCounter,uint *word, uint *fedIndex) { 
  
 
  cout<<"Inside RawToDigi , total words: "<<wordCounter<<endl;
  int nBlocks = MAX_FED; // = MAX_FED
  int threads = 512; //
  fedIndex[nBlocks] = wordCounter;
  // for debugging 
  uint *fedId;
  int mSize = totalModule*sizeof(int);
  { 
    uint eventSize = wordCounter*sizeof(uint);
	  // initialize moduleStart & moduleEnd with some constant(-1)
	  // number just to check if it updated in kernel or not
    cudaMemset(mIndexStart_d, -1, mSize);
    cudaMemset(mIndexEnd_d, -1, mSize);
    cudaMemcpy(word_d, word, eventSize, cudaMemcpyHostToDevice);
    cudaMemcpy(fedIndex_d, fedIndex, (MAX_FED+1)*sizeof(uint), cudaMemcpyHostToDevice); 
    // for debugging 
    cudaMallocManaged((void**)&fedId, eventSize);
    // Launch rawToDigi kernel
    RawToDigi_kernel<<<nBlocks,threads>>>(Map,word_d, fedIndex_d, xx_d, yy_d, RawId_d,
                                          moduleId_d, mIndexStart_d, mIndexEnd_d, adc_d, layer_d);
    cudaDeviceSynchronize();

    cudaMemcpy(yy  , yy_d,    eventSize,    cudaMemcpyDeviceToHost);
    cudaMemcpy(adc_h, adc_d, eventSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(moduleId, moduleId_d, eventSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(xx,    xx_d,   eventSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(mIndexStart, mIndexStart_d, mSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(mIndexEnd,   mIndexEnd_d,   mSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(word, layer_d, eventSize, cudaMemcpyDeviceToHost);

     // apply the correction to the moduleStart & moduleEnd
     // if module contains only one pixel then either moduleStart 
     // or moduleEnd is not updated(remains 9999) in RawToDigi kernel
     // ex. moduleStart[1170] =9999 & moduleEnd[1170] = 34700
     // because of 1 pixel moduleStart[1170] didn't update
     // as per the if condition
    
    for(int i=0;i<totalModule;i++) {
      // if module is empty then index are not updated in kernel
      if(mIndexStart[i]==-1 && mIndexEnd[i]==-1) {
        mIndexStart[i]=0;
        mIndexEnd[i]=0;
      }
      else if(mIndexStart[i]==-1) {
        mIndexStart[i] = mIndexEnd[i];
      }
      else if(mIndexEnd[i]==-1) {
        mIndexEnd[i] = mIndexStart[i];
      }
    }
    //copy te data back to the device memory
    cudaMemcpy(mIndexStart_d, mIndexStart, mSize, cudaMemcpyHostToDevice);
    cudaMemcpy(mIndexEnd_d,   mIndexEnd,   mSize, cudaMemcpyHostToDevice);
    
  }
  //static int eventno = 0;
  //ofstream outFile;
  //outFile.open("fedId_moduleId_afterBugFix.txt", ios::out | ios::app);  
  //for(uint i=0; i<wordCounter;i++) {
    //if(RawId[i]!=0)
    //outFile <<setw(6)<<moduleId[i]+1200<<setw(14)<<word[i]<<setw(6)<<xx[i]<<setw(6)<<yy[i]<<setw(6)<<adc_h[i]<<endl;
    //outFile<<setw(4)<<word[i]<<setw(8)<<moduleId[i]<<setw(4)<<xx[i]<<setw(4)<<yy[i]<<endl;
     //cout<<"ww: "<<setw(10)<<RawId[i]<<"  xx: "<<setw(3)<<xx[i]<<"  yy: "<<setw(3)<<yy[i]<<endl;
  //}
  //outFile.close();
  cudaFree(fedId);
  //cout<<"RawToDigi Kernel executed successfully!\n";
  
  /*ofstream outFile; 
  outFile.open("InputForCluster.txt");
  outFile<<"FedId   "<<"  wwIndex   "<<"  moduleId  "<<" RawId  "<<endl;
  for(uint i=0; i<wordCounter;i++) {
   // if(RawId[i]!=0)
    outFile <<setw(10)<<xx[i]<<"\t\t"<<setw(10)<<yy[i]<<endl;
    //outFile<<setw(4)<<xx[i]<<setw(12)<<yy[i]<<setw(12)<<moduleId[i]<<setw(16)<<RawId[i]<<endl;
    //cout<<"ww: "<<setw(10)<<RawId[i]<<"  xx: "<<setw(3)<<xx[i]<<"  yy: "<<setw(3)<<yy[i]<<endl;
  }
  */
  //ofstream mse("ModuleStartEndIndex.txt");
  //for(int i=0;i<totalModule;i++) {
   // mse<<mIndexStart[i]<<"\t\t"<<mIndexEnd[i]<<endl;
    //cout<<mIndexStart[i]<<"\t\t"<<mIndexEnd[i]<<endl;
  //}
  //mse.close();
  //outFile.close();                                           
  //++eventno;
   
  //cout<<"Calling pixel cluster"<<endl;
  // kernel to apply adc threashold on the channel
  ADCThreshold adcThreshold;
  uint numThreads = 512;
  uint numBlocks = wordCounter/512 +1;
  applyADCthreshold_kernel<<<numBlocks, numThreads>>>(xx_d, yy_d,layer_d,adc_d,wordCounter,adcThreshold, xx_adc, yy_adc);
  cudaDeviceSynchronize();
  // only for testing purpose these memcpy used
  //cudaMemcpy(yy  , yy_adc,    wordCounter*sizeof(uint),    cudaMemcpyDeviceToHost);
  //cudaMemcpy(adc_h, adc_d, wordCounter*sizeof(uint), cudaMemcpyDeviceToHost);
  //cudaMemcpy(moRawToDigiOutput_after_adcThreshold.txtduleId, moduleId_d, eventSize, cudaMemcpyDeviceToHost);
  //cudaMemcpy(xx,    xx_adc, wordCounter*sizeof(uint), cudaMemcpyDeviceToHost);
  //for(uint i=0; i<wordCounter;i++) {
      //if(RawId[i]!=0)
      //cout <<"\t\t"<<adc_h[i]<<"\t\t"<<xx[i]<<"\t\t"<<yy[i]<<endl;
  //}
  // call to pixelClusterizer kernel from here
  //PixelCluster_Wrapper(xx_adc , yy_adc, adc_d,wordCounter, mIndexStart_d, mIndexEnd_d,xx,yy);
}
