/*Copyright 2017 Sushil Dubey, Shashi Dugad, TIFR

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/


#ifndef RAWTODIGI_CPU_GPU_H
#define RAWTODIGI_CPU_GPU_H

// wrapper function to call RawToDigi on the GPUvfrom host side
void RawToDigi_Cluster_CPE_wrapper (const unsigned int wordCounterGPU,
     unsigned int *word,const unsigned int fedCounter,unsigned int *fedIndex, unsigned int *eventIndex);

void initCablingMap();
void initDeviceMemory();
void freeMemory();

// reference cmssw/RecoLocalTracker/SiPixelClusterizer
// all are runtime const, should be specified in python _cfg.py
struct ADCThreshold {
  const int thePixelThreshold = 1000; // default Pixel threshold in electrons
  const int theSeedThreshold = 1000; //seed thershold in electrons not used in our algo
  const float theClusterThreshold = 4000; // Cluster threshold in electron
  const int ConversionFactor = 65;  // adc to electron conversion factor
  
  // following are the default value
  // not present in _cfg.py
  const int theStackADC_  = 255; // the maximum adc count for stack layer
  const int theFirstStack_ = 5; // the index of the fits stack layer
  const double theElectronPerADCGain_ = 135; //ADC to electron conversion

};
#endif
