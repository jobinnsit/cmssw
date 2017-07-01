#ifndef PIXEL_CLUSTER_GPU_H
#define PIXEL_CLUSTER_GPU_H

typedef unsigned int uint;
// wrapper function to CPE kernel
void CPE_wrapper(const uint total_cluster, const uint64 *ClusterId, const uint *Index, const uint *xx, const uint *yy,
                 const uint *adc );
             
//const uint MAX_EVENT_SIZE  = 200000;   // Number of pixel hit in one events
const uint MAX_MODULE_SIZE = 1856;     // Number of module in one event
const uint NO_THREADS      = 416;      // x:0-159 y:0-415
const uint MAX_X           = 160;      // origin is shifted by 1
const uint MAX_Y           = 416; 

//auxilary memory required for the clustering
// since our sorting is not in-place sort.
uint *Index, *d_xx1, *d_yy1, *d_ADC1; 
uint64 *d_gClusterId, *d_gClusterId1; //uint64
#endif