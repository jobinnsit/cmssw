/*
 * CPU clustering code
 * Input: xx[] yy[], moduleStartIndex[], moduleEndIndex[]
 * Output: clusterId[] xx[] yy[]
 * finds the cluster
 * Note: this input is taken from GPU RawToDigi code
 * it genererate the clusterid in similar fashion to that of GPU cluster
 * clusterId <- f(eventno, moduleno, x_strip, y_strip)
 * */

#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <malloc.h>
#include <time.h>
#include <algorithm>
#include <iomanip>

#include "EventInfoGPU.h"
#include "PixelClusterUtil.h"

const int MAX_FED = 150;
const int MAX_WORD = 2000;

constexpr int SIZE = MAX_FED*MAX_WORD*NEVENT*sizeof(int);

using namespace std;
void  getClusters(uint event,uint module_no);
int   getSubRectangles(uint event, uint module_no,int itrn, uint64 clusterId, int &icount, int &ndivision, uint64 *new_rect_idx);

uint64 getClusterId(uint event, uint module, uint y, uint x) {
  uint64 clusterId = ((uint64)event << EVENT_shift) | (module << MODULE_shift) | (y << YCOR_shift) | (x << XCOR_shift);
  return clusterId;
}

                          
struct PixelCluster{
  //Max. # of allowed PixelHits in a Module=2000
  //Max. # of allowed clusters  in a Module=200
  int nPixHits;          //# of Pixel Hits
  int nClusters;         //# of Clusters
  int xx[2000];          //x-cordiante of each pixel in module
  int yy[2000];          //y-cordiante of each pixel in module
  uint64 rect_idxPt[2000];  //Cluster # associated with each pixel
  uint64 rect_idx[200];           //Cluster #
  int cluster_Xlocation[200];  //x-cordiante of bottom left of cluster
  int cluster_Ylocation[200];  //y-cordiante of bottom left of cluster
  int cluster_Xwidth[200];     //x-Width of cluster
  int cluster_Ywidth[200];     //y-Width of cluster
  int cluster_size[200];       //Number Pixel Hits in cluster
};
struct PixelCluster gCluster;

int main(int argc, char **argv) {
  ofstream outputFile, CPUTIME;
  outputFile.open("Cluster_CPU.txt");
  outputFile <<"ClusterId\t\t"<<" xx\t "<<" yy"<<endl;
  
  int *xx, *yy, mIndexStart[NMODULE*NEVENT+1], mIndexEnd[NMODULE*NEVENT+1];
  xx = (int*)malloc(SIZE);
  yy = (int*)malloc(SIZE);
  int module_size =0 , moduleBegin=0, moduleEnd=0;

  ifstream ifile;
  ifile.open("Cluster_Input_CPU.txt");
  string str;
  getline(ifile, str);
  uint x, y, i=0;
  while(!ifile.eof()) {
    ifile>> x>> y;
    xx[i] =x;
    yy[i] =y; 
    i++;
  }
  ifile.close();
  uint wordCounter = i;
  i=0;
  uint mId, start , end;
  ifile.open("ModuleStartEndIndex.txt");
  getline(ifile, str);
  while(!ifile.eof()) {
    ifile>>mId>>start>>end;
    mIndexEnd[i] = end;
    mIndexStart[i] = start;
      //cout<<"i: "<<i<<" end: "<<end<<" start: "<<start<<endl;
    i++;
  }
  ifile.close();
  cout<<"Total Event: "<<NEVENT<<"  Total hits: "<<wordCounter<<" Total Module: "<<i-1<<endl;
  for(uint Event=0;Event<NEVENT;Event++) {
    for(int ml=0;ml<NMODULE;ml++) {          // for each module form the cluster
      // intialization
      for(int i=0; i<2000; i++){
        gCluster.xx[i]=0; gCluster.yy[i]=0;
        gCluster.rect_idxPt[i]=0;
      }

      for(int i=0; i<200; i++) {
        gCluster.rect_idx[i]=0;
      }
      // feed the data for each module
      moduleBegin = mIndexStart[NMODULE*Event+ml];
      moduleEnd   = mIndexEnd[NMODULE*Event+ml];
      //skip empty and bad module
	    if(moduleBegin==0 && moduleEnd==0) {continue;}
      // if module contains only on pixel
   	  if( moduleEnd==moduleBegin ) {
        uint64 gClusterId;
        gClusterId = getClusterId(Event, ml,yy[moduleBegin], xx[moduleBegin]);
        outputFile <<setw(20)<<gClusterId<<setw(6)<<xx[moduleBegin]<<setw(6)<<yy[moduleBegin]<<endl;
        continue;
	    }
    
      module_size = moduleEnd-moduleBegin+1;
      for(int j=0;j<module_size;j++) {
        gCluster.xx[j] = xx[moduleBegin+j]; 
        gCluster.yy[j] = yy[moduleBegin+j];
      }

      gCluster.nPixHits=module_size;
      // call cluster function to find all the cluster in the module
      getClusters(Event, ml);
      // store the output
      int pix_no =0;
      for (int i=0;i<module_size;i++) {
        outputFile <<setw(20)<<gCluster.rect_idxPt[i]<<setw(6)<<gCluster.xx[i]<<setw(6)<<gCluster.yy[i]<<endl;
      }
    } // end of for(int ml=0;ml<module.size();ml++)
  } // end of event
  free(xx);
  free(yy);
  return 0;
}

void getClusters(uint event, uint module_no) {
  int  nrectangles, nPoints;
  uint64 new_rect_idx[200];

  //Initialise local rectangles (clusters)
  for(int i=0; i<200; i++){new_rect_idx[i]=0;}       

  nPoints = gCluster.nPixHits;

  //Cluster # is ONE for all hits in the first iteration
  for(int i=0; i<nPoints; i++){gCluster.rect_idxPt[i]=1;}
  //To begin with, there is only ONE cluster that includes entire frame
  gCluster.rect_idx[0]=1;
  gCluster.nClusters = 1;

  for(int itrn=0; itrn<2; itrn++){
    nrectangles=gCluster.nClusters;
    if(nrectangles>999) nrectangles=999;
    int inew_count=0; int ndivision=0;
    //New IDX to be copied to the Global IDX later
    for(int i=0; i<200; i++){new_rect_idx[i]=0;} //Max allowed clusters=200

    for(int inr=0; inr<nrectangles; inr++){
      int iDivide = 0;
      uint64 clust_id = gCluster.rect_idx[inr];
      iDivide= getSubRectangles(event, module_no,itrn, clust_id, inew_count, ndivision, new_rect_idx);
    }//for(int inr=0; inr<nrectangles; inr++)

    //Copying Local IDX to the Global IDX
    for(int ict=0; ict<inew_count; ict++) {
      gCluster.rect_idx[ict] = new_rect_idx[ict];
    }
    gCluster.nClusters = inew_count;
    if(ndivision==0)break;
  } // end of for(int itrn=0; itrn<3; itrn++)

}

int getSubRectangles(uint event, uint module_no, int itrn, uint64 cluster_id, int &icount, int &ndivision, uint64 *new_rect_idx){
  bool xb[162],    yb[418];   //one element each in array at start and end added
  uint64 rectangle[50][50];
  int  iAction=0;             //Flag if rectangle can be further subdivided
  //int  rect_lidx[200];       //Max allowed clusters=200
  int  lxx[2000],   lyy[2000], lToGlobal[2000],   nlPoints=0, nTPoints=0;
  int  xstart[50],  xend[50],  ystart[50],        yend[50];
  int  nxs=0,       nxe=0,     nys=0,             nye=0;

  nTPoints = gCluster.nPixHits;

  //Intialise lxx, lyy and lToGlobal to unphysical value
  for(int i=0; i<nTPoints+2; i++){lxx[i]=-1; lyy[i]=-1; lToGlobal[i]=-1;}

  //Transform Module cordinates to the respective rectangle corrdinates 
  for(int i=0; i<nTPoints; i++){
    if(gCluster.rect_idxPt[i] != cluster_id) continue;
    lxx[nlPoints]=gCluster.xx[i];   lyy[nlPoints]=gCluster.yy[i];
    lToGlobal[nlPoints]=i;            nlPoints++;
  }

  if(nlPoints==0){iAction=-9; return iAction;}

  //Initialise projected x and y-cordinate bollean array
  for(int i=0; i<418; i++){yb[i]=0;}
  for(int i=0; i<162; i++){xb[i]=0;}

  //Initialise start and end location of each strip on x and y-axis
  for(int i=0; i<50; i++){
    xstart[i]=0; xend[i]=0; ystart[i]=0; yend[i]=0;
  }

  //Take projection of all points on X (and Y) axis
  for(int i=0; i<nlPoints; i++){
    int ix = lxx[i];  int iy = lyy[i];
    xb[ix]=1;         yb[iy]=1;
  }
  xb[0] =0; yb[0] =0;

  //Get the Start and End location of projected x-cordinates of all hits
  for(int i=1; i<162; i++){
    if(!xb[i])continue;
    if(!xb[i-1]){xstart[nxs]=i; nxs++;}
    if(!xb[i+1]){xend[nxe]=i;   nxe++;}
  }

  //Get the Start and End location of projected y-cordinates of all hits
  for(int i=1; i<418; i++){
    if(!yb[i])continue;
    if(!yb[i-1]){ystart[nys]=i; nys++;}
    if(!yb[i+1]){yend[nye]=i;   nye++;}
  }
  
  //Rectangels cannot be further suvdivided
  if(nxs==1 && nys==1){
    iAction=-1;
    
    uint64 new_cluster_id= 0;
    if(itrn==0) {
      if(xstart[0]!=0 && ystart[0]!=0) //cc4
        new_cluster_id = getClusterId(event, module_no,ystart[0], xstart[0]);
    }
    else {
      new_cluster_id = cluster_id;
    } 
    new_rect_idx[icount]=new_cluster_id; 
    for(int iPt=0; iPt<nlPoints; iPt++){
      int globalID = lToGlobal[iPt];
      if(lxx[iPt]==0 || lyy[iPt]==0) continue;
      gCluster.rect_idxPt[globalID] = new_cluster_id; 
    } 
    icount++;
    return iAction;
  }

  iAction=1;
  //Initialise boolean rectangles[][]
  for(int i=0; i<50;   i++){
    for(int j=0; j<50; j++){rectangle[i][j]=0;}
  }
     
  uint64 new_cluster_id = 0;
  for(int iPt=0; iPt<nlPoints; iPt++){
    int globalID = lToGlobal[iPt];
    int ix  = lxx[iPt]; int iy  = lyy[iPt];
    int irx = -1;       int iry = -1;
    if(ix==0 || iy==0) continue; // avoid bad pixel
    for(int i=0; i<nxs; i++){
      if(xstart[i]<= ix && xend[i]>=ix){irx=i;  break;}
    }
    for(int i=0; i<nys; i++){
      if(ystart[i]<= iy && yend[i]>=iy){iry=i; break;}
    }
    
    if(!rectangle[irx][iry]){
      new_rect_idx[icount]=new_cluster_id; 
      if(xstart[irx]!=0 && ystart[iry]!=0) // ignore the bad pixel
        new_cluster_id = getClusterId(event, module_no,ystart[iry], xstart[irx]);

      rectangle[irx][iry]  = new_cluster_id;
      icount++;  ndivision++;  //Number of actual subdivisions
      new_rect_idx[icount-1] = rectangle[irx][iry];
    }
    if(ix!=0 && iy!=0)
      gCluster.rect_idxPt[globalID]   = rectangle[irx][iry];  //Global variable
  }//for(int iPt=0; iPt<nlPoints; iPt++){

  return iAction;
}




