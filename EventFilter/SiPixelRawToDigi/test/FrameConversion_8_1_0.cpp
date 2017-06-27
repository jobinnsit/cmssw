FrameConversion::FrameConversion(bool bpix, int side, int rocIdInDetUnit) {
   int slopeRow =0;
   int slopeCol = 0;
   int  rowOffset = 0;
   int  colOffset = 0; 
  
   if (bpix ) { // bpix 
     
     if (side==-1) {  // -Z side
 
       if (rocIdInDetUnit <8) {
     slopeRow = 1;
     slopeCol = -1;
     rowOffset = 0;
     colOffset = (8-rocIdInDetUnit)*LocalPixel::numColsInRoc-1;  
       } else {
     slopeRow = -1;
     slopeCol = 1;      
     rowOffset = 2*LocalPixel::numRowsInRoc-1;
     colOffset = (rocIdInDetUnit-8)*LocalPixel::numColsInRoc;
       } // if roc
       
     } else {  // +Z side
       
       if (rocIdInDetUnit <8) {
     slopeRow = -1;
     slopeCol = 1;
     rowOffset = 2*LocalPixel::numRowsInRoc-1;
     colOffset = rocIdInDetUnit * LocalPixel::numColsInRoc; 
       } else {
     slopeRow = 1;
     slopeCol = -1;
     rowOffset = 0;
     colOffset = (16-rocIdInDetUnit)*LocalPixel::numColsInRoc-1; 
       }
       
     } // end if +-Z
 
 
 
   } else { // fpix 
 
     // for fpix follow Urs's code for pilot blade
     // no difference between panels
     if(side==-1) { // pannel 1
       if (rocIdInDetUnit < 8) {
     slopeRow = 1;
     slopeCol = -1;
     rowOffset = 0;
     colOffset = (8-rocIdInDetUnit)*LocalPixel::numColsInRoc-1;
       } else {
     slopeRow = -1;
     slopeCol = 1;
     rowOffset = 2*LocalPixel::numRowsInRoc-1;
     colOffset = (rocIdInDetUnit-8)*LocalPixel::numColsInRoc;
       }
     } else { // pannel 2 
       if (rocIdInDetUnit < 8) {
     slopeRow = 1;
     slopeCol = -1;
     rowOffset = 0;
     colOffset = (8-rocIdInDetUnit)*LocalPixel::numColsInRoc-1;
       } else {
     slopeRow = -1;
     slopeCol = 1;
     rowOffset = 2*LocalPixel::numRowsInRoc-1;
     colOffset = (rocIdInDetUnit-8)*LocalPixel::numColsInRoc;
       }
       
     } // side 
 
   } // bpix/fpix
 
   theRowConversion      = LinearConversion(rowOffset,slopeRow);
   theCollumnConversion =  LinearConversion(colOffset, slopeCol);
   
 }