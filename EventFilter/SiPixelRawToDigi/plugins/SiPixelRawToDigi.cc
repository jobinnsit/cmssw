// Skip FED40 pilot-blade
// Include parameter driven interface to SiPixelQuality for study purposes
// exclude ROC(raw) based on bad ROC list in SiPixelQuality
// enabled by: process.siPixelDigis.UseQualityInfo = True (BY DEFAULT NOT USED)
// 20-10-2010 Andrew York (Tennessee)
// Jan 2016 Tamas Almos Vami (Tav) (Wigner RCP) -- Cabling Map label option

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"

#include "DataFormats/SiPixelRawData/interface/SiPixelRawDataError.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDRawData.h"

#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/DetId/interface/DetIdCollection.h"
#include "CondFormats/SiPixelObjects/interface/SiPixelFedCablingMap.h"
#include "CondFormats/SiPixelObjects/interface/SiPixelFedCablingTree.h"
#include "EventFilter/SiPixelRawToDigi/interface/PixelDataFormatter.h"

#include "CondFormats/SiPixelObjects/interface/SiPixelQuality.h"

#include "EventFilter/SiPixelRawToDigi/interface/PixelUnpackingRegions.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "TH1D.h"
#include "TFile.h"
#include "SiPixelRawToDigi.h"
#include <string>
#include <chrono>
#include <iostream>
#include <fstream>
// for GPU
// device memory intialization for RawTodigi
#include "RawToDigiMem.h"
// device memory initialization for CPE
#include "CPEGPUMem.h"
//device memory initialization for clustering
#include "PixelClusterMem.h"

using namespace std;

// -----------------------------------------------------------------------------
SiPixelRawToDigi::SiPixelRawToDigi( const edm::ParameterSet& conf ) 
  : config_(conf), 
    badPixelInfo_(0),
    regions_(0),
    hCPU(0), hDigi(0)
{

  includeErrors = config_.getParameter<bool>("IncludeErrors");
  useQuality = config_.getParameter<bool>("UseQualityInfo");
  if (config_.exists("ErrorList")) {
    tkerrorlist = config_.getParameter<std::vector<int> > ("ErrorList");
  }
  if (config_.exists("UserErrorList")) {
    usererrorlist = config_.getParameter<std::vector<int> > ("UserErrorList");
  }
  tFEDRawDataCollection = consumes <FEDRawDataCollection> (config_.getParameter<edm::InputTag>("InputLabel"));

  //start counters
  ndigis = 0;
  nwords = 0;

  // Products
  produces< edm::DetSetVector<PixelDigi> >();
  if(includeErrors){
    produces< edm::DetSetVector<SiPixelRawDataError> >();
    produces<DetIdCollection>();
    produces<DetIdCollection>("UserErrorModules");
  }

  // regions
  if (config_.exists("Regions")) {
    if(config_.getParameter<edm::ParameterSet>("Regions").getParameterNames().size() > 0)
    {
      regions_ = new PixelUnpackingRegions(config_, consumesCollector());
    }
  }

  // Timing
  bool timing = config_.getUntrackedParameter<bool>("Timing",false);
  if (timing) {
    theTimer.reset( new edm::CPUTimer );
    hCPU = new TH1D ("hCPU","hCPU",100,0.,0.050);
    hDigi = new TH1D("hDigi","hDigi",50,0.,15000.);
  }

  // Control the usage of pilot-blade data, FED=40
  usePilotBlade = false; 
  if (config_.exists("UsePilotBlade")) {
    usePilotBlade = config_.getParameter<bool> ("UsePilotBlade");
    if(usePilotBlade) edm::LogInfo("SiPixelRawToDigi")  << " Use pilot blade data (FED 40)";
  }

  // Control the usage of phase1
  usePhase1 = false;
  if (config_.exists("UsePhase1")) {
    usePhase1 = config_.getParameter<bool> ("UsePhase1");
    if(usePhase1) edm::LogInfo("SiPixelRawToDigi")  << " Use pilot blade data (FED 40)";
  }
  //CablingMap could have a label //Tav
  cablingMapLabel = config_.getParameter<std::string> ("CablingMapLabel");
  
  //GPU specific
  const int MAX_FED  = 150;
  const int MAX_WORD = 2000;
  word = (unsigned int*)malloc(MAX_FED*MAX_WORD*sizeof(unsigned int));
  fedIndex =(unsigned int*)malloc(2*(MAX_FED+1)*sizeof(unsigned int));
  // allocate memory for RawToDigi on GPU
  initDeviceMemory();
  // allocate memory for CPE on GPU
  initDeviceMemCPE();
  // allocate auxilary memory for clustering
  initDeviceMemCluster();
}


// -----------------------------------------------------------------------------
SiPixelRawToDigi::~SiPixelRawToDigi() {
  edm::LogInfo("SiPixelRawToDigi")  << " HERE ** SiPixelRawToDigi destructor!";

  if (regions_) delete regions_;

  if (theTimer) {
    TFile rootFile("analysis.root", "RECREATE", "my histograms");
    hCPU->Write();
    hDigi->Write();
  }
  free(word);
  free(fedIndex);
  // free device memory used for RawToDigi on GPU
  freeMemory(); 
  // free device memory used for CPE on GPU
  freeDeviceMemCPE();
  // free auxilary memory used for clustering
  freeDeviceMemCluster();
}

void
SiPixelRawToDigi::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<bool>("IncludeErrors",true);
  desc.add<bool>("UseQualityInfo",false);
  {
    std::vector<int> temp1;
    temp1.reserve(1);
    temp1.push_back(29);
    desc.add<std::vector<int> >("ErrorList",temp1)->setComment("## ErrorList: list of error codes used by tracking to invalidate modules");
  }
  {
    std::vector<int> temp1;
    temp1.reserve(1);
    temp1.push_back(40);
    desc.add<std::vector<int> >("UserErrorList",temp1)->setComment("## UserErrorList: list of error codes used by Pixel experts for investigation");
  }
  desc.add<edm::InputTag>("InputLabel",edm::InputTag("siPixelRawData"));
  {
    edm::ParameterSetDescription psd0;
    psd0.addOptional<std::vector<edm::InputTag>>("inputs");
    psd0.addOptional<std::vector<double>>("deltaPhi");
    psd0.addOptional<std::vector<double>>("maxZ");
    psd0.addOptional<edm::InputTag>("beamSpot");
    desc.add<edm::ParameterSetDescription>("Regions",psd0)->setComment("## Empty Regions PSet means complete unpacking");
  }
  desc.addUntracked<bool>("Timing",false);
  desc.add<bool>("UsePilotBlade",false)->setComment("##  Use pilot blades");
  desc.add<bool>("UsePhase1",false)->setComment("##  Use phase1");
  desc.add<std::string>("CablingMapLabel","")->setComment("CablingMap label"); //Tav
  desc.addOptional<bool>("CheckPixelOrder");  // never used, kept for back-compatibility
  descriptions.add("siPixelRawToDigi",desc);
}

// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
void SiPixelRawToDigi::produce( edm::Event& ev,
                              const edm::EventSetup& es) 
{
  const uint32_t dummydetid = 0xffffffff;
  debug = edm::MessageDrop::instance()->debugEnabled;

// initialize cabling map or update if necessary
  if (recordWatcher.check( es )) {
    // cabling map, which maps online address (fed->link->ROC->local pixel) to offline (DetId->global pixel)
    edm::ESTransientHandle<SiPixelFedCablingMap> cablingMap;
    es.get<SiPixelFedCablingMapRcd>().get( cablingMapLabel, cablingMap ); //Tav
    fedIds   = cablingMap->fedIds();
    cabling_ = cablingMap->cablingTree();
    LogDebug("map version:")<< cabling_->version();
  }
// initialize quality record or update if necessary
  if (qualityWatcher.check( es )&&useQuality) {
    // quality info for dead pixel modules or ROCs
    edm::ESHandle<SiPixelQuality> qualityInfo;
    es.get<SiPixelQualityRcd>().get( qualityInfo );
    badPixelInfo_ = qualityInfo.product();
    if (!badPixelInfo_) {
      edm::LogError("SiPixelQualityNotPresent")<<" Configured to use SiPixelQuality, but SiPixelQuality not present"<<endl;
    }
  }

  edm::Handle<FEDRawDataCollection> buffers;
  ev.getByToken(tFEDRawDataCollection, buffers);

// create product (digis & errors)
  auto collection = std::make_unique<edm::DetSetVector<PixelDigi>>();
  // collection->reserve(8*1024);
  auto errorcollection = std::make_unique<edm::DetSetVector<SiPixelRawDataError>>();
  auto tkerror_detidcollection = std::make_unique<DetIdCollection>();
  auto usererror_detidcollection = std::make_unique<DetIdCollection>();

  //PixelDataFormatter formatter(cabling_.get()); // phase 0 only
  PixelDataFormatter formatter(cabling_.get(), usePhase1); // for phase 1 & 0

  formatter.setErrorStatus(includeErrors);

  if (useQuality) formatter.setQualityStatus(useQuality, badPixelInfo_);

  if (theTimer) theTimer->start();
  bool errorsInEvent = false;
  PixelDataFormatter::DetErrors nodeterrors;

  if (regions_) {
    regions_->run(ev, es);
    formatter.setModulesToUnpack(regions_->modulesToUnpack());
    LogDebug("SiPixelRawToDigi") << "region2unpack #feds: "<<regions_->nFEDs();
    LogDebug("SiPixelRawToDigi") << "region2unpack #modules (BPIX,EPIX,total): "<<regions_->nBarrelModules()<<" "<<regions_->nForwardModules()<<" "<<regions_->nModules();
  }
  // GPU specific 
  // data extraction for RawToDigi GPU
  unsigned int wordCounterGPU =0;
  unsigned int fedCounter =0;
  const unsigned int MAX_FED = 150;
  ErrorChecker errorcheck;
  for (auto aFed = fedIds.begin(); aFed != fedIds.end(); ++aFed) {
    int fedId = *aFed;
    // cout<<"FedId: "<<fedId<<endl;
    // for GPU
    // first 150 index stores the fedId and next 150 will store the
    // start index of word in that fed
    fedIndex[fedCounter] = fedId-1200;
    fedIndex[MAX_FED + fedCounter] = wordCounterGPU; // MAX_FED = 150
    fedCounter++;

    //get event data for this fed
    const FEDRawData& rawData = buffers->FEDData( fedId );
    //GPU specific 
    PixelDataFormatter::Errors errors;
    int nWords = rawData.size()/sizeof(Word64);
    if(nWords==0) {
      word[wordCounterGPU++] =0;
      continue;
    }  

    // check CRC bit
    const Word64* trailer = reinterpret_cast<const Word64* >(rawData.data())+(nWords-1);  
    if(!errorcheck.checkCRC(errorsInEvent, fedId, trailer, errors)) {
      word[wordCounterGPU++] =0;
      continue;
    }

    // check headers
    const Word64* header = reinterpret_cast<const Word64* >(rawData.data()); header--;
    bool moreHeaders = true;
    while (moreHeaders) {
      header++;
      //LogTrace("")<<"HEADER:  " <<  print(*header);
      bool headerStatus = errorcheck.checkHeader(errorsInEvent, fedId, header, errors);
      moreHeaders = headerStatus;
    }

    // check trailers
    bool moreTrailers = true;
    trailer++;
    while (moreTrailers) {
      trailer--;
      //LogTrace("")<<"TRAILER: " <<  print(*trailer);
      bool trailerStatus = errorcheck.checkTrailer(errorsInEvent, fedId, nWords, trailer, errors);
      moreTrailers = trailerStatus;
    }

    // data words
    //LogTrace("")<<"data words: "<< (trailer-header-1);

    const  Word32 * bw =(const  Word32 *)(header+1);
    const  Word32 * ew =(const  Word32 *)(trailer);
    if ( *(ew-1) == 0 ) { ew--; }
    for (auto ww = bw; ww < ew; ++ww) {
      //LogTrace("")<<"DATA: " <<  print(*word);
      word[wordCounterGPU++] = *ww;
    }
  }  // end of for loop
  

  // original for loop
 /*
  for (auto aFed = fedIds.begin(); aFed != fedIds.end(); ++aFed) {
    int fedId = *aFed;
    //cout<<"FedId: "<<fedId<<endl;
    // for GPU
    // first 150 index stores the fedId and next 150 will store the
    // start index of word in that fed
    fedIndex[fedCounter] = fedId-1200;
    fedIndex[MAX_FED + fedCounter] = wordCounterGPU; // MAX_FED = 150
    fedCounter++;
    if(!usePilotBlade && (fedId==40) ) continue; // skip pilot blade data

    if (regions_ && !regions_->mayUnpackFED(fedId)) continue;

    if(debug) LogDebug("SiPixelRawToDigi")<< " PRODUCE DIGI FOR FED: " <<  fedId << endl;

    PixelDataFormatter::Errors errors;

    //get event data for this fed
    const FEDRawData& fedRawData = buffers->FEDData( fedId );
    //GPU specific
    
    //convert data to digi and strip off errors
    formatter.interpretRawData( errorsInEvent, fedId, fedRawData, *collection, errors,word, wordCounterGPU);
   
    //pack errors into collection
    if(includeErrors) {
      typedef PixelDataFormatter::Errors::iterator IE;
      for (IE is = errors.begin(); is != errors.end(); is++) {
        uint32_t errordetid = is->first;
        if (errordetid==dummydetid) {           // errors given dummy detId must be sorted by Fed
          nodeterrors.insert( nodeterrors.end(), errors[errordetid].begin(), errors[errordetid].end() );
        } 
        else {
          edm::DetSet<SiPixelRawDataError>& errorDetSet = errorcollection->find_or_insert(errordetid);
          errorDetSet.data.insert(errorDetSet.data.end(), is->second.begin(), is->second.end());
          // Fill detid of the detectors where there is error AND the error number is listed
          // in the configurable error list in the job option cfi.
          // Code needs to be here, because there can be a set of errors for each 
          // entry in the for loop over PixelDataFormatter::Errors
          if(!tkerrorlist.empty() || !usererrorlist.empty()){
            DetId errorDetId(errordetid);
            edm::DetSet<SiPixelRawDataError>::const_iterator itPixelError=errorDetSet.begin();
            for(; itPixelError!=errorDetSet.end(); ++itPixelError){
              // fill list of detIds to be turned off by tracking
              if(!tkerrorlist.empty()) {
                std::vector<int>::iterator it_find = find(tkerrorlist.begin(), tkerrorlist.end(), itPixelError->getType());
               if(it_find != tkerrorlist.end()){
              tkerror_detidcollection->push_back(errordetid);
               }
            }
              // fill list of detIds with errors to be studied
          if(!usererrorlist.empty()) {
            std::vector<int>::iterator it_find = find(usererrorlist.begin(), usererrorlist.end(), itPixelError->getType());
            if(it_find != usererrorlist.end()){
            usererror_detidcollection->push_back(errordetid);
            }
          }
         }
       }
      }
     }
    } // end of includeErrors
     
  }  // end of for loop
 */ 
  /*
  if(includeErrors) {
    edm::DetSet<SiPixelRawDataError>& errorDetSet = errorcollection->find_or_insert(dummydetid);
    errorDetSet.data = nodeterrors;
  }
  if (errorsInEvent) LogDebug("SiPixelRawToDigi") << "Error words were stored in this event";

  if (theTimer) {
    theTimer->stop();
    LogDebug("SiPixelRawToDigi") << "TIMING IS: (real)" << theTimer->realTime() ;
    ndigis += formatter.nDigis();
    nwords += formatter.nWords();
    LogDebug("SiPixelRawToDigi") << " (Words/Digis) this ev: "
         <<formatter.nWords()<<"/"<<formatter.nDigis() << "--- all :"<<nwords<<"/"<<ndigis;
    hCPU->Fill( theTimer->realTime() ); 
    hDigi->Fill(formatter.nDigis());
  }

  //send digis and errors back to framework 
  ev.put(std::move(collection));
  if(includeErrors){
    ev.put(std::move(errorcollection));
    ev.put(std::move(tkerror_detidcollection));
    ev.put(std::move(usererror_detidcollection), "UserErrorModules");
  }
  */
  //GPU specific
 
  static int eventCount=0;
  eventCount++;
  //cout<<"Event: "<<setw(4)<<eventCount<<"  Total Hits: "<<setw(8)<<wordCounterGPU<<endl;
  // following code is to extract the Raw data put it to ascii file and pass it to
  // standalone RawToDigi
  /*
  ofstream fedEventFile("fedCount_EventFile.dat", ios::out | ios::app);
  fedEventFile<<setw(6)<<fedCounter<<setw(8)<<wordCounterGPU<<endl;
  fedEventFile.close();
  
  // to store the fedId and there index
  ofstream fedIndexFile("fedIndexFile_HT_after.dat", ios::out | ios::app);
  for(unsigned int i=0;i<150*2;i++) {
   // cout<<"fedId: "<<i<<"   Index: "<<fedIndex[150+i]<<endl;
   fedIndexFile<<fedIndex[i]<<endl;
   //fedIndex[i] =0;
  }
  ofstream dataFile("wordDataFile_HT_debug_after.dat", ios::out | ios::app);
  for(unsigned int i=0;i<wordCounterGPU;i++) {
    dataFile<<word[i]<<endl;
  }
  fedIndexFile.close();
  dataFile.close();
  */
  // RawToDigi -> clustering -> CPE
  RawToDigi_kernel_wrapper (wordCounterGPU, word, fedCounter,fedIndex);

	wordCounterGPU = 0;
  fedCounter =0;
} // end of produce function

//define as runnable module
DEFINE_FWK_MODULE(SiPixelRawToDigi);

