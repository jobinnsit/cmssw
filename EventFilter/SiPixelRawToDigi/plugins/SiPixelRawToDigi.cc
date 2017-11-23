// Skip FED40 pilot-blade
// Include parameter driven interface to SiPixelQuality for study purposes
// exclude ROC(raw) based on bad ROC list in SiPixelQuality
// enabled by: process.siPixelDigis.UseQualityInfo = True (BY DEFAULT NOT USED)
// 20-10-2010 Andrew York (Tennessee)
// Jan 2016 Tamas Almos Vami (Tav) (Wigner RCP) -- Cabling Map label option
// Jul 2017 Viktor Veszpremi -- added PixelFEDChannel

#include <string>
#include <vector>
#include <iterator>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <map>

#include "SiPixelRawToDigi.h"

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

#include "DataFormats/SiPixelDetId/interface/PixelFEDChannel.h"
#include "EventFilter/SiPixelRawToDigi/interface/PixelUnpackingRegions.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "TH1D.h"
#include "TFile.h"


using namespace std;

// -----------------------------------------------------------------------------
SiPixelRawToDigi::SiPixelRawToDigi( const edm::ParameterSet& conf ) 
  : config_(conf), 
    badPixelInfo_(nullptr),
    regions_(nullptr),
    hCPU(nullptr), hDigi(nullptr)
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
    produces<edmNew::DetSetVector<PixelFEDChannel> >();
  }

  // regions
  if (config_.exists("Regions")) {
    if(!config_.getParameter<edm::ParameterSet>("Regions").getParameterNames().empty())
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
    if(usePhase1) edm::LogInfo("SiPixelRawToDigi")  << " Using phase1";
  }
  //CablingMap could have a label //Tav
  cablingMapLabel = config_.getParameter<std::string> ("CablingMapLabel");

  // to generate cabling map required for gpu raw to digi,can be moved in config file
  cablingMapGPU = true;

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

void SiPixelRawToDigi::generateCablingMapGPU(const string file) {
  typedef std::map<SiPixelRawToDigi::Key,DetId> MAP;
  MAP theMap;
  std::map<SiPixelRawToDigi::Key,DetId>::iterator it2;

  const string FED = "FED:";
  const string LNK = "LNK:";
  const string ROCs = "ROCs:";
  const string BARREL = "barrel";
  const string ENDCAP = "endcap";
  const string header = " FedID    LINK  idinLNK   B_F      Rawid     idinDU";
  char  buffer[1024];

  string buf; // Have a buffer string
  string fedstr,    lnkstr;
  string idinLNK_str, BigId_str, idinDU_str;
  int    idinLNK, idinDU;

  int PixRawID[30000];
  for(int i=0; i<30000; i++){PixRawID[i]=0;}

  //Part 1 of the code
  ifstream inMap;  inMap.open("Map_Phase1_d4.txt");
  ofstream outMap; outMap.open("FED_Map_Phase1.dat");
  outMap << header << endl;

  int FedRawIdx=0;
  while(!inMap.eof()) {
    inMap >> buf;
    if(buf==FED) {inMap >> fedstr;}
    if(buf==LNK) {inMap >> lnkstr;}
    if(buf=="========") {
      //======== PixelROC  endcap  (344014096) idInDU: 0 idInLk: 1
      string pixelRoc, type;
      inMap >> pixelRoc >> type >> BigId_str >> idinDU_str >> idinDU >> idinLNK_str >> idinLNK;
      char temp1[12];
      strcpy(temp1, BigId_str.c_str());
      char temp2[12];
      strcpy(temp2,"ABCDEFGHI");
      for(int i=1;i<=9;i++) {temp2[i-1]= temp1[i];}
      string newBigId_str(temp2);
      int B_F; // B_F = 1 if barrel, B_F =2 if endcap
      if(type=="endcap") {B_F = 2;}
      else {B_F = 1;}

      int i1 = boost::lexical_cast<int>(fedstr);
      int i2 = boost::lexical_cast<int>(lnkstr);
      int i3 = idinLNK;    int i4 = B_F;
      int i5=0; sscanf(temp2,"%d", &i5);  //rawid
      int i6 = idinDU;
      sprintf(buffer, "  %04d     %02d     %02d       %1d     %9d     %02d \n", i1, i2, i3, i4, i5, i6);
      outMap << buffer;
      PixRawID[FedRawIdx] = i5; FedRawIdx++;
    }
  }//while(!inMap.eof()) {
  inMap.close();
  outMap << endl; outMap.close();

  //Part 2 of the Code
  sort(PixRawID, PixRawID+FedRawIdx);
  ofstream outMap2;
  outMap2.open("RawID_ModuleID_Phase1.dat");
  const string header2 = "    RawID      ModuleID";
  outMap2 << header2 << endl;
  std::map<int, int> Map;
  std::map<int, int>::iterator it;

  int num_prev=-1; int imod=0;
  for(int i=0; i<30000; i++){
    if(num_prev!=PixRawID[i] && PixRawID[i]>0){
      sprintf(buffer, "  %9d     %04d \n", PixRawID[i], imod);
      Map.insert(std::pair<int, int>(PixRawID[i], imod));
      outMap2 << buffer;
      num_prev=PixRawID[i]; imod++;
    }
  }
  outMap2.close();

  //Part 3 of the Code
  int FedID, Link, idInLnk, B_F, rawId, idInDU, module;
  string str, str1="    ModuleID";

  ifstream mapFile;   mapFile.open("FED_Map_Phase1.dat");
  ofstream mapFile2;  mapFile2.open("FED_Map_With_ModuleID_Phase1.dat");

  getline(mapFile, str);
  mapFile2 << str+str1 <<endl;

  while(!mapFile.eof()) {
    mapFile>> FedID >> Link >> idInLnk >> B_F >> rawId >> idInDU;
    it = Map.find(rawId);
    if(it !=Map.end()) {module = it->second;}
    else {module = -999;}
    int i1=FedID; int i2=Link; int i3=idInLnk;  int i4=B_F; int i5=rawId; int i6=idInDU; int i7=module;
    sprintf(buffer, "  %04d     %02d     %02d       %1d     %9d     %02d        %04d\n", i1, i2, i3, i4, i5, i6, i7);
    mapFile2 << buffer;
  }
  mapFile.close();  mapFile2.close();


  //Part 4 (Final Stage) of the Code
  ifstream fedFile;
  fedFile.open("FED_Map_With_ModuleID_Phase1.dat");
  ofstream  arrayFile;
  arrayFile.open("Pixel_Phase1_Raw2Digi_GPU_Cabling_Map.dat");

  getline(fedFile, str);
  unsigned int uFedID,  uLink, uidinLnk,  uRawId, uidinDU, uModuleId;
  while(!fedFile.eof()){
    fedFile >> uFedID >> uLink >> uidinLnk >> B_F >> uRawId >> uidinDU >> uModuleId;
    Key key(uFedID, uLink, uidinLnk);
    DetId detId(uRawId, uidinDU, uModuleId);
    theMap.insert(std::pair<Key, DetId>(key,detId));
  }
  fedFile.close();

  const string header3 = "  Index      FedID  Link  idinLNK   B_F      RawID     idinDU    ModuleID";
  arrayFile<< header3 <<endl;
  int idxF=0;
  for(unsigned int fed=1200;fed<=1338;fed++) {
    for(unsigned int link=1;link<=48;link++) {
      for(unsigned int roc=1;roc<=8;roc++) {
        idxF++;
        Key key(fed, link, roc);
        it2 = theMap.find(key);
        if(it2 !=theMap.end()) {
           uRawId = it2->second.RawId;
           uModuleId = it2->second.ModuleId;
           uidinDU   = it2->second.idinDU;
           int b_f=2;                //Forward
           if(uModuleId<1184)b_f=1;  //Barrel
           sprintf(buffer, "  %06d     %04d    %02d     %02d       %1d     %9d     %02d        %04d\n", 
           idxF, fed, link, roc, b_f, uRawId, uidinDU, uModuleId);
           arrayFile << buffer;
        }
        else {
          //This combination of FED, Link and ROC is not present 
          //in the cabling map. Therefore, fill with dummy 9999
          int i1=9999, i2=99, i3=0;
          sprintf(buffer, "  %06d     %04d    %02d     %02d       %1d     %9d     %02d        %04d\n", 
          idxF, i1, i2, i2, i3, i1, i2, i1);
          arrayFile << buffer;
        }
      }  //for(unsigned int roc=1;roc<=8;roc++) {
    }   //for(unsigned int link=1;link<=48;link++) {
  }    //for(unsigned int fed=1200;fed<=1303;fed++) {

  //Given FedId, Link and idinLnk; use the following formula
  //to get the RawId and idinDU 
  //index = (FedID-1200) * MAX_LINK* MAX_ROC + (Link-1)* MAX_ROC + idinLnk;  
  //where, MAX_LINK = 48, MAX_ROC = 8 for Phase1 as mentioned Danek's email
  //FedID varies between 1200 to 1338 (In total 108 FED's)
  //Link varies between 1 to 48
  //idinLnk varies between 1 to 8

  cout << "Cabling map for GPU generated Succesfully!" << endl;
  arrayFile.close();
}


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
    if(cablingMapGPU) {
      const string file = "Map_Phase1_d4.txt";
      ofstream ofile(file);
      ofile << cabling_->print(4); // print(depth), depth = 4 prints feds, links, roc
      ofile.close();
      generateCablingMapGPU(file);
    }
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
  auto disabled_channelcollection = std::make_unique<edmNew::DetSetVector<PixelFEDChannel> >();

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

  for (auto aFed = fedIds.begin(); aFed != fedIds.end(); ++aFed) {
    int fedId = *aFed;

    if(!usePilotBlade && (fedId==40) ) continue; // skip pilot blade data

    if (regions_ && !regions_->mayUnpackFED(fedId)) continue;

    if(debug) LogDebug("SiPixelRawToDigi")<< " PRODUCE DIGI FOR FED: " <<  fedId << endl;

    PixelDataFormatter::Errors errors;

    //get event data for this fed
    const FEDRawData& fedRawData = buffers->FEDData( fedId );

    //convert data to digi and strip off errors
    formatter.interpretRawData( errorsInEvent, fedId, fedRawData, *collection, errors);

    //pack errors into collection
    if(includeErrors) {
      typedef PixelDataFormatter::Errors::iterator IE;
      for (IE is = errors.begin(); is != errors.end(); is++) {
	uint32_t errordetid = is->first;
	if (errordetid==dummydetid) {           // errors given dummy detId must be sorted by Fed
	  nodeterrors.insert( nodeterrors.end(), errors[errordetid].begin(), errors[errordetid].end() );
	} else {
	  edm::DetSet<SiPixelRawDataError>& errorDetSet = errorcollection->find_or_insert(errordetid);
	  errorDetSet.data.insert(errorDetSet.data.end(), is->second.begin(), is->second.end());
	  // Fill detid of the detectors where there is error AND the error number is listed
	  // in the configurable error list in the job option cfi.
	  // Code needs to be here, because there can be a set of errors for each 
	  // entry in the for loop over PixelDataFormatter::Errors

	  std::vector<PixelFEDChannel> disabledChannelsDetSet;

	  for (auto const& aPixelError : errorDetSet) {
	    // For the time being, we extend the error handling functionality with ErrorType 25
	    // In the future, we should sort out how the usage of tkerrorlist can be generalized
	    if (aPixelError.getType()==25) {
	      assert(aPixelError.getFedId()==fedId);
	      const sipixelobjects::PixelFEDCabling* fed = cabling_->fed(fedId);
	      if (fed) {
		cms_uint32_t linkId = formatter.linkId(aPixelError.getWord32());
		const sipixelobjects::PixelFEDLink* link = fed->link(linkId);
		if (link) {
		  // The "offline" 0..15 numbering is fixed by definition, also, the FrameConversion depends on it
		  // in contrast, the ROC-in-channel numbering is determined by hardware --> better to use the "offline" scheme
		  PixelFEDChannel ch = {fed->id(), linkId, 25, 0};
		  for (unsigned int iRoc=1; iRoc<=link->numberOfROCs(); iRoc++) {
		    const sipixelobjects::PixelROC * roc = link->roc(iRoc);
		    if (roc->idInDetUnit()<ch.roc_first) ch.roc_first=roc->idInDetUnit();
		    if (roc->idInDetUnit()>ch.roc_last) ch.roc_last=roc->idInDetUnit();
		  }
		  disabledChannelsDetSet.push_back(ch);
		}
	      }
	    } else {
	      // fill list of detIds to be turned off by tracking
	      if(!tkerrorlist.empty()) {
		std::vector<int>::iterator it_find = find(tkerrorlist.begin(), tkerrorlist.end(), aPixelError.getType());
		if(it_find != tkerrorlist.end()){
		  tkerror_detidcollection->push_back(errordetid);
		}
	      }
	    }

	    // fill list of detIds with errors to be studied
	    if(!usererrorlist.empty()) {
	      std::vector<int>::iterator it_find = find(usererrorlist.begin(), usererrorlist.end(), aPixelError.getType());
	      if(it_find != usererrorlist.end()){
		usererror_detidcollection->push_back(errordetid);
	      }
	    }

	  } // loop on DetSet of errors

	  if (!disabledChannelsDetSet.empty()) {
	    disabled_channelcollection->insert(errordetid, disabledChannelsDetSet.data(), disabledChannelsDetSet.size());
	  }

	} // if error assigned to a real DetId
      } // loop on errors in event for this FED
    } // if errors to be included in the event
  } // loop on FED data to be unpacked

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
    ev.put(std::move(disabled_channelcollection));
  }
}
