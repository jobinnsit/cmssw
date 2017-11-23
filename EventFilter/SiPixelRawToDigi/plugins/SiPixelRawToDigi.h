#ifndef SiPixelRawToDigi_H
#define SiPixelRawToDigi_H

/** \class SiPixelRawToDigi_H
 *  Plug-in module that performs Raw data to digi conversion 
 *  for pixel subdetector
 */

#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CondFormats/DataRecord/interface/SiPixelFedCablingMapRcd.h"
#include "CondFormats/DataRecord/interface/SiPixelQualityRcd.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/CPUTimer.h"

class SiPixelFedCablingTree;
class SiPixelFedCabling;
class SiPixelQuality;
class TH1D;
class PixelUnpackingRegions;

class SiPixelRawToDigi : public edm::stream::EDProducer<> {
public:

  /// ctor
  explicit SiPixelRawToDigi( const edm::ParameterSet& );

  /// dtor
  ~SiPixelRawToDigi() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  // genrate the cabling map for GPU RawToDigi
  void generateCablingMapGPU(const std::string file);

  /// get data, convert to digis attach againe to Event
  void produce( edm::Event&, const edm::EventSetup& ) override;


private:

  edm::ParameterSet config_;
  std::unique_ptr<SiPixelFedCablingTree> cabling_;
  const SiPixelQuality* badPixelInfo_;
  PixelUnpackingRegions* regions_;
  edm::EDGetTokenT<FEDRawDataCollection> tFEDRawDataCollection; 

  TH1D *hCPU, *hDigi;
  std::unique_ptr<edm::CPUTimer> theTimer;
  bool includeErrors;
  bool useQuality;
  bool debug;
  std::vector<int> tkerrorlist;
  std::vector<int> usererrorlist;
  std::vector<unsigned int> fedIds;
  edm::ESWatcher<SiPixelFedCablingMapRcd> recordWatcher;
  edm::ESWatcher<SiPixelQualityRcd> qualityWatcher;
  edm::InputTag label;
  int ndigis;
  int nwords;
  bool usePilotBlade;
  bool usePhase1;
  std::string cablingMapLabel;
  bool cablingMapGPU;

  class Key {
    public:
      Key(unsigned int _fed, unsigned int _link, unsigned int _roc): 
      fed(_fed), link(_link), roc(_roc) {}
      //bool operator < (const Key & other) const;
      bool operator < (const Key & other) const {
        if (this->fed < other.fed) return true;
        if (this->fed > other.fed) return false;
        if (this->link < other.link) return true;
        if (this->link > other.link) return false;
        if (this->roc < other.roc) return true;
        if (this->roc > other.roc) return false; 
        return false;
      }    
    private:
      unsigned int fed;
      unsigned int link;
      unsigned int roc;
  };
  struct DetId {
    unsigned int RawId;
    unsigned int idinDU;
    unsigned int ModuleId;
    DetId(unsigned int _RawId,unsigned int _idinDU, unsigned int _ModuleId): 
    RawId(_RawId), idinDU(_idinDU), ModuleId(_ModuleId) {}
  };
};
#endif
