// -*- C++ -*-
//
// Package:    HSCPAnalysis/muSimHits
// Class:      muSimHits
// 
/**\class muSimHits muSimHits.cc HSCPAnalysis/muSimHits/plugins/muSimHits.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Gabriel Ramirez Sanchez
//         Created:  Thu, 22 Jun 2017 08:45:07 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include <DataFormats/MuonDetId/interface/GEMDetId.h>
#include <DataFormats/MuonDetId/interface/DTLayerId.h>
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/DetId/interface/DetId.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class muSimHits : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit muSimHits(const edm::ParameterSet&);
      ~muSimHits();

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      
      edm::EDGetTokenT<edm::PSimHitContainer> simHitRPCtoken_;
      edm::EDGetTokenT<edm::PSimHitContainer> simHitDTtoken_;
      edm::EDGetTokenT<edm::PSimHitContainer> simHitCSCtoken_;
      edm::EDGetTokenT<edm::PSimHitContainer> simHitGEMtoken_;
      edm::EDGetTokenT<std::vector<reco::GenParticle>> genParToken_;
      int pId_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
muSimHits::muSimHits(const edm::ParameterSet& iConfig)
: simHitRPCtoken_(consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("simHitRPCLabel"))),
  simHitDTtoken_(consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("simHitDTLabel"))),
  simHitCSCtoken_(consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("simHitCSCLabel"))),
  simHitGEMtoken_(consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("simHitGEMLabel"))),
  genParToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticleLabel"))),
  pId_(iConfig.getParameter<int>("particleId"))
{
   //now do what ever initialization is needed
   usesResource("TFileService");

}


muSimHits::~muSimHits()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
muSimHits::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;

  Handle<std::vector<reco::GenParticle>> genParticlesHandle;
  iEvent.getByToken(genParToken_, genParticlesHandle);

  Handle<PSimHitContainer> simHitsRPCHandle;
  iEvent.getByToken(simHitRPCtoken_, simHitsRPCHandle);

  Handle<PSimHitContainer> simHitsDTHandle;
  iEvent.getByToken(simHitDTtoken_, simHitsDTHandle);

  Handle<PSimHitContainer> simHitsCSCHandle;
  iEvent.getByToken(simHitCSCtoken_, simHitsCSCHandle);

  Handle<PSimHitContainer> simHitsGEMHandle;
  iEvent.getByToken(simHitGEMtoken_, simHitsGEMHandle);

  if(genParticlesHandle.isValid())
  {
    for(auto& genParticle : *genParticlesHandle)
    {
      cout << genParticle.p4().Eta() << endl; 
    }
  }
  
  
  if(simHitsRPCHandle.isValid())
  {
    for(auto& simHitRPC : *simHitsRPCHandle)
    {
      
      int particleId = simHitRPC.particleType();
      if( particleId != pId_) continue;
      
      DetId detId = DetId(simHitRPC.detUnitId());
      if(detId.subdetId() == MuonSubdetId::RPC)
      {
        RPCDetId rollId(detId);
        cout << simHitRPC.particleType() << endl;
        cout << rollId.station() << endl;
        cout << "is iRPC: " << roll->isIRPC();
      }
    }
  }

  if(simHitsDTHandle.isValid())
  {
    for(auto& simHitDT : *simHitsDTHandle)
    {
      int particleId = simHitDT.particleType();
      if(particleId != pId_) continue;

      DetId detId= DetId(simHitDT.detUnitId());
      if(detId.subdetId() == MuonSubdetId::DT)
      {
        DTLayerId layerId(detId);
        cout << "DT: " << dtLayerId.station() <<endl;
      }
    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
muSimHits::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
muSimHits::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(muSimHits);
