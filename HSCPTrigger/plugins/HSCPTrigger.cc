// -*- C++ -*-
//
// Package:    test/HSCPTrigger
// Class:      HSCPTrigger
// 
/**\class HSCPTrigger HSCPTrigger.cc test/HSCPTrigger/plugins/HSCPTrigger.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Gabriel Ramirez Sanchez
//         Created:  Tue, 04 Apr 2017 12:00:17 GMT
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

#include "FWCore/Framework/interface/ESHandle.h"

#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "SimDataFormats/RPCDigiSimLink/interface/RPCDigiSimLink.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/DetSetVector.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include <Geometry/RPCGeometry/interface/RPCRoll.h>
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/CommonDetUnit/interface/GeomDet.h>
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/RPCDigi/interface/RPCDigi.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "TMath.h"
#include "TH1.h"
#include "TTree.h"
#include <vector>
#include "TVector3.h"
#include "TLorentzVector.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class HSCPTrigger : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit HSCPTrigger(const edm::ParameterSet&);
      ~HSCPTrigger();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
   	UInt_t getTriggerBits(const edm::Event &, std::vector<std::string>);
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      Float_t betaGen;
		Float_t etaGen;
		
		TH1D* fHbeta_pas;
		TH1D* fHbeta_tot;
		
		std::vector<std::string> triggNames_;
      edm::EDGetTokenT<std::vector<reco::GenParticle>> genParToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerToken_;

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
HSCPTrigger::HSCPTrigger(const edm::ParameterSet& iConfig)
: 	triggNames_(iConfig.getParameter<std::vector<std::string>>("triggNames")),
	genParToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticlesLabel"))),
	triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerLabel")))
{
   //now do what ever initialization is needed
   usesResource("TFileService");

}


HSCPTrigger::~HSCPTrigger()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
UInt_t HSCPTrigger::getTriggerBits(const edm::Event& iEvent, std::vector<std::string> TestFilterNames_)
{
   UInt_t trigger=0;
   edm::Handle<edm::TriggerResults> hltTriggerResultHandle;
   iEvent.getByToken(triggerToken_,hltTriggerResultHandle);
   
   if(hltTriggerResultHandle.isValid())
   {
      const edm::TriggerNames & triggerNames = iEvent.triggerNames(*hltTriggerResultHandle);
      for(unsigned i=0 ; i< TestFilterNames_.size();i++)
      {
	      unsigned int bit = triggerNames.triggerIndex(edm::InputTag(TestFilterNames_[i].c_str()).label().c_str());
      	if(bit < hltTriggerResultHandle->size() && hltTriggerResultHandle->accept(bit) && !hltTriggerResultHandle->error(bit))
         {
         	//trigger += 1 << bit;
            trigger += 1 << i;
            break;
         }
      }
   } else std::cout << "HSCPTrigger::getTriggerBits: **** No triggers found ****" << std::endl;
   return trigger;
}


// ------------ method called for each event  ------------
void
HSCPTrigger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   //edm::Handle<edm::TriggerResults> trigResults; //our trigger result object
	
	edm::Handle<vector<reco::GenParticle>> genParHandle;
	iEvent.getByToken(genParToken_,genParHandle);
   // Access Trigger Results
   UInt_t trig = getTriggerBits(iEvent,triggNames_);
   
	double c = 29.9792458; // speed of light in cm/ns
	if(genParHandle.isValid())
	{
		for (auto& genPar : *genParHandle)
		{	
			if(genPar.pdgId() == 1000015 and genPar.status() == 1)
			{
				double momentum =  genPar.p4().P();
				double energy = genPar.p4().E();
				double beta = momentum / energy;
				cout << "Beta: " <<  beta << "trig: " << trig << endl;
				fHbeta_tot->Fill(beta);
				if(trig) fHbeta_pas->Fill(beta);
			}
		}
	}
	
}


// ------------ method called once each job just before starting event loop  ------------
void 
HSCPTrigger::beginJob()
{
	edm::Service<TFileService> fs;
	fHbeta_pas = fs->make<TH1D>("fHbeta_pas", "Beta of triggered sTau", 50, 0., 1. );
	fHbeta_tot = fs->make<TH1D>("fHbeta_tot", "Beta of tot sTau", 50, 0., 1. );
	return;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HSCPTrigger::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HSCPTrigger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HSCPTrigger);
