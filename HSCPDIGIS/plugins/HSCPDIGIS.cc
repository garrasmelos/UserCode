// -*- C++ -*-
//
// Package:    test/HSCPDIGIS
// Class:      HSCPDIGIS
// 
/**\class HSCPDIGIS HSCPDIGIS.cc test/HSCPDIGIS/plugins/HSCPDIGIS.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Gabriel Ramirez Sanchez
//         Created:  Tue, 14 Feb 2017 12:25:08 GMT
//
//


// system include files
#include <memory>
#include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"



#include "SimDataFormats/RPCDigiSimLink/interface/RPCDigiSimLink.h"
#include "DataFormats/Common/interface/DetSetVector.h"

//
// class declaration
//

class HSCPDIGIS : public edm::EDAnalyzer {
   public:
      explicit HSCPDIGIS(const edm::ParameterSet&);
      ~HSCPDIGIS();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::DetSetVector<RPCDigiSimLink>> digisToken_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
using namespace edm;
using namespace std;
//
// constructors and destructor
//
HSCPDIGIS::HSCPDIGIS(const edm::ParameterSet& iConfig)
: digisToken_(consumes<edm::DetSetVector<RPCDigiSimLink>>(iConfig.getParameter<edm::InputTag>("digisLabel")))
{
   //now do what ever initialization is needed

}


HSCPDIGIS::~HSCPDIGIS()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HSCPDIGIS::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // using namespace edm;
   Handle<edm::DetSetVector<RPCDigiSimLink>> digislink;
   iEvent.getByToken(digisToken_, digislink);

   for(edm::DetSetVector<RPCDigiSimLink>::const_iterator itlink = digislink->begin();itlink!=digislink->end();itlink++)
   {
      for(edm::DetSet<RPCDigiSimLink>::const_iterator itdigi= itlink->data.begin(); itdigi != itlink->data.end();itdigi++)
      {
         int particleId =itdigi->getParticleType();
         if(particleId == 1000015)std::cout << "Particle ID: " << particleId << std::endl;
      }
   }
   

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
HSCPDIGIS::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HSCPDIGIS::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
HSCPDIGIS::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
HSCPDIGIS::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
HSCPDIGIS::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
HSCPDIGIS::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HSCPDIGIS::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HSCPDIGIS);
