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

#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "SimDataFormats/RPCDigiSimLink/interface/RPCDigiSimLink.h"
#include "DataFormats/Common/interface/DetSetVector.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include <Geometry/RPCGeometry/interface/RPCRoll.h>
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/CommonDetUnit/interface/GeomDet.h>

#include "TMath.h"
#include "TH1.h"

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
      TH1D* fHisttofBarrel1_in;
      TH1D* fHisttofBarrel1_out;
      TH1D* fHisttofBarrel2_in;
      TH1D* fHisttofBarrel2_out;
      TH1D* fHisttofBarrel3;
      TH1D* fHisttofBarrel4;
      TH1D* fHisttofEndCapF1;
      TH1D* fHisttofEndCapF2;
      TH1D* fHisttofEndCapF3;
      TH1D* fHisttofEndCapF4;
      TH1D* fHisttofEndCapB1;
      TH1D* fHisttofEndCapB2;
      TH1D* fHisttofEndCapB3;
      TH1D* fHisttofEndCapB4; 
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
: fHisttofBarrel1_in(0), fHisttofBarrel1_out(0), fHisttofBarrel2_in(0), fHisttofBarrel2_out(0), fHisttofBarrel3(0), fHisttofBarrel4(0), fHisttofEndCapF1(0), fHisttofEndCapF2(0), fHisttofEndCapF3(0), fHisttofEndCapF4(0), fHisttofEndCapB1(0), fHisttofEndCapB2(0), fHisttofEndCapB3(0), fHisttofEndCapB4(0),
   digisToken_(consumes<edm::DetSetVector<RPCDigiSimLink>>(iConfig.getParameter<edm::InputTag>("digisLabel")))
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
   //edm::ESHandle<RPCGeometry> rpcGeo;
   //iSetup.get<MuonGeometryRecord>().get(rpcGeo);
   double tof=0;
   for(edm::DetSetVector<RPCDigiSimLink>::const_iterator itlink = digislink->begin();itlink!=digislink->end();itlink++)
   {
      for(edm::DetSet<RPCDigiSimLink>::const_iterator itdigi= itlink->data.begin(); itdigi != itlink->data.end();itdigi++)
      {
         int particleId =itdigi->getParticleType();
         if(TMath::Abs(particleId) == 1000015)
         {
            DetId theDetId = itdigi->getDetUnitId();
            RPCDetId rollId(theDetId);
            tof= itdigi->getTimeOfFlight();
            //const rpcRoll* roll = rpcGeo->roll(rpcId);
            if(rollId.station()==1 && rollId.region()==0 && rollId.layer()==1) fHisttofBarrel1_in->Fill(tof);
            if(rollId.station()==1 && rollId.region()==0 && rollId.layer()==2) fHisttofBarrel1_out->Fill(tof);
            if(rollId.station()==2 && rollId.region()==0 && rollId.layer()==1) fHisttofBarrel2_in->Fill(tof);
            if(rollId.station()==2 && rollId.region()==0 && rollId.layer()==2) fHisttofBarrel2_out->Fill(tof);
            if(rollId.station()==3 && rollId.region()==0) fHisttofBarrel3->Fill(tof);
            if(rollId.station()==4 && rollId.region()==0) fHisttofBarrel4->Fill(tof);
            if(rollId.station()==1 && rollId.region()==1) fHisttofEndCapF1->Fill(tof);
            if(rollId.station()==2 && rollId.region()==1) fHisttofEndCapF2->Fill(tof);
            if(rollId.station()==3 && rollId.region()==1) fHisttofEndCapF3->Fill(tof);
            if(rollId.station()==4 && rollId.region()==1) fHisttofEndCapF4->Fill(tof);
            if(rollId.station()==1 && rollId.region()==-1) fHisttofEndCapB1->Fill(tof);
            if(rollId.station()==2 && rollId.region()==-1) fHisttofEndCapB2->Fill(tof);
            if(rollId.station()==3 && rollId.region()==-1) fHisttofEndCapB3->Fill(tof);
            if(rollId.station()==4 && rollId.region()==-1) fHisttofEndCapB4->Fill(tof);
         }
      }
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
HSCPDIGIS::beginJob()
{
   Service<TFileService> fs;
   fHisttofBarrel1_in = fs->make<TH1D>("HistotofBarrelStation1_in", "toF Barrel station 1 in (DIGIS)",150,1.,150.);
   fHisttofBarrel1_out = fs->make<TH1D>("HistotofBarrelStation1_out","tof Barrel Station 1 out (DIGIS)",150,0.,150.);
   fHisttofBarrel2_in = fs->make<TH1D>("HistotofBarrelStation2_in","tof Barrel Station 2 in (DIGIS)",150,0.,150.);
   fHisttofBarrel2_out = fs->make<TH1D>("HistotofBarrelStation2_out","tof Barrel Station 2 out (DIGIS)",150,0.,150.);
   fHisttofBarrel3 = fs->make<TH1D>("HistotofBarrelStation3","tof Barrel Station 3 (DIGIS)",150,0.,150.);
   fHisttofBarrel4 = fs->make<TH1D>("HistotofBarrelStation4","tof Barrel Station 4 (DIGIS)",150,0.,150.);
   fHisttofEndCapF1 = fs->make<TH1D>("HistotofEndCapStationF1","tof EndCap Station 1(Forward) (DIGIS)",150,0.,150.);
   fHisttofEndCapF2 = fs->make<TH1D>("HistotofEndCapStationF2","tof EndCap Station 2(Forward) (DIGIS)",150,0.,150.);
   fHisttofEndCapF3 = fs->make<TH1D>("HistotofEndCapStationF3","tof EndCap Station 3(Forward) (DIGIS)",150,0.,150.);
   fHisttofEndCapF4 = fs->make<TH1D>("HistotofEndCapStationF4","tof EndCap Station 4(Forward) (DIGIS)",150,0.,150.);
   fHisttofEndCapB1 = fs->make<TH1D>("HistotofEndCapStationB1","tof EndCap Station 1(Backward) (DIGIS)",150,0.,150.);
   fHisttofEndCapB2 = fs->make<TH1D>("HistotofEndCapStationB2","tof EndCap Station 2(Backward) (DIGIS)",150,0.,150.);
   fHisttofEndCapB3 = fs->make<TH1D>("HistotofEndCapStationB3","tof EndCap Station 3(Backward) (DIGIS)",150,0.,150.);
   fHisttofEndCapB4 = fs->make<TH1D>("HistotofEndCapStationB4","tof EndCap Station 4(Backward) (DIGIS)",150,0.,150.);
   return;
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
