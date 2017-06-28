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
#include "FWCore/Framework/interface/ESHandle.h"

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
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/RPCDigi/interface/RPCDigi.h"

#include "TMath.h"
#include "TH1.h"
#include "TTree.h"
#include "TVector3.h"

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
      TH1D* fHistDiff2_3Tof;
      TH1D* fHistDiff3_4Tof; 
      TH1D* fHistDiff2_3Bx;
      TH1D* fHistDiff3_4Bx;
      TTree* hscpTree;
      Int_t bx;
      Float_t tof;
      UInt_t station;
      Int_t region;
      UInt_t layer;
      UInt_t roll;
      
      
      edm::EDGetTokenT<edm::DetSetVector<RPCDigiSimLink>> digisToken_;
      edm::EDGetTokenT<RPCDigiCollection> digis2Token_;
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
: fHistDiff2_3Tof(0), fHistDiff3_4Tof(0), fHistDiff2_3Bx(0), fHistDiff3_4Bx(0),
   digisToken_(consumes<edm::DetSetVector<RPCDigiSimLink>>(iConfig.getParameter<edm::InputTag>("digisLabel"))),
   digis2Token_(consumes<RPCDigiCollection>(iConfig.getParameter<edm::InputTag>("digis2Label")))
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
   
   Handle<RPCDigiCollection> digis;
   iEvent.getByToken(digis2Token_, digis);
   
   
   for(RPCDigiCollection::DigiRangeIterator detUnitIt = digis->begin(); detUnitIt != digis->end() ; detUnitIt++)
   {
   	const RPCDetId& id = (*detUnitIt).first;
   	const RPCDigiCollection::Range& range = (*detUnitIt).second;
   	cout << "######## STATION: " << id.station() << endl;
   	
   	for(RPCDigiCollection::const_iterator digi_it = range.first ; digi_it != range.second ; ++digi_it)
   	{
   		
   		cout << "######## BX: " << digi_it->bx() << endl;
   	}
   }

   tof=0;
   bx=0;
   station=0;
   region=0;
   layer=0;
   roll=0;
   if(digislink->size()>0 && digislink.isValid())
   {
   for(edm::DetSetVector<RPCDigiSimLink>::const_iterator itlink = digislink->begin();itlink!=digislink->end();itlink++)
   {
      for(edm::DetSet<RPCDigiSimLink>::const_iterator itdigi= itlink->data.begin(); itdigi != itlink->data.end();itdigi++)
      {
         int particleId =itdigi->getParticleType();
		 //bx = itdigi->getBx();
		 //hscpTree->Fill();
         //int pId = 13; //Muons
         int pId = 1000015; //sTau
         if(TMath::Abs(particleId) == pId && itdigi->getTrackId()==1)
         {
            DetId theDetId = itdigi->getDetUnitId();
            RPCDetId rollId(theDetId);
            tof= itdigi->getTimeOfFlight();
            //const rpcRoll* roll = rpcGeo->roll(rpcId);
            //if(rollId.station()==1 && rollId.region()==0 && rollId.layer()==1) fHisttofBarrel1_in->Fill(tof);
     
            bx = itdigi->getBx();
            station = rollId.station();
            region = rollId.region();
            layer = rollId.layer();
            roll = rollId.roll();
            cout << "FILLING TREE" << endl;
            hscpTree->Fill();
            cout << "TrackId: " << itdigi->getTrackId() << endl;
            if(rollId.station()==3 && rollId.region()==-1 )
            {
            	for(edm::DetSetVector<RPCDigiSimLink>::const_iterator itlink2 = digislink->begin();itlink2!=digislink->end();itlink2++)
            	{
            		for(edm::DetSet<RPCDigiSimLink>::const_iterator itdigi2= itlink2->data.begin(); itdigi2 != itlink2->data.end();itdigi2++)
            		{
            			DetId theDetId2 = itdigi2->getDetUnitId();
            			RPCDetId rollId2(theDetId2);
            			if(rollId2.station()==2 && rollId2.region()==-1 && rollId2.ring()==rollId.ring() && rollId.sector()==rollId2.sector())
            			{            			
            				//fHistDiffTof->Fill(itdigi2->getBx()- itdigi->getBx());
            				fHistDiff2_3Tof->Fill(itdigi->getTimeOfFlight()- itdigi2->getTimeOfFlight());
            				fHistDiff2_3Bx->Fill(itdigi->getBx()- itdigi2->getBx());
            			}  
            			if(rollId2.station()==4 && rollId2.region()==-1 && rollId2.ring()==rollId.ring() && rollId.sector()==rollId2.sector())
            			{
								//fHistDiffTof->Fill(itdigi2->getBx()- itdigi->getBx());
            				fHistDiff3_4Tof->Fill(itdigi2->getTimeOfFlight()- itdigi->getTimeOfFlight());
            				fHistDiff3_4Bx->Fill(itdigi2->getBx()- itdigi->getBx());
            			}          		
            		}
            	}
            }
         }
      }
   }
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
HSCPDIGIS::beginJob()
{
   Service<TFileService> fs;
   fHistDiff2_3Tof = fs->make<TH1D>("fHistDiff2_3Tof","fHistDiff2_3Tof",100,-50.,50.);
   fHistDiff3_4Tof = fs->make<TH1D>("fHistDiff3_4Tof","fHistDiff3_4Tof",100,-50.,50.);
	fHistDiff2_3Bx = fs->make<TH1D>("fHistDiff2_3Bx","fHistDiff2_3Bx",100,-50.,50.);
   fHistDiff3_4Bx = fs->make<TH1D>("fHistDiff3_4Bx","fHistDiff3_4Bx",100,-50.,50.);
   hscpTree = fs->make<TTree>("hscpTree","Tree of HSCP particles");
   hscpTree->Branch("bx",			&bx,			"bx/I");
   hscpTree->Branch("tof",			&tof,			"tof/F");
   hscpTree->Branch("station",	&station,	"station/i");
   hscpTree->Branch("region",		&region,		"region/I");
   hscpTree->Branch("layer",		&layer,		"layer/i");
   hscpTree->Branch("roll",		&roll,		"roll/i");
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
