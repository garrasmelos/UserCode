// -*- C++ -*-
//
// Package:    test/HSCPRecHits
// Class:      HSCPRecHits
// 
/**\class HSCPRecHits HSCPRecHits.cc test/HSCPRecHits/plugins/HSCPRecHits.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Gabriel Ramirez Sanchez
//         Created:  Fri, 10 Mar 2017 16:04:25 GMT
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
#include <vector>
#include "TVector3.h"
//
// class declaration
//

class HSCPRecHits : public edm::EDAnalyzer {
   public:
      explicit HSCPRecHits(const edm::ParameterSet&);
      ~HSCPRecHits();

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

      TTree* rechitTree;
      Int_t bunchX;
      UInt_t stationhit;
      UInt_t layerhit;
      edm::EDGetTokenT<edm::DetSetVector<RPCDigiSimLink>> digisToken_;
      edm::EDGetTokenT<RPCRecHitCollection> recHitToken_;
      
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
HSCPRecHits::HSCPRecHits(const edm::ParameterSet& iConfig)
:  digisToken_(consumes<edm::DetSetVector<RPCDigiSimLink>>(iConfig.getParameter<edm::InputTag>("digisLabel"))),
	recHitToken_(consumes<RPCRecHitCollection>(iConfig.getParameter<edm::InputTag>("recHitLabel")))
{
   //now do what ever initialization is needed

}


HSCPRecHits::~HSCPRecHits()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HSCPRecHits::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   Handle<edm::DetSetVector<RPCDigiSimLink>> digislink;
   iEvent.getByToken(digisToken_, digislink);
   
   Handle<RPCRecHitCollection> rechits;
   iEvent.getByToken(recHitToken_,rechits);
   
   edm::ESHandle<RPCGeometry> rpcGeo;
   iSetup.get<MuonGeometryRecord>().get(rpcGeo);
   
   vector<vector<TVector3>> aPos; //array of positions.
   vector<vector<unsigned int>> aStation;
   vector<vector<unsigned int>> aLayer;
   vector<vector<int>> aRegion;
   vector<vector<int>> aBx;
   unsigned int sStation;
   unsigned int sLayer;
   int sRegion;
   int sBx; 
   
   vector<RPCDetId> detIds;
   for(edm::DetSetVector<RPCDigiSimLink>::const_iterator itlink = digislink->begin();itlink!=digislink->end();itlink++)
   {
      for(edm::DetSet<RPCDigiSimLink>::const_iterator itdigi= itlink->data.begin(); itdigi != itlink->data.end();itdigi++)
      {
         int particleId =itdigi->getParticleType();
         if(TMath::Abs(particleId) == 1000015 && itdigi->getTrackId()==1)
         {
         	DetId theDetId = itdigi->getDetUnitId();
            RPCDetId rollId(theDetId);
            detIds.push_back(rollId);
            cout << rollId << " is a HSCP." << endl;
            
         }
      }
   }
   cout << "We have " << detIds.size() << " HSCP digis" << endl;
   
   for(RPCRecHitCollection::const_iterator rechit_it = rechits->begin(); rechit_it != rechits->end() ; rechit_it++)
   {
   	//cout << "Bx (RecHit): " << rechit_it->BunchX() << endl;
   	sBx = rechit_it->BunchX();
   	
   	//cout << "RPCId: " << rechit_it->rpcId() << endl ;
   	RPCDetId idRoll(rechit_it->rpcId());
   	bool isHSCP= false;
   	for(unsigned int i=0; i<detIds.size();i++)
   	{
   		if(detIds[i]==rechit_it->rpcId()) isHSCP = true;
   	}
   	if(isHSCP)
   	{
   	 	cout << idRoll << " is a HSCP." << endl;
   	}else
   	{
   	 	cout << idRoll << " is not a HSCP." << endl;
   	}
   	
   	
   	
   	sStation = idRoll.station();
      	sLayer = idRoll.layer();
   	sRegion = idRoll.region();
      	//cout << "Global position: " << rechit_it->globalPosition().x() << endl;
   	LocalPoint lPos = rechit_it->localPosition();
   	const RPCRoll* roll = rpcGeo->roll(idRoll);
   	const BoundPlane& rollSurface = roll->surface();
   	GlobalPoint gPos = rollSurface.toGlobal(lPos);
   	
   	//cout << "Global Point: " << gPos.x() << " " << gPos.y() << " " << gPos.z() << endl;
   	TVector3 pos(gPos.x(),gPos.x(),gPos.z());
   	double phi1 = pos.Phi();
   	double theta1 = pos.Theta();
   	bool found = false;
   	double dRmax = 1.;
   	if(aPos.size()==0)
   	{
   		vector<TVector3> vPos;
   		vPos.push_back(pos);
   		aPos.push_back(vPos);
   		
   		vector<unsigned int> vStation;
   		vStation.push_back(sStation);
   		aStation.push_back(vStation);
   		
   		vector<unsigned int> vLayer;
   		vLayer.push_back(sLayer);
   		aLayer.push_back(vLayer);

         	vector<int> vRegion;
         	vRegion.push_back(sRegion);
         	aRegion.push_back(vRegion);
   		
   		vector<int> vBx;
   		vBx.push_back(sBx);
   		aBx.push_back(vBx);
   	}else
   	{
   		int ntks = aPos.size();
   		for(int i=0; i < ntks ; i++)
   		{
   			double dTheta = theta1 - aPos[i][0].Theta();
   			double dPhi = phi1 - aPos[i][0].Phi();
   			double dR= TMath::Sqrt(dTheta*dTheta+dPhi*dPhi);
   			if (dR < dRmax)
   			{
   				aStation[i].push_back(sStation);
   				aLayer[i].push_back(sLayer);
               			aRegion[i].push_back(sRegion);
               			aBx[i].push_back(sBx);
   				aPos[i].push_back(pos);
   				
               			found= true;
   			}
   		}
   		if(!found)
   		{
   			vector<TVector3> vPos;
   			vPos.push_back(pos);
   			aPos.push_back(vPos); 
   			
   			vector<unsigned int> vStation;
   			vStation.push_back(sStation);
   			aStation.push_back(vStation);
   			
   			vector<unsigned int> vLayer;
   			vLayer.push_back(sLayer);
   			aLayer.push_back(vLayer);

            		vector<int> vRegion;
            		vRegion.push_back(sRegion);
            		aRegion.push_back(vRegion);
   		
   			vector<int> vBx;
   			vBx.push_back(sBx);
   			aBx.push_back(vBx);  			
   		}	
   	}
   }
   
   //cout<< "Number of tracks found: " << aPos.size() << endl;
   for(unsigned int j=0; j< aPos.size();j++)
   {	
      	int nhits=aPos[j].size();
   	//cout << "Track " << j << " has " <<nhits << " hits." << endl;
   	if(nhits > 2 && aBx[j][0]>0 &&aRegion[j][0] == 0)
   	{
   		bool timeCorr= true;
         /*for(unsigned int l=1 ; l < aBx[j].size();l++)
         {
            if(aBx[j][l]<aBx[j][l-1]) timeCorr=false;  
         }
         */
         	if(timeCorr){
         	for(int k=0; k<nhits;k++)
   		{
   			bunchX= aBx[j][k];
            		stationhit = aStation[j][k];
   			if(aRegion[j][k]==0)
            		{
               			if(aStation[j][k] == 1 && aLayer[j][k]==1 )stationhit = 1;
               			if(aStation[j][k] == 1 && aLayer[j][k]==2 )stationhit = 2;
               			if(aStation[j][k] == 2 && aLayer[j][k]==1 )stationhit = 3;
               			if(aStation[j][k] == 2 && aLayer[j][k]==2 )stationhit = 4;
               			if(aStation[j][k] == 3 && aLayer[j][k]==1 )stationhit = 5;
               			if(aStation[j][k] == 4 && aLayer[j][k]==1 )stationhit = 6;
           		}
   			layerhit = aLayer[j][k];
			rechitTree->Fill();
   		}
        	}
   	}
   }


}


// ------------ method called once each job just before starting event loop  ------------
void 
HSCPRecHits::beginJob()
{
   Service<TFileService> fs;
   rechitTree = fs->make<TTree>("rechitTree","Tree of Bx in rechit collection");
   rechitTree->Branch("bunchX",		&bunchX,	"bunchX/I");
   rechitTree->Branch("stationhit",	&stationhit,	"stationhit/i");
   rechitTree->Branch("layerhit",	&layerhit,	"layerhit/i");
   return;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HSCPRecHits::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
HSCPRecHits::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
HSCPRecHits::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
HSCPRecHits::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
HSCPRecHits::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HSCPRecHits::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HSCPRecHits);
