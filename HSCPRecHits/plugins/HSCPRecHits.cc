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

#define  Nhltpaths 3


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
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
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

class HSCPRecHits : public edm::EDAnalyzer {
   public:
      explicit HSCPRecHits(const edm::ParameterSet&);
      ~HSCPRecHits();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      std::string arrayHLTpathsNames[Nhltpaths];
      bool filled;

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
      TTree* effTree;
      
      Int_t bunchX;
      UInt_t stationhit;
      UInt_t layerhit;
      UInt_t isChosen;
      
      UInt_t triggL1;
      Float_t beta;
      Float_t betaSim;
		Float_t betaGen;
		Float_t etaGen;     
      Float_t timeOfFlight;
      
		edm::EDGetTokenT<GenEventInfoProduct> genEventInfo_;
		edm::EDGetTokenT<edm::HepMCProduct> hepEventInfo_;
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
:  genEventInfo_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo"))),
	hepEventInfo_(consumes<edm::HepMCProduct>(iConfig.getParameter<edm::InputTag>("hepEventInfo"))),
	digisToken_(consumes<edm::DetSetVector<RPCDigiSimLink>>(iConfig.getParameter<edm::InputTag>("digisLabel"))),
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
   
   Handle<GenEventInfoProduct> GetInfoHandle;
   iEvent.getByToken(genEventInfo_, GetInfoHandle); 
   
   Handle< HepMCProduct > EvtHandle;
   iEvent.getByToken(hepEventInfo_, EvtHandle);
   const HepMC::GenEvent* Evt = EvtHandle->GetEvent();

   edm::Handle<edm::TriggerResults> trigResults; //our trigger result object
   edm::InputTag trigResultsTag("TriggerResults","","DIGI2RAW"); //make sure have correct process on MC

   // Access Trigger Results
   edm::Handle<edm::TriggerResults> hltTriggerResultHandle;
   iEvent.getByLabel(trigResultsTag,hltTriggerResultHandle);

   //cout<<"Number of paths in Trigger:"<<hltTriggerResultHandle->size()<<endl;

   if(filled==false){
       const edm::TriggerNames & triggerNames = iEvent.triggerNames(*hltTriggerResultHandle);
       //std::cout<<triggerNames.triggerName(1)<<std::endl; 
       //std::cout<<arrayHLTpathsNames.size()<<std::endl;
       for(int i=0;i<Nhltpaths;i++){
	   arrayHLTpathsNames[i]=triggerNames.triggerName(i);
       }
       filled=true;
       //std::cout<<arrayHLTpathsNames.size()<<std::endl;
   }

   


   int hltCount= hltTriggerResultHandle->size();
   bool allHLTResults[Nhltpaths] = { false };
   int bin;
   triggL1 = 0;
   if(!hltTriggerResultHandle.isValid()) {
       //std::cout << "invalid handle for HLT TriggerResults" << std::endl;
   }else{
       for(bin = 0 ; bin < hltCount ; bin++){
	   allHLTResults[bin] = hltTriggerResultHandle->accept(bin);
	   //std::cout<<"bit:"<<i<<" "<<allHLTResults[i]<<std::endl;
       }
       triggL1 = hltTriggerResultHandle->accept(1);
   }

   for(bin = 0 ; bin < hltCount ; bin++){
       if(allHLTResults[bin]) cout<<"bin:"<<bin<<" arrayHLTpathsNames:"<<arrayHLTpathsNames[bin]<<endl;
   }



   vector<vector<TVector3>> aPos; //array of positions.
   vector<vector<unsigned int>> aStation;
   vector<vector<unsigned int>> aLayer;
   vector<vector<int>> aRegion;
   vector<vector<int>> aBx;
   vector<vector<double>> aTof;
   
   
   unsigned int sStation;
   unsigned int sLayer;
   int sRegion;
   int sBx; 
   int sTof;
  
   
   double c = 29.979;
	
	vector<double> tof;
   vector<RPCDetId> detIds;
   
   HepMC::GenParticle* sTau=0;
   
   for (HepMC::GenEvent::vertex_const_iterator vit = Evt->vertices_begin(); vit != Evt->vertices_end(); vit++)
   {
    	for(HepMC::GenVertex::particles_out_const_iterator pout=(*vit)->particles_out_const_begin(); pout!=(*vit)->particles_out_const_end(); pout++)
    	{
      	if( (*pout)->pdg_id() == 1000015 &&(*pout)->status()==1)
      	{
      		sTau = (*pout);
      		TLorentzVector sTau_p4;
      		sTau_p4.SetPxPyPzE(sTau->momentum().px(),sTau->momentum().py(),sTau->momentum().pz(),sTau->momentum().e());
      		double sTauP = sTau_p4.P();
				double sTauMass= sTau_p4.M();
				betaGen= sqrt(sTauP*sTauP/(sTauP*sTauP+sTauMass*sTauMass));
				etaGen = sTau_p4.Eta();
				
      	}
      }
   }
   for(edm::DetSetVector<RPCDigiSimLink>::const_iterator itlink = digislink->begin();itlink!=digislink->end();itlink++)
   {
      for(edm::DetSet<RPCDigiSimLink>::const_iterator itdigi= itlink->data.begin(); itdigi != itlink->data.end();itdigi++)
      {
         int particleId =itdigi->getParticleType();
         //cout << "Process type: "<< itdigi->getProcessType() << endl;
         //if(TMath::Abs(particleId) == 1000015 && itdigi->getTrackId()==1)
         if(particleId == 1000015 )
         {
         	DetId theDetId = itdigi->getDetUnitId();
            RPCDetId rollId(theDetId);
            detIds.push_back(rollId);
            tof.push_back(itdigi->getTimeOfFlight());
            //cout<< "Track ID/particleID :" << itdigi->getTrackId() << "/" << itdigi->getParticleType() <<  endl;
            //cout << rollId << " is a HSCP." << endl;
            
         }
      }
   }
   //cout << "We have " << detIds.size() << " HSCP digis" << endl;
   
   for(RPCRecHitCollection::const_iterator rechit_it = rechits->begin(); rechit_it != rechits->end() ; rechit_it++)
   {
   	//cout << "Bx (RecHit): " << rechit_it->BunchX() << endl;
   	sBx = rechit_it->BunchX();
   	
   	
   	//cout << "RPCId: " << rechit_it->rpcId() << endl ;
   	RPCDetId idRoll(rechit_it->rpcId());
   	bool isHSCP= false;
   	for(unsigned int i=0; i<detIds.size();i++)
   	{
   		if(detIds[i]==rechit_it->rpcId()) 
   		{
   			sTof= tof[i];
   			isHSCP = true;
   			break;
   		}
   	}
   	if(!isHSCP)
   	{
   		continue;
   	}

   	sStation = idRoll.station();
      sLayer = idRoll.layer();
   	sRegion = idRoll.region();
      	//cout << "Global position: " << rechit_it->globalPosition().x() << endl;
   	LocalPoint lPos = rechit_it->localPosition();
   	const RPCRoll* roll = rpcGeo->roll(idRoll);
   	const BoundPlane& rollSurface = roll->surface();
   	GlobalPoint gPos = rollSurface.toGlobal(lPos);
   	TVector3 pos(gPos.x(),gPos.x(),gPos.z());
   	//cout << "Global Point: " << gPos.x() << " " << gPos.y() << " " << gPos.z() << endl;
   	//cout << "Distance: " << pos.Mag() << endl;
   	
   	double phi1 = pos.Phi();
   	double theta1 = pos.Theta();
   	bool found = false;
   	double dRmax = 0.2;
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
   		
   		vector<double> vTof;
   		vTof.push_back(sTof);
   		aTof.push_back(vTof);
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
            	aTof[i].push_back(sTof);
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
   			
   			vector<double> vTof;
   			vTof.push_back(sTof);
   			aTof.push_back(vTof);			
   		}	
   	}
   }
   
   cout<< "Number of tracks found: " << aPos.size() << endl;
   int theTrack=-1;
   unsigned int theTrackSize=0;
   for(unsigned int j=0; j< aPos.size();j++)
   {	
   	if (aPos[j].size() > theTrackSize)
   	{
   	 	theTrackSize = aPos[j].size();
   	 	theTrack = j;
   	}
   }
   if(theTrack!=-1)
   {	
      int nhits=aPos[theTrack].size();
      cout << "nhits: " << nhits << endl;
      double betaSum = 0.;
      double betaSimSum = 0.;
   	//if(nhits > 2 &&aRegion[theTrack][0] == 0)
   	int isTrack=1;
      for(unsigned int m=0 ; m < aBx[theTrack].size();m++)
      {
         if(nhits<3 && aRegion[theTrack][m]!=0)  isTrack=0;
      }
   	if(isTrack)
   	{
   		isChosen=1;
   		
         for(unsigned int l=0 ; l < aBx[theTrack].size();l++)
         {
            if(aBx[theTrack][l]<6 )  isChosen=0;
         }
         
            for(int ihit=0; ihit<nhits;ihit++)
  		   {
            bunchX= aBx[theTrack][ihit];
            timeOfFlight = aTof[theTrack][ihit];
            stationhit = aStation[theTrack][ihit];
           	Float_t d = aPos[theTrack][ihit].Mag();
           	//cout <<  aPos[theTrack][ihit].Mag() << endl;
           	betaSimSum += d/(timeOfFlight*c);
            betaSum += d/(bunchX*c+d);
            
           	if(aRegion[theTrack][ihit]==0)
           	{
              	if(aStation[theTrack][ihit] == 1 && aLayer[theTrack][ihit]==1 )stationhit = 1;
              	if(aStation[theTrack][ihit] == 1 && aLayer[theTrack][ihit]==2 )stationhit = 2;
              	if(aStation[theTrack][ihit] == 2 && aLayer[theTrack][ihit]==1 )stationhit = 3;
              	if(aStation[theTrack][ihit] == 2 && aLayer[theTrack][ihit]==2 )stationhit = 4;
              	if(aStation[theTrack][ihit] == 3 && aLayer[theTrack][ihit]==1 )stationhit = 5;
              	if(aStation[theTrack][ihit] == 4 && aLayer[theTrack][ihit]==1 )stationhit = 6;
        		}
  				layerhit = aLayer[theTrack][ihit];
				rechitTree->Fill();
  			}
			//cout << "########################## " << betaSum << "########################" << nhits << endl;
        	beta = betaSum/nhits;
        	betaSim = betaSimSum/nhits;
      	effTree->Fill();
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
   rechitTree->Branch("timeOfFlight",		&timeOfFlight,	"timeOfFlight/f");
   rechitTree->Branch("stationhit",	&stationhit,	"stationhit/i");
   rechitTree->Branch("layerhit",	&layerhit,	"layerhit/i");

   effTree = fs->make<TTree>("effTree","effTree");
   effTree->Branch("isChosen",&isChosen, "isChosen/i");
   effTree->Branch("beta",&beta,"beta/f");
   effTree->Branch("betaSim",&betaSim,"betaSim/f");
   effTree->Branch("betaGen",&betaGen,"betaGen/f");
   effTree->Branch("etaGen",&etaGen,"etaGen/f");
   effTree->Branch("triggL1", &triggL1, "triggL1/i");
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
