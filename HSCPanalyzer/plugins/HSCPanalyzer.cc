// -*- C++ -*-
//
// Package:    test/HSCPanalyzer
// Class:      HSCPanalyzer
// 
/**\class HSCPanalyzer HSCPanalyzer.cc test/HSCPanalyzer/plugins/HSCPanalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Gabriel Ramirez Sanchez
//         Created:  Thu, 09 Feb 2017 14:57:36 GMT
//
//


// system include files
#include <memory>

// user include files
#include <iostream>

//#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "FastSimulation/Tracking/test/FastTrackAnalyzer.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include <DataFormats/RPCRecHit/interface/RPCRecHit.h>
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/RPCDigi/interface/RPCDigi.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include <DataFormats/GeometrySurface/interface/LocalError.h>
#include <DataFormats/GeometryVector/interface/LocalPoint.h>
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "DataFormats/DetId/interface/DetId.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include <Geometry/RPCGeometry/interface/RPCRoll.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/RPCGeometry/interface/RPCGeomServ.h>
#include <Geometry/CommonDetUnit/interface/GeomDet.h>
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"


#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
//
// class declaration
//

class HSCPanalyzer : public edm::EDAnalyzer {
   public:
      explicit HSCPanalyzer(const edm::ParameterSet&);
      ~HSCPanalyzer();
      edm::ESHandle <RPCGeometry> rpcGeo;
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      virtual void beginRun(const edm::Run&, const edm::EventSetup&);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      TH1D* fHistSTauMass;
      TH1D* fHistSTauEta;
      TH1D* fHistSTauBarEta;
      TH1D* fHistSTauBeta;
      TH1D* fHistSTauPhi;
      TH1D* fHistAngle;
      TH1D* fHistDeltaPhi;
      TH2D* fHistSTauEtaBeta;
      TH2D* fHistEta;
      TH1D* fHistdtof;
      TH1D* fHistdtofHits;  
      TH1D* fHistAcc;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
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
HSCPanalyzer::HSCPanalyzer(const edm::ParameterSet& iConfig)
  : fHistSTauMass(0), fHistSTauEta(0),fHistSTauBarEta(0),  fHistSTauBeta(0), fHistSTauPhi(0), fHistAngle(0), fHistDeltaPhi(0), fHistSTauEtaBeta(0), fHistEta(0), fHistdtof(0), fHistdtofHits(0), fHistAcc(0)
{
   //now do what ever initialization is needed

}


HSCPanalyzer::~HSCPanalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HSCPanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   Handle< GenEventInfoProduct > GetInfoHandle;
   iEvent.getByLabel( "generator", GetInfoHandle);

   vector< Handle<PSimHitContainer> > theSimHitContainers;
   iEvent.getManyByType(theSimHitContainers);
   cout <<"The number of SimHits in this event are: " << theSimHitContainers.size() << endl;


   //double qScale = GetInfoHandle->qScale();
   //cout << "qScale = " << qScale << endl;

   Handle< HepMCProduct > EvtHandle;
   iEvent.getByLabel("generator", EvtHandle);
   const HepMC::GenEvent* Evt = EvtHandle->GetEvent();

   HepMC::GenParticle* sTau=0;
   HepMC::GenParticle* sTauBar=0;
   TLorentzVector sTau_p4;
   HepMC::ThreeVector sTau_p3;
   TLorentzVector sTauBar_p4;
   double dtof=0;
   double sTauP=0;
   double sTauMass=0;
   double sTauBeta=0;
   double angle=0;
   double deltaPhi=0;
   const double pi=3.1415926535;
   const double d=3.5;
   const double c=0.0003; //light speed in m/ps

   for (HepMC::GenEvent::vertex_const_iterator vit = Evt->vertices_begin(); vit != Evt->vertices_end(); vit++){
      for(HepMC::GenVertex::particles_out_const_iterator pout=(*vit)->particles_out_const_begin(); pout!=(*vit)->particles_out_const_end(); pout++){
         if( (*pout)->pdg_id() == 1000015 ){
       sTau = (*pout);
    }
    if( (*pout)->pdg_id() == -1000015 ){
       sTauBar = (*pout);
    }
    if(sTauBar!=0 && sTau!=0){
            sTau_p4.SetPxPyPzE(sTau->momentum().px(),sTau->momentum().py(),sTau->momentum().pz(),sTau->momentum().e());
            sTau_p3 = HepMC::ThreeVector((sTau->momentum().px()),(sTau->momentum().py()),(sTau->momentum().pz()));
            sTauP = sTau_p3.r();
            sTauMass= sTau_p4.M();
            sTauBeta= sqrt(sTauP*sTauP/(sTauP*sTauP+sTauMass*sTauMass));
            dtof =(d/c)*(1/sTauBeta-1);
       fHistdtof->Fill(dtof);
            fHistSTauMass->Fill(sTauMass);
            fHistSTauEta->Fill(sTau_p4.Eta());
            fHistSTauEtaBeta->Fill(sTauBeta,abs(sTau_p4.Eta()));
       fHistSTauPhi->Fill(sTau_p4.Phi());
            fHistSTauBeta->Fill(sTauBeta);
            if(abs(sTau_p4.Eta())>1.6){
          fHistAcc->Fill(1.);
          if(dtof>50) fHistAcc->Fill(2.);
               if(dtof>100) fHistAcc->Fill(3.);
               if(dtof>200) fHistAcc->Fill(4.);
               if(dtof>1000) fHistAcc->Fill(5.);
               if(dtof>2000) fHistAcc->Fill(6.);
               if(dtof>4000) fHistAcc->Fill(7.);
               if(dtof>10000) fHistAcc->Fill(8.);
               if(dtof>25000) fHistAcc->Fill(9.);
            }
       cout << sTau->pdg_id()<<"\t"<<sTau_p4.Phi()<<"\t"<<sTau_p4.Eta() << endl;

            sTauBar_p4.SetPxPyPzE(sTauBar->momentum().px(),sTauBar->momentum().py(),sTauBar->momentum().pz(),sTauBar->momentum().e());
            fHistSTauBarEta->Fill(sTauBar_p4.Eta());
            fHistEta->Fill(sTau_p4.Eta(),sTauBar_p4.Eta());

       angle = sTauBar_p4.Angle(sTau_p4.Vect());
       fHistAngle->Fill(angle);
       if(sTau_p4.Phi()<0){
          if(sTauBar_p4.Phi()<0)deltaPhi = abs(sTau_p4.Phi()-sTauBar_p4.Phi());
          else deltaPhi = abs(sTau_p4.Phi()+2*pi-sTauBar_p4.Phi());
              } else {
               if(sTauBar_p4.Phi()<0)deltaPhi = abs(sTauBar_p4.Phi()+2*pi-sTau_p4.Phi());
               else deltaPhi = abs(sTau_p4.Phi()-sTauBar_p4.Phi());
        }
        fHistDeltaPhi->Fill(deltaPhi);
        break;
    }
      }
      if (sTauBar!=0 && sTau!=0){
         break;
      }

   }

   vector<PSimHit> theSimHits;
   for(int i=0; i < int(theSimHitContainers.size());i++){
      theSimHits.insert(theSimHits.end(),theSimHitContainers.at(i)->begin(), theSimHitContainers.at(i)->end());
   }

   for(vector<PSimHit>::const_iterator iHit = theSimHits.begin(); iHit != theSimHits.end(); iHit++){
   //   if(abs((*iHit).particleType())== 13) cout << "Muon here!"<< endl;
   //}
      if((int)abs((*iHit).particleType())== 13) {
         DetId theDetUnitId = DetId((*iHit).detUnitId());
         DetId simdetid= DetId((*iHit).detUnitId());
         if(simdetid.det()==DetId::Muon) cout << "Muon detector" << endl;
         if(simdetid.det()==DetId::Muon &&  simdetid.subdetId()== MuonSubdetId::RPC){

            RPCDetId rollId(theDetUnitId);
            RPCGeomServ rpcsrv(rollId);

            cout << "here" << endl;

            const RPCRoll * rollasociated = rpcGeo->roll(rollId);
            const BoundPlane & RPCSurface = rollasociated->surface();

            //GlobalPoint SimHitInGlobal = RPCSurface.toGlobal((*iHit).localPosition());
           //if(rollId.region()==0) continue; //skip barrel
           //if(rollId.ring()!=1) continue; //leave only ring
            fHistdtofHits->Fill((*iHit).timeOfFlight());
            if(rollId.station()==1) {cout << "hit in station 1. tof is: " << (*iHit).timeOfFlight()<< endl; }
         }
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
  HSCPanalyzer::beginRun(const edm::Run& run, const edm::EventSetup& iSetup)
 {
      iSetup.get<MuonGeometryRecord>().get(rpcGeo);

 }
// ------------ method called once each job just before starting event loop  ------------
void 
HSCPanalyzer::beginJob()
{
   Service<TFileService> fs;
   fHistSTauMass = fs->make<TH1D>("HistSTauMass", "sTau inv. mass", 100, 60., 3000. );
   fHistSTauEta = fs->make<TH1D>("HistSTauEta", "sTau Eta", 150, -4., 4.);
   fHistSTauEtaBeta = fs->make<TH2D>("HistSTauEtaBeta", "sTau Eta vs Beta", 150, 0., 1.,150,0.,4.);
   fHistEta = fs->make<TH2D>("HistEta","sTau Eta vs sTauBar Eta",150, -4., 4.,150, -4., 4.);
   fHistSTauBarEta = fs->make<TH1D>("HistSTauBarEta","sTau bar Eta",150,-4.,4.);
   fHistSTauBeta = fs->make<TH1D>("HistSTauBeta", "sTau #Beta", 150, 0, 1.0);
   fHistSTauPhi = fs->make<TH1D>("HistSTauPhi","sTau #phi", 150, -3.15,3.15);
   fHistAngle = fs->make<TH1D>("HistAngle", "sTau-sTau bar angle",150, 0.,3.15);
   fHistDeltaPhi = fs->make<TH1D>("HistDeltaPhi","abs(Phi(sTau)-Phi(sTauBar))",150,0.,3.15);
   fHistdtof = fs->make<TH1D>("Histdtof","dtof",150, 0.,40000.);
   fHistdtofHits = fs->make<TH1D>("HistdtofHits","dtofHits",150,0.0,400.);
   fHistAcc = fs->make<TH1D>("HistoAcc","Acceptance ring 1",10,0,10);
   return ;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HSCPanalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------

// ------------ method called when ending the processing of a run  ------------
/*
void 
HSCPanalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
HSCPanalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
HSCPanalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HSCPanalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HSCPanalyzer);
