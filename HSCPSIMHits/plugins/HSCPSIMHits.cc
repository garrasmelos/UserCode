// -*- C++ -*-
//
// Package:    test/HSCPSIMHits
// Class:      HSCPSIMHits
// 
/**\class HSCPSIMHits HSCPSIMHits.cc test/HSCPSIMHits/plugins/HSCPSIMHits.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Gabriel Ramirez Sanchez
//         Created:  Tue, 14 Feb 2017 12:08:09 GMT
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

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "FastSimulation/Tracking/test/FastTrackAnalyzer.h"

#include "DataFormats/Common/interface/Handle.h"
#include <DataFormats/RPCRecHit/interface/RPCRecHit.h>
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include <DataFormats/GeometrySurface/interface/LocalError.h>
#include <DataFormats/GeometryVector/interface/LocalPoint.h>
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include <Geometry/RPCGeometry/interface/RPCRoll.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/RPCGeometry/interface/RPCGeomServ.h>
#include <Geometry/CommonDetUnit/interface/GeomDet.h>
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"



#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
//
// class declaration
//

class HSCPSIMHits : public edm::EDAnalyzer {
   public:
      explicit HSCPSIMHits(const edm::ParameterSet&);
      ~HSCPSIMHits();
      edm::ESHandle <RPCGeometry> rpcGeo;
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      virtual void beginRun(const edm::Run&, const edm::EventSetup&);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
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
      TH1D* fHistdtofBarrel1;
      TH1D* fHistdtofBarrel2;
      TH1D* fHistdtofBarrel3;
      TH1D* fHistdtofBarrel4;
      TH1D* fHistdtofEndCapF1;
      TH1D* fHistdtofEndCapF2;
      TH1D* fHistdtofEndCapF3;
      TH1D* fHistdtofEndCapF4;
      TH1D* fHistdtofEndCapB1;
      TH1D* fHistdtofEndCapB2;
      TH1D* fHistdtofEndCapB3;
      TH1D* fHistdtofEndCapB4;
      edm::EDGetTokenT<std::vector<PSimHit>> hitsToken_;
};

//
// constants, enums and typedefs
//
using namespace edm;
using namespace std;
//
// static data member definitions
//

//
// constructors and destructor
//
HSCPSIMHits::HSCPSIMHits(const edm::ParameterSet& iConfig)
: fHistSTauMass(0), fHistSTauEta(0),fHistSTauBarEta(0),  fHistSTauBeta(0), fHistSTauPhi(0), fHistAngle(0), fHistDeltaPhi(0), fHistSTauEtaBeta(0), fHistEta(0), fHistdtof(0), fHistdtofHits(0), fHistAcc(0), fHistdtofBarrel1(0), fHistdtofBarrel2(0), fHistdtofBarrel3(0), fHistdtofBarrel4(0), fHistdtofEndCapF1(0), fHistdtofEndCapF2(0), fHistdtofEndCapF3(0), fHistdtofEndCapF4(0), fHistdtofEndCapB1(0), fHistdtofEndCapB2(0), fHistdtofEndCapB3(0), fHistdtofEndCapB4(0),
   hitsToken_(consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("hitsLabel")))
{
   //now do what ever initialization is needed

}


HSCPSIMHits::~HSCPSIMHits()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HSCPSIMHits::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   Handle< GenEventInfoProduct > GetInfoHandle;
   iEvent.getByLabel( "generator", GetInfoHandle);

   Handle< HepMCProduct > EvtHandle;
   iEvent.getByLabel("generator", EvtHandle);
   const HepMC::GenEvent* Evt = EvtHandle->GetEvent();

   Handle< std::vector<PSimHit>> hits;
   iEvent.getByToken(hitsToken_, hits);

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
   for(vector<PSimHit>::const_iterator iHit = hits->begin(); iHit != hits->end(); iHit++){
      if((int)abs((*iHit).particleType())== 1000015) {
         DetId theDetUnitId = DetId((*iHit).detUnitId());
         DetId simdetid= DetId((*iHit).detUnitId());

         if(simdetid.det()==DetId::Muon &&  simdetid.subdetId()== MuonSubdetId::RPC){
            RPCDetId rollId(theDetUnitId);
            RPCGeomServ rpcsrv(rollId);
            const RPCRoll * rollasociated = rpcGeo->roll(rollId);
            const BoundPlane & RPCSurface = rollasociated->surface();
            //GlobalPoint SimHitInGlobal = RPCSurface.toGlobal((*iHit).localPosition());
            if(rollId.station()==1 && rollId.region()==1) fHistdtofHits->Fill((*iHit).timeOfFlight());
            if(rollId.station()==1 && rollId.region()==0) fHistdtofBarrel1->Fill((*iHit).timeOfFlight());
            if(rollId.station()==2 && rollId.region()==0) fHistdtofBarrel2->Fill((*iHit).timeOfFlight());
            if(rollId.station()==3 && rollId.region()==0) fHistdtofBarrel3->Fill((*iHit).timeOfFlight());
            if(rollId.station()==4 && rollId.region()==0) fHistdtofBarrel4->Fill((*iHit).timeOfFlight());
            if(rollId.station()==1 && rollId.region()==1) fHistdtofEndCapF1->Fill((*iHit).timeOfFlight());
            if(rollId.station()==2 && rollId.region()==1) fHistdtofEndCapF2->Fill((*iHit).timeOfFlight());
            if(rollId.station()==3 && rollId.region()==1) fHistdtofEndCapF3->Fill((*iHit).timeOfFlight());
            if(rollId.station()==4 && rollId.region()==1) fHistdtofEndCapF4->Fill((*iHit).timeOfFlight());
            if(rollId.station()==1 && rollId.region()==-1) fHistdtofEndCapB1->Fill((*iHit).timeOfFlight());
            if(rollId.station()==2 && rollId.region()==-1) fHistdtofEndCapB2->Fill((*iHit).timeOfFlight());
            if(rollId.station()==3 && rollId.region()==-1) fHistdtofEndCapB3->Fill((*iHit).timeOfFlight());
            if(rollId.station()==4 && rollId.region()==-1) fHistdtofEndCapB4->Fill((*iHit).timeOfFlight());
         }
      }
   }
}
// ------------ method called once each job just before starting event loop  ------------
void
HSCPSIMHits::beginRun(const edm::Run& run, const edm::EventSetup& iSetup)
{
   iSetup.get<MuonGeometryRecord>().get(rpcGeo);
}

// ------------ method called once each job just before starting event loop  ------------
void 
HSCPSIMHits::beginJob()
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
   fHistdtofBarrel1 = fs->make<TH1D>("HistodtofBarrelStation1","Dtof Barrel Station 1",150,0.,400.);
   fHistdtofBarrel2 = fs->make<TH1D>("HistodtofBarrelStation2","Dtof Barrel Station 2",150,0.,400.);
   fHistdtofBarrel3 = fs->make<TH1D>("HistodtofBarrelStation3","Dtof Barrel Station 3",150,0.,400.);
   fHistdtofBarrel4 = fs->make<TH1D>("HistodtofBarrelStation4","Dtof Barrel Station 4",150,0.,400.);
   fHistdtofEndCapF1 = fs->make<TH1D>("HistodtofEndCapStationF1","Dtof EndCap Station 1(Forward)",150,0.,400.);
   fHistdtofEndCapF2 = fs->make<TH1D>("HistodtofEndCapStationF2","Dtof EndCap Station 2(Forward)",150,0.,400.);
   fHistdtofEndCapF3 = fs->make<TH1D>("HistodtofEndCapStationF3","Dtof EndCap Station 3(Forward)",150,0.,400.);
   fHistdtofEndCapF4 = fs->make<TH1D>("HistodtofEndCapStationF4","Dtof EndCap Station 4(Forward)",150,0.,400.);
   fHistdtofEndCapB1 = fs->make<TH1D>("HistodtofEndCapStationB1","Dtof EndCap Station 1(Backward)",150,0.,400.);
   fHistdtofEndCapB2 = fs->make<TH1D>("HistodtofEndCapStationB2","Dtof EndCap Station 2(Backward)",150,0.,400.);
   fHistdtofEndCapB3 = fs->make<TH1D>("HistodtofEndCapStationB3","Dtof EndCap Station 3(Backward)",150,0.,400.);
   fHistdtofEndCapB4 = fs->make<TH1D>("HistodtofEndCapStationB4","Dtof EndCap Station 4(Backward)",150,0.,400.);
   return;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HSCPSIMHits::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
HSCPSIMHits::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
HSCPSIMHits::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
HSCPSIMHits::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
HSCPSIMHits::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HSCPSIMHits::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HSCPSIMHits);
