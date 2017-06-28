// system include files
#include <memory>
#include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

//#include "FastSimulation/Tracking/test/FastTrackAnalyzer.h"

#include "DataFormats/Common/interface/Handle.h"
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
  TH1D* fHistSTauP;
  TH1D* fHistSTauEta;
  TH1D* fHistSTauBarEta;
  TH1D* fHistSTauBeta;
  TH1D* fHistSTauPhi;
  TH1D* fHistAngle;
  TH1D* fHistDeltaPhi;
  TH2D* fHistSTauEtaBeta;
  TH2D* fHistEta;
  TH1D* fHisttof;
  TH1D* fHisttofHits;
  TH1D* fHistAcc;
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
  edm::EDGetTokenT<std::vector<PSimHit>> hitsToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParToken_;
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
HSCPSIMHits::HSCPSIMHits(const edm::ParameterSet& iConfig): 
  fHistSTauMass(0), 
  fHistSTauP(0),
  fHistSTauEta(0),
  fHistSTauBarEta(0),  
  fHistSTauBeta(0), 
  fHistSTauPhi(0), 
  fHistAngle(0), 
  fHistDeltaPhi(0), 
  fHistSTauEtaBeta(0), 
  fHistEta(0), 
  fHisttof(0), 
  fHisttofHits(0), 
  fHistAcc(0), 
  fHisttofBarrel1_in(0), 
  fHisttofBarrel1_out(0), 
  fHisttofBarrel2_in(0), 
  fHisttofBarrel2_out(0), 
  fHisttofBarrel3(0), 
  fHisttofBarrel4(0), 
  fHisttofEndCapF1(0), 
  fHisttofEndCapF2(0), 
  fHisttofEndCapF3(0), 
  fHisttofEndCapF4(0), 
  fHisttofEndCapB1(0), 
  fHisttofEndCapB2(0), 
  fHisttofEndCapB3(0), 
  fHisttofEndCapB4(0),
  hitsToken_(consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("hitsLabel"))),
  genParToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleLabel")))
{ 
}


HSCPSIMHits::~HSCPSIMHits()
{
  
}
void
HSCPSIMHits::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  Handle<reco::GenParticleCollection> genParticleHandle;
  iEvent.getByToken(genParToken_,genParticleHandle); 
      
  Handle< std::vector<PSimHit>> hits;
  iEvent.getByToken(hitsToken_, hits);
  
  TLorentzVector sTau_p4;
  TLorentzVector sTauBar_p4;
  double tof=0;
  double sTauP=0;
  double sTauMass=0;
  double sTauBeta=0;
  double angle=0;
  double deltaPhi=0;
  const double pi=3.1415926535;
  const double d=3.5;
  const double c=0.0003; //light speed in m/ps

  if(genParticleHandle.isValid())
  {
    for(auto& pout : *genParticleHandle)
    {
      if( pout.pdgId() == 1000015 ){
        sTau_p4.SetPxPyPzE(pout.p4().Px(),pout.p4().Py(),pout.p4().Pz(),pout.p4().E());
      }
      if( pout.pdgId() == -1000015 ){
        sTauBar_p4.SetPxPyPzE(pout.p4().Px(),pout.p4().Py(),pout.p4().Pz(),pout.p4().E());
      }
    }
  }
  
  sTauP = sTau_p4.P();
  sTauMass= sTau_p4.M();   
  sTauBeta= sqrt(sTauP*sTauP/(sTauP*sTauP+sTauMass*sTauMass));
  tof =(d/c)*(1/sTauBeta); 
 
  fHisttof->Fill(tof);
  fHistSTauMass->Fill(sTauMass);
  fHistSTauP->Fill(sTauP);
  fHistSTauEta->Fill(sTau_p4.Eta());
  fHistSTauEtaBeta->Fill(sTauBeta,abs(sTau_p4.Eta()));
  fHistSTauPhi->Fill(sTau_p4.Phi());
  fHistSTauBeta->Fill(sTauBeta);
  if(abs(sTau_p4.Eta())>1.6){
    fHistAcc->Fill(1.);
    if(tof>50) fHistAcc->Fill(2.);
    if(tof>100) fHistAcc->Fill(3.);
    if(tof>200) fHistAcc->Fill(4.);
    if(tof>1000) fHistAcc->Fill(5.);
    if(tof>2000) fHistAcc->Fill(6.);
    if(tof>4000) fHistAcc->Fill(7.);
    if(tof>10000) fHistAcc->Fill(8.);
    if(tof>25000) fHistAcc->Fill(9.);
  }
  fHistSTauBarEta->Fill(sTauBar_p4.Eta());
  fHistEta->Fill(sTau_p4.Eta(),sTauBar_p4.Eta());
  //TODO: Correct dR implementation. 
  angle = sTauBar_p4.Angle(sTau_p4.Vect());
  fHistAngle->Fill(angle);
  fHistDeltaPhi->Fill(deltaPhi);
 
  for(auto& iHit : *hits)
  {
    if((int)abs(iHit.particleType())== 1000015) {
      DetId theDetUnitId = DetId(iHit.detUnitId());
      DetId simdetid= DetId(iHit.detUnitId());
      
      if(simdetid.det()==DetId::Muon &&  simdetid.subdetId()== MuonSubdetId::RPC){
    	RPCDetId rollId(theDetUnitId);
    	RPCGeomServ rpcsrv(rollId);
    	//const RPCRoll * rollasociated = rpcGeo->roll(rollId);
    	//const BoundPlane & RPCSurface = rollasociated->surface();
    	//GlobalPoint SimHitInGlobal = RPCSurface.toGlobal((*iHit).localPosition());
    	if(rollId.station()==1 && rollId.region()==1) fHisttofHits->Fill(iHit.timeOfFlight());
    	if(rollId.station()==1 && rollId.region()==0 && rollId.layer()==1) fHisttofBarrel1_in->Fill(iHit.timeOfFlight());
    	if(rollId.station()==1 && rollId.region()==0 && rollId.layer()==2) fHisttofBarrel1_out->Fill(iHit.timeOfFlight());
    	if(rollId.station()==2 && rollId.region()==0 && rollId.layer()==1) fHisttofBarrel2_in->Fill(iHit.timeOfFlight());
    	if(rollId.station()==2 && rollId.region()==0 && rollId.layer()==2) fHisttofBarrel2_out->Fill(iHit.timeOfFlight());
    	if(rollId.station()==3 && rollId.region()==0) fHisttofBarrel3->Fill(iHit.timeOfFlight());
    	if(rollId.station()==4 && rollId.region()==0) fHisttofBarrel4->Fill(iHit.timeOfFlight());
    	if(rollId.station()==1 && rollId.region()==1) fHisttofEndCapF1->Fill(iHit.timeOfFlight());
    	if(rollId.station()==2 && rollId.region()==1) fHisttofEndCapF2->Fill(iHit.timeOfFlight());
    	if(rollId.station()==3 && rollId.region()==1) fHisttofEndCapF3->Fill(iHit.timeOfFlight());
    	if(rollId.station()==4 && rollId.region()==1) fHisttofEndCapF4->Fill(iHit.timeOfFlight());
    	if(rollId.station()==1 && rollId.region()==-1) fHisttofEndCapB1->Fill(iHit.timeOfFlight());
    	if(rollId.station()==2 && rollId.region()==-1) fHisttofEndCapB2->Fill(iHit.timeOfFlight());
    	if(rollId.station()==3 && rollId.region()==-1) fHisttofEndCapB3->Fill(iHit.timeOfFlight());
    	if(rollId.station()==4 && rollId.region()==-1) fHisttofEndCapB4->Fill(iHit.timeOfFlight());
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
  fHistSTauP = fs->make<TH1D>("HistSTauP", "sTau P", 300, 0., 3000. );
  fHistSTauEta = fs->make<TH1D>("HistSTauEta", "sTau Eta", 150, -4., 4.);
  fHistSTauEtaBeta = fs->make<TH2D>("HistSTauEtaBeta", "sTau Eta vs Beta", 150, 0., 1.,150,0.,4.);
  fHistEta = fs->make<TH2D>("HistEta","sTau Eta vs sTauBar Eta",150, -4., 4.,150, -4., 4.);
  fHistSTauBarEta = fs->make<TH1D>("HistSTauBarEta","sTau bar Eta",150,-4.,4.);
  fHistSTauBeta = fs->make<TH1D>("HistSTauBeta", "sTau #Beta", 150, 0, 1.0);
  fHistSTauPhi = fs->make<TH1D>("HistSTauPhi","sTau #phi", 150, -3.15,3.15);
  fHistAngle = fs->make<TH1D>("HistAngle", "sTau-sTau bar angle",150, 0.,3.15);
  fHistDeltaPhi = fs->make<TH1D>("HistDeltaPhi","abs(Phi(sTau)-Phi(sTauBar))",150,0.,3.15);
  fHisttof = fs->make<TH1D>("Histtof","tof",150, 0.,40000.);
  fHisttofHits = fs->make<TH1D>("HisttofHits","tofHits",150,0.0,400.);
  fHistAcc = fs->make<TH1D>("HistoAcc","Acceptance ring 1",10,0,10);
  fHisttofBarrel1_in = fs->make<TH1D>("HistotofBarrelStation1_in","Dtof Barrel Station 1 in",150,0.,400.);
  fHisttofBarrel1_out = fs->make<TH1D>("HistotofBarrelStation1_out","Dtof Barrel Station 1 out",150,0.,400.);
  fHisttofBarrel2_in = fs->make<TH1D>("HistotofBarrelStation2_in","Dtof Barrel Station 2 in",150,0.,400.);
  fHisttofBarrel2_out = fs->make<TH1D>("HistotofBarrelStation2_out","Dtof Barrel Station 2 out",150,0.,400.);
  fHisttofBarrel3 = fs->make<TH1D>("HistotofBarrelStation3","Dtof Barrel Station 3",150,0.,400.);
  fHisttofBarrel4 = fs->make<TH1D>("HistotofBarrelStation4","Dtof Barrel Station 4",150,0.,400.);
  fHisttofEndCapF1 = fs->make<TH1D>("HistotofEndCapStationF1","Dtof EndCap Station 1(Forward)",150,0.,400.);
  fHisttofEndCapF2 = fs->make<TH1D>("HistotofEndCapStationF2","Dtof EndCap Station 2(Forward)",150,0.,400.);
  fHisttofEndCapF3 = fs->make<TH1D>("HistotofEndCapStationF3","Dtof EndCap Station 3(Forward)",150,0.,400.);
  fHisttofEndCapF4 = fs->make<TH1D>("HistotofEndCapStationF4","Dtof EndCap Station 4(Forward)",150,0.,400.);
  fHisttofEndCapB1 = fs->make<TH1D>("HistotofEndCapStationB1","Dtof EndCap Station 1(Backward)",150,0.,400.);
  fHisttofEndCapB2 = fs->make<TH1D>("HistotofEndCapStationB2","Dtof EndCap Station 2(Backward)",150,0.,400.);
  fHisttofEndCapB3 = fs->make<TH1D>("HistotofEndCapStationB3","Dtof EndCap Station 3(Backward)",150,0.,400.);
  fHisttofEndCapB4 = fs->make<TH1D>("HistotofEndCapStationB4","Dtof EndCap Station 4(Backward)",150,0.,400.);
  return;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HSCPSIMHits::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(HSCPSIMHits);
