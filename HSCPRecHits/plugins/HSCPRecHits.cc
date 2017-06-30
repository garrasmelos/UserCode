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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include <Geometry/RPCGeometry/interface/RPCRoll.h>
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/CommonDetUnit/interface/GeomDet.h>
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/RPCDigi/interface/RPCDigi.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
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
      
      std::vector<double>  doFit(std::vector<TVector3> POS,std::vector<double> TIME);
      UInt_t getTriggerBits(const edm::Event& iEvent, std::vector<std::string> TestFilterNames_);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
    TTree* tree_;
    const static unsigned short rpcHits_N=100;
    unsigned short b_rpcHits_n;
    unsigned int b_trig;
    double b_rpcBeta, b_rpcBetaErr, b_fitSlope, b_genBeta, b_t0, b_t0_bx, b_genEta; 
    double b_rpcHitTime[rpcHits_N], b_rpcHitPos[rpcHits_N],b_rpcHitTimeErr[rpcHits_N], b_rpcHitPosErr[rpcHits_N];
      
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genParToken_;
    edm::EDGetTokenT<RPCRecHitCollection> recHitToken_;
    edm::EDGetTokenT<edm::PSimHitContainer> simHitToken_;  
    edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
    std::vector<std::string> triggNames_;
    int pId_;
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
:  genParToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticlesLabel"))),
  recHitToken_(consumes<RPCRecHitCollection>(iConfig.getParameter<edm::InputTag>("recHitLabel"))),
  simHitToken_(consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("simHitLabel"))),
  triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerLabel"))),
  triggNames_(iConfig.getParameter<std::vector<std::string>>("triggNames")),
  pId_(iConfig.getParameter<int>("particleId"))
{
   //now do what ever initialization is needed
}
HSCPRecHits::~HSCPRecHits()
{
}

//
// member functions
//
std::vector<double>
HSCPRecHits::doFit(std::vector<TVector3> POS,std::vector<double> TIME)
{ 
  double s = 0;
  double sx = 0;
  double sy = 0;
  double sxy = 0;
  double sxx = 0;
  double syy = 0;
  double ssxy = 0;
  double ssxx = 0;
  double ssyy = 0;
  double a=0;
  double aStdErr=0;
  double b=0;
  double bStdErr=0;
  
  const double n = POS.size();
  for(int i = 0 ; i < n;i++)
  {
    sx+= POS[i].Mag();
    sy+= TIME[i];
    sxy+= POS[i].Mag()*TIME[i];
    sxx+= POS[i].Mag()*POS[i].Mag();
    syy+= TIME[i]*TIME[i];
  }
  ssxy = sxy-sx*sy/n;
  ssxx = sxx-sx*sx/n;
  ssyy = syy-sy*sy/n;
  b = ssxy/ssxx;
  a = (1/n)*(sy-b*sx);
  s = TMath::Sqrt((ssyy - b*ssxy)/(n-2));
  aStdErr = s * TMath::Sqrt((1/n)+(sx*sx/n*n*ssxx));
  bStdErr = s/TMath::Sqrt(ssxx);
  
  vector<double> parameters;
  parameters.push_back(a);
  parameters.push_back(aStdErr);
  parameters.push_back(b);
  parameters.push_back(bStdErr);
  
  return parameters;
}

UInt_t HSCPRecHits::getTriggerBits(const edm::Event& iEvent, std::vector<std::string> TestFilterNames_)
{
  UInt_t trigger=0;
  edm::Handle<edm::TriggerResults> hltTriggerResultHandle;
  iEvent.getByToken(triggerToken_,hltTriggerResultHandle);
   
  if(hltTriggerResultHandle.isValid())
  {
    const edm::TriggerNames & triggerNames = iEvent.triggerNames(*hltTriggerResultHandle);
    //if ( triggerNames.size()>0 ) for (unsigned int l=0; l < triggerNames.size() ;l++) std::cout << triggerNames.triggerName(l) << std::endl;
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
  }// else std::cout << "HSCPTrigger::getTriggerBits: **** No triggers found ****" << std::endl;
   return trigger;
}	
// ------------ method called for each event  ------------
void
HSCPRecHits::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  Handle<PSimHitContainer> simhits;
  iEvent.getByToken(simHitToken_,simhits); 

  Handle<RPCRecHitCollection> rechits;
  iEvent.getByToken(recHitToken_,rechits);
  
  edm::ESHandle<RPCGeometry> rpcGeo;
  iSetup.get<MuonGeometryRecord>().get(rpcGeo);
   
  edm::Handle<vector<reco::GenParticle>> genParHandle;
  iEvent.getByToken(genParToken_,genParHandle);

  vector<TVector3> vPos; //array of positions.
  vector<TVector3> vPosSim;
  vector<double> vTime;
  vector<double> vBx;
  
  const double c = 29.979;
  double betaGen=0, etaGen=0;
  b_rpcBeta=b_rpcBetaErr=b_t0=b_t0_bx=b_fitSlope=-10000.;
  
  TVector3 hscpDir(0,0,0);
 
  b_trig = getTriggerBits(iEvent,triggNames_);
  
  if(genParHandle.isValid())
  {
    for (auto& pout : *genParHandle)
    { 
      if(pout.pdgId() == pId_ && pout.status()==1)
      {
        double sTauP = pout.p4().P();
        double sTauMass= pout.p4().M();
        betaGen= sqrt(sTauP*sTauP/(sTauP*sTauP+sTauMass*sTauMass));
        hscpDir.SetXYZ(pout.p4().Px(),pout.p4().Py(),pout.p4().Pz());       
        etaGen= hscpDir.Eta();
      }
    }
  }
  if(rechits.isValid() && hscpDir.Mag()!=0.)
  { 
    int j=0;
    for(RPCRecHitCollection::const_iterator rechit_it = rechits->begin(); rechit_it != rechits->end() ; rechit_it++)
    {
      double time = rechit_it->time();
      double timeErr = rechit_it->timeError();
      int bx = rechit_it->BunchX();
      RPCDetId idRoll(rechit_it->rpcId());
      LocalPoint lPos = rechit_it->localPosition();
      const RPCRoll* roll = rpcGeo->roll(idRoll);
      double stripLength = roll->specificTopology().stripLength();
      const BoundPlane& rollSurface = roll->surface();
      GlobalPoint gPos = rollSurface.toGlobal(lPos);
      TVector3 pos(gPos.x(),gPos.y(),gPos.z());
      
      //Create a point to compute the error in the position.
      LocalPoint lPosEdge(lPos.x(),stripLength/2,0.);
      GlobalPoint gPosEdge = rollSurface.toGlobal(lPosEdge);
      TVector3 posEdge(gPosEdge.x(),gPosEdge.y(),gPosEdge.z());


      double dr = hscpDir.DeltaR(pos);
      //Filling arrays with time and position of the rpdHits found inside a cone around de 
      //direction of the Gen particle. 
      if(dr < 0.5)
      {
        b_rpcHitTime[j] = time;
        b_rpcHitPos[j] = pos.Mag();
        b_rpcHitTimeErr[j] =  timeErr;//stripLength/c;
        b_rpcHitPosErr[j] = TMath::Abs(pos.Mag()-posEdge.Mag());
        j++;
        vPos.push_back(pos);
        vBx.push_back(bx*25.);
        vTime.push_back(time);
      }
    }
    if(j>20)cout << j << endl;
    b_rpcHits_n = j;
  }
  int nSimHits=0;
  if(simhits.isValid())
  {
    for(auto& simhit_it : *simhits)
    {
      if(simhit_it.particleType() != pId_) continue;
      RPCDetId idRollSim(simhit_it.detUnitId());
      LocalPoint lPosSim = simhit_it.localPosition();
      const RPCRoll* rollSim = rpcGeo->roll(idRollSim);
      const BoundPlane& rollSurfaceSim = rollSim->surface();
      GlobalPoint gPosSim = rollSurfaceSim.toGlobal(lPosSim);
      TVector3 posSim(gPosSim.x(), gPosSim.y(), gPosSim.z()); 
      
      double drSim = hscpDir.DeltaR(posSim);
      
      if(drSim < 0.5) vPosSim.push_back(posSim);
      nSimHits++;
    }
  }
  vector<double> params;
  vector<double> params_bx;
  if (vPos.size() > 2)  
  {
    params = doFit(vPos,vTime);

    b_genBeta = betaGen;
    b_genEta = etaGen;
    b_t0 = params[0];
    b_rpcBeta = 1/((params[2]*c)+1);
    b_rpcBetaErr = (c*params[3])/((params[2]*c+1)*(params[2]*c+1));
    b_fitSlope = params[2];
    
    //Fit using bx
    params_bx = doFit(vPos,vBx);
    b_t0_bx = params_bx[0];
  }
  if(b_rpcHits_n < 20) tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
HSCPRecHits::beginJob()
{
  Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree","tree");
  tree_->Branch("rpcHits_n", &b_rpcHits_n, "rpcHits_n/s");
  tree_->Branch("rpcBeta", &b_rpcBeta, "rpcBeta/d");
  tree_->Branch("rpcBetaErr", &b_rpcBetaErr, "rpcBetaErr/d");
  tree_->Branch("rpcHitTime", b_rpcHitTime, "rpcHitTime[rpcHits_n]/D");
  tree_->Branch("rpcHitTimeErr", b_rpcHitTimeErr, "rpcHitTimeErr[rpcHits_n]/D");
  tree_->Branch("rpcHitPos", b_rpcHitPos, "rpcHitPos[rpcHits_n]/D");
  tree_->Branch("rpcHitPosErr", b_rpcHitPosErr, "rpcHitPosErr[rpcHits_n]/D");
  tree_->Branch("t0", &b_t0, "t0/d");
  tree_->Branch("t0_bx", &b_t0_bx, "t0_bx/d");
  tree_->Branch("genBeta", &b_genBeta, "genBeta/d");
  tree_->Branch("genEta", &b_genEta, "genEta/d");
  tree_->Branch("fitSlope", &b_fitSlope, "fitSlope/d");
  tree_->Branch("trig", &b_trig, "trig/I");
  return;
}
DEFINE_FWK_MODULE(HSCPRecHits);
