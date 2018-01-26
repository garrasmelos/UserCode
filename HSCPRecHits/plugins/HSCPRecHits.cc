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
      void clearVectors();
      void findRpcHits(const edm::Event&);
      void findRpcHits2(const edm::Event&);
      void findSimHits(const edm::Event&);
      TVector3 getGPos(const RPCRecHit&);
      bool isFirstHit(int s, int l, int r);
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
      TTree* tree2_;
      TTree* tree3_;
      const double c = 29.979;
      const double nMinHits = 2;
      const static unsigned short rpcHits_N=100;
      unsigned short b_rpcHits_n,b_tracks_N;
      unsigned int b_trig;
    
      int b_charge;
      double b_rpcBeta, b_rpcBetaErr, b_fitSlope, b_genBeta, b_t0, b_t0_bx, b_rpcBeta_bx, b_genEta, b_genP, b_mass; 
      double b_rpcHitTime[rpcHits_N], b_rpcHitPos[rpcHits_N],b_rpcHitTimeErr[rpcHits_N], b_rpcHitPosErr[rpcHits_N];
      double b_dR, b_gp_x;
      TVector3 hscpDir;
    
      std::vector<TVector3> vPos; //array of positions.
      std::vector<TVector3> vPosSim;
      std::vector<double> vTime;
      std::vector<double> vBx;
        
      edm::ESHandle<RPCGeometry> rpcGeo;
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
: genParToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticlesLabel"))),
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

// ------------ method called for each event  ------------
void
HSCPRecHits::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<vector<reco::GenParticle>> genParHandle;
  iEvent.getByToken(genParToken_,genParHandle);

  iSetup.get<MuonGeometryRecord>().get(rpcGeo);

  b_trig = getTriggerBits(iEvent,triggNames_);
  // Getting the generated particle intial direction. 
  if(genParHandle.isValid())
  {
    for (auto& pout : *genParHandle)
    { 
      b_rpcBeta=b_rpcBetaErr=b_t0=b_t0_bx=b_rpcBeta_bx=b_fitSlope=-10000.;
      //cout << pout.pdgId() << endl;
      //cout << "Status" << pout.status() << endl;
      
      if(TMath::Abs(pout.pdgId()) == pId_ && pout.status()==1)
      {
        double sTauP = pout.p4().P();
        double sTauMass= pout.p4().M();
        b_genBeta= sqrt(sTauP*sTauP/(sTauP*sTauP+sTauMass*sTauMass));
        hscpDir.SetXYZ(pout.p4().Px(),pout.p4().Py(),pout.p4().Pz());       
        b_genEta= hscpDir.Eta();
        b_charge = pout.charge();
        b_genP = sTauP;
        //cout << etaGen << endl;
        //cout << pout.charge() << endl;
        //Finding the rpcHits and simHits inside the cone define 
        //by a dR of 0.5 around the generated particle dir.
        //cout << "hscp x " << hscpDir.x() << endl;

        if(hscpDir.Mag()!=0.)
        { 
          //int j=0;
          clearVectors();
          //cout << "Searching for HSCP hits..." << endl;
          findRpcHits(iEvent);
          //findRpcHits2(iEvent);
          findSimHits(iEvent);
          int j = vPos.size();
          b_rpcHits_n = j;
          //cout << j <<  " hits were found." << endl;
        }
        vector<double> params,params_bx;
        if (vPos.size() >= nMinHits)  
        {
          params = doFit(vPos,vTime);

          b_t0 = params[0];
          b_rpcBeta = 1/((params[2]*c)+1);
          double gamma = 1/TMath::Sqrt(1-b_rpcBeta*b_rpcBeta);
          b_mass = sTauP/(b_rpcBeta*gamma);
          if(vPos.size()==2) b_rpcBetaErr = 0;
          else b_rpcBetaErr = (c*params[3])/((params[2]*c+1)*(params[2]*c+1));
          
          //if (etaGen > 1.8) cout << (c*params[3])/((params[2]*c+1)*(params[2]*c+1)) << endl;
          b_fitSlope = params[2];
          
          //Fit using bx
          params_bx = doFit(vPos,vBx);
          b_t0_bx = params_bx[0];
          b_rpcBeta_bx = 1/((params_bx[2]*c)+1);
        }
        if(b_rpcHits_n < rpcHits_N) tree_->Fill();
      }
    }
  }
  findRpcHits2(iEvent); 
}

//
// member functions
//
void HSCPRecHits::findRpcHits(const edm::Event& event)
{
  edm::Handle<RPCRecHitCollection> rechits;
  event.getByToken(recHitToken_,rechits);
  int j = 0;
  if(rechits.isValid())
  {  
    for(RPCRecHitCollection::const_iterator rechit_it = rechits->begin(); rechit_it != rechits->end() ; rechit_it++)
    {
      double time = rechit_it->time();
      double timeErr = rechit_it->timeError();
      int bx = rechit_it->BunchX();
      LocalPoint lPos = rechit_it->localPosition(); 
      TVector3 pos = getGPos(*rechit_it);
      
      //Create a point to compute the error in the position.
      RPCDetId idRoll(rechit_it->rpcId());
      const RPCRoll* roll = rpcGeo->roll(idRoll);
      double stripLength = roll->specificTopology().stripLength();
      const BoundPlane& rollSurface = roll->surface();
      LocalPoint lPosEdge(lPos.x(),stripLength/2,0.);
      GlobalPoint gPosEdge = rollSurface.toGlobal(lPosEdge);
      TVector3 posEdge(gPosEdge.x(),gPosEdge.y(),gPosEdge.z());
      double dr = hscpDir.DeltaR(pos);
      b_dR = dr;
      b_gp_x = pos.x();
      tree2_->Fill();
      //cout << dr << endl;    //Filling arrays with time and position of the rpdHits found inside a cone around de 
      //cout << pos.Mag() << endl;
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
  }
  return;
}
TVector3 HSCPRecHits::getGPos(const RPCRecHit& hit)
{
    RPCDetId idRoll(hit.rpcId());
    LocalPoint lPos = hit.localPosition(); const RPCRoll* roll = rpcGeo->roll(idRoll);
    const BoundPlane& rollSurface = roll->surface();
    GlobalPoint gPos = rollSurface.toGlobal(lPos);
    TVector3 position(gPos.x(),gPos.y(),gPos.z());
    return position;
}
bool HSCPRecHits::isFirstHit(int s, int l, int r)
{
  bool firstHit=false;  
  
  int tl =4; //total layers.
  int cl =s; //Current layer.
  
  if(TMath::Abs(r)==0) 
  {
    tl = 6;
    if (s == 1) cl = s + l;
    else if (s == 2) cl = s + l + 1;
    else if (s > 2 ) cl = s + 2;
  }
  if(cl <= (tl-nMinHits) ) firstHit=true;
  return firstHit;  
}
void HSCPRecHits::findRpcHits2(const edm::Event& event)
{
  edm::Handle<RPCRecHitCollection> rechits;
  event.getByToken(recHitToken_,rechits);
  if(rechits->size()<nMinHits) return;
  vector<vector<TVector3>> p;
  vector<vector<double>> t;
  vector<vector<int>> bx;
  vector< unsigned int> trlist;

  unsigned int nH = 0;
  for(auto& rpcHit : *rechits ) 
  {
    RPCDetId id(rpcHit.rpcId());
    int station = id.station();
    int layer = id.layer();
    int region = id.region();
    if(isFirstHit(station, layer, region))
    {
      TVector3 position1 = getGPos(rpcHit);
      unsigned int j =0;
      for(auto& rpcHit2 : *rechits)
      {
        TVector3 position2 = getGPos(rpcHit2);
        double deltaR = TMath::Abs(position2.DeltaR(position1));
        if(deltaR < 0.5)
        { 
          trlist.push_back(j);
          j =0;
          break;
        }
        //else if (j== nvectors-1) trlist.push_back(1000);
        j++;
      }
    }
    nH++;
  }
      
  //cout << "New event" << endl;  
  int nTracks=0;
  for(unsigned int s=0; s<6; s++)
  {
    vector<double> vT;
    unsigned int t=0;
    for(auto& rpcHit3 : *rechits)
    {
      if(0<trlist.size())
      {
        if(s == trlist[t])
        {
          vT.push_back(rpcHit3.time());
        } 
      }
      t++;
    } 
    if(vT.size()>= nMinHits)
    {
      //cout << vT.size() << endl; 
      nTracks++;
    }
  }
  b_tracks_N = nTracks;
  tree3_->Fill();
  return; 

}
void HSCPRecHits::findSimHits(const edm::Event& event)
{
  Handle<PSimHitContainer> simhits;
  event.getByToken(simHitToken_,simhits); 

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
      //cout << drSim << endl;
      
      if(drSim < 0.5) vPosSim.push_back(posSim);
      nSimHits++;
    }
  }
}
void HSCPRecHits::clearVectors()
{
  vPos.clear();
  vPosSim.clear();
  vTime.clear();
  vBx.clear();
  return;
}
std::vector<double> HSCPRecHits::doFit(std::vector<TVector3> POS,std::vector<double> TIME)
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
// ------------ method called once each job just before starting event loop  ------------
void 
HSCPRecHits::beginJob()
{
  Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree","tree");
  tree2_ = fs->make<TTree>("tree2","tree2");
  tree3_ = fs->make<TTree>("tree3","tree3");
  tree_->Branch("rpcHits_n", &b_rpcHits_n, "rpcHits_n/s");
  tree_->Branch("rpcBeta", &b_rpcBeta, "rpcBeta/d");
  tree_->Branch("rpcBetaErr", &b_rpcBetaErr, "rpcBetaErr/d");
  tree_->Branch("rpcHitTime", b_rpcHitTime, "rpcHitTime[rpcHits_n]/D");
  tree_->Branch("rpcHitTimeErr", b_rpcHitTimeErr, "rpcHitTimeErr[rpcHits_n]/D");
  tree_->Branch("rpcHitPos", b_rpcHitPos, "rpcHitPos[rpcHits_n]/D");
  tree_->Branch("rpcHitPosErr", b_rpcHitPosErr, "rpcHitPosErr[rpcHits_n]/D");
  tree_->Branch("t0", &b_t0, "t0/d");
  tree_->Branch("t0_bx", &b_t0_bx, "t0_bx/d");
  tree_->Branch("rpcBeta_bx", &b_rpcBeta_bx, "rpcBeta_bx/d");
  tree_->Branch("genBeta", &b_genBeta, "genBeta/d");
  tree_->Branch("genEta", &b_genEta, "genEta/d");
  tree_->Branch("genP", &b_genP, "genP/d");
  tree_->Branch("mass", &b_mass, "mass/d");
  tree_->Branch("fitSlope", &b_fitSlope, "fitSlope/d");
  tree_->Branch("trig", &b_trig, "trig/I");
  tree_->Branch("charge", &b_charge, "charge/I");
  tree2_->Branch("dR",&b_dR,"dR/d");
  tree2_->Branch("x",&b_gp_x,"x/d");
  tree3_->Branch("tracks_N", &b_tracks_N, "tracks_N/s");
  return;
}
DEFINE_FWK_MODULE(HSCPRecHits);
