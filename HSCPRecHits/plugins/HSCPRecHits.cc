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
    TH1D* fHnHits;
    TH1D* fHnSimHits;
    TH1D* fHnHitsDiff;
    TH1D* fHt0;
    TH1D* fHt0err;
    TH1D* fHbetaGen;
    TH1D* fHb;
    TH1D* fHberr;
    TH1D* fHbetaRPC;
    TH1D* fHetaGen;
    TH1D* fHbeta0hits;
    TH1D* fHeta0hits;
    TH1D* fHbeta1hits;
    TH1D* fHeta1hits;
    TH1D* fHbeta2hits;
    TH1D* fHeta2hits;
    TH1D* fHbeta_pas;
    TH1D* fHeta_pas;
    TH1D* fHbeta_pas_SlopeCut;
    TH1D* fHbeta_pas_BetaErrorCut;
    TH1D* fHbeta_pas_BetaRelError;
    TH1D* fHbeta_tot;
    TH1D* fHeta_tot;
    TH2D* fHbetaEta_pas;
    TH1D* fHres;
    TH1D* fHdR;
    
      
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genParToken_;
    edm::EDGetTokenT<RPCRecHitCollection> recHitToken_;
    edm::EDGetTokenT<edm::PSimHitContainer> simHitToken_;  
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
:  fHnHits(0), fHnSimHits(0), fHnHitsDiff(0), fHt0(0), fHt0err(0), fHbetaGen(0), fHb(0), fHberr(0), fHbetaRPC(0), fHetaGen(0), fHbeta0hits(0), fHeta0hits(0), fHbeta1hits(0), fHeta1hits(0),fHbeta2hits(0), fHeta2hits(0),fHbeta_pas(0), fHeta_pas(0), fHbeta_pas_SlopeCut(0), fHbeta_pas_BetaErrorCut(0), fHbeta_pas_BetaRelError(0),  fHbeta_tot(0), fHeta_tot(0),  fHbetaEta_pas(0), fHres(0), fHdR(0),
  genParToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticlesLabel"))),
  recHitToken_(consumes<RPCRecHitCollection>(iConfig.getParameter<edm::InputTag>("recHitLabel"))),
  simHitToken_(consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("simHitLabel")))
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
  vector<double> vTime;
  
  const double c = 29.979;
  double betaGen=0, etaGen=0;
  TVector3 hscpDir;
   
  if(genParHandle.isValid())
  {
    for (auto& pout : *genParHandle)
    { 
      if( pout.pdgId() == 1000015 && pout.status()==1)
      {
        double sTauP = pout.p4().P();
        double sTauMass= pout.p4().M();
        betaGen= sqrt(sTauP*sTauP/(sTauP*sTauP+sTauMass*sTauMass));
        hscpDir.SetXYZ(pout.p4().Px(),pout.p4().Py(),pout.p4().Pz());       
        etaGen= hscpDir.Eta();
      }
    }
  }
  if(rechits.isValid())
  {
    for(RPCRecHitCollection::const_iterator rechit_it = rechits->begin(); rechit_it != rechits->end() ; rechit_it++)
    {
      double time = rechit_it->time();
      RPCDetId idRoll(rechit_it->rpcId());
      LocalPoint lPos = rechit_it->localPosition();
      const RPCRoll* roll = rpcGeo->roll(idRoll);
      const BoundPlane& rollSurface = roll->surface();
      GlobalPoint gPos = rollSurface.toGlobal(lPos);
      TVector3 pos(gPos.x(),gPos.y(),gPos.z());
      
      double dr = hscpDir.DeltaR(pos);
      fHdR->Fill(dr);

      if(dr < 0.5)
      {
        vPos.push_back(pos);
        vTime.push_back(time);
      }
    }
  }
  int nSimHits=0;
  if(simhits.isValid())
  {
    for(auto& simhit : *simhits)
    {
      const int pid = simhit.particleType();
      if ( pid == 1000015) nSimHits++;
    }
    fHnSimHits->Fill(nSimHits);
    if(vPos.size() < 3) fHnHitsDiff->Fill(nSimHits-vPos.size());
  }
  fHbetaGen->Fill(betaGen);
  fHetaGen->Fill(etaGen);

  vector<double> params;
  if (vPos.size() == 0) 
  {
    fHbeta0hits->Fill(betaGen);
    fHeta0hits->Fill(etaGen);
  }
  if (vPos.size() == 1)
  {
    fHbeta1hits->Fill(betaGen);
    fHeta1hits->Fill(etaGen);
  }
  if (vPos.size() == 2)
  {
    fHbeta2hits->Fill(betaGen);
    fHeta2hits->Fill(etaGen);
  }
  if (vPos.size() > 2)  
  {
    fHbeta_tot->Fill(betaGen);
    fHeta_tot->Fill(etaGen);
    params = doFit(vPos,vTime);
    fHt0->Fill(params[0]);
    fHt0err->Fill(params[1]);
    fHb->Fill(params[2]);
    fHberr->Fill(params[3]);

    double betaRPC = 1/((params[2]*c)+1);
    double betaRPCerr = (c*params[3])/((params[2]*c+1)*(params[2]*c+1));
    double errFitCut = 0.3;
    fHbetaRPC->Fill(betaRPC);
    
    if (TMath::Abs(betaRPCerr/betaRPC) < errFitCut && params[2]>0. ) //
    {
      double res = (betaGen-betaRPC)/betaGen;
      fHres->Fill(res);
      fHbeta_pas->Fill(betaGen);
      fHeta_pas->Fill(etaGen);
      fHbetaEta_pas->Fill(betaGen,etaGen);
    }
    //if ( params[2]>0. || (params[2]+params[3])>0. || (params[2]-params[3])>0. )
    if ( params[2]>0.)
    {
      fHbeta_pas_SlopeCut->Fill(betaGen);
    }
    if (TMath::Abs(betaRPCerr/betaRPC) < errFitCut)
    {
      fHbeta_pas_BetaErrorCut->Fill(betaGen);
    }
    fHbeta_pas_BetaRelError->Fill(betaRPCerr/betaRPC);
  }
  fHnHits->Fill(vPos.size());
}


// ------------ method called once each job just before starting event loop  ------------
void 
HSCPRecHits::beginJob()
{
  Service<TFileService> fs;
  fHnHits = fs->make<TH1D>("fHnHits", "n hits in sTau direction", 20, 0., 20. );
  fHnSimHits = fs->make<TH1D>("fHnSimHits", "n stau simhits ", 20, 0., 20. );
  fHnHitsDiff = fs->make<TH1D>("fHnHitsDiff", "Difference between the # of stau simhits and rechits ", 20, 0., 20. );
  fHt0 = fs->make<TH1D>("fHt0","initial time",100,-25.,25.);
  fHt0err = fs->make<TH1D>("fHt0err","initial time error",100,-25.,25.);
  fHbetaGen = fs->make<TH1D>("fHbetaGen","beta Gen",150,-1.,2.);
  fHetaGen = fs->make<TH1D>("fHetaGen","eta Gen",50,-3.15,3.15);
  fHb = fs->make<TH1D>("fHb","b",80,-2.,2.);
  fHberr = fs->make<TH1D>("fHberr","b error",60,-3.,3.);
  fHbetaRPC = fs->make<TH1D>("fHbetaRPC","betaRPC",80,-2.,2.);
  fHbeta0hits = fs->make<TH1D>("fHbeta0hits","beta for stau without rechits",50,0.,1.);
  fHeta0hits = fs->make<TH1D>("fHeta0hits","eta for stau without rechits",50,-3.15,3.15);
  fHbeta1hits = fs->make<TH1D>("fHbeta1hits","beta for stau with 1 rechits",50,0.,1.);
  fHeta1hits = fs->make<TH1D>("fHeta1hits","eta for stau with 1 rechits",50,-3.15,3.15);
  fHbeta2hits = fs->make<TH1D>("fHbeta2hits","beta for stau with 2 rechits",50,0.,1.);
  fHeta2hits = fs->make<TH1D>("fHeta2hits","eta for stau with 2 rechits",50,-3.15,3.15);
  fHbeta_pas = fs->make<TH1D>("fHbeta_pas","beta generated pas",50,0.,1.);
  fHbeta_pas_SlopeCut = fs->make<TH1D>("fHbeta_pas_SlopeCut","beta generated pas Slope slope cut",50,0.,1.);
  fHbeta_pas_BetaErrorCut = fs->make<TH1D>("fHbeta_pas_BetaErrorCut","beta generated pas BetaErrorCut",50,0.,1.);
  fHbeta_pas_BetaRelError = fs->make<TH1D>("fHbeta_pas_BetaRelError","beta generated relative error",100,-2.,2.);
  fHbeta_tot = fs->make<TH1D>("fHbeta_tot","betaGen afteer 3 hits cut",50,0.,1.);
  fHeta_tot = fs->make<TH1D>("fHeta_tot","etaGen after 3 hits cut",50,-3.15,3.15);
  fHeta_pas = fs->make<TH1D>("fHeta_pas","etaGen after all cuts",50,-3.15,3.15);
  fHbetaEta_pas = fs->make<TH2D>("fHbetaEta_pas","fHbetaEta_pas", 50,0.,1., 50, -3.15,3.15);
  fHres = fs->make<TH1D>("fHres","Beta resolution",60,-3.,3.);
  fHdR = fs->make<TH1D>("fHdR","fHdR",50,0.,5.);
  return;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HSCPRecHits::endJob() 
{
}
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
