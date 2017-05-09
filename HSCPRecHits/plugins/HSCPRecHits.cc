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
      
      std::vector<double>  doFit(std::vector<TVector3> POS,std::vector<int> BX);

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
		TH1D* fHt0;
		TH1D* fHbeta;
		TH1D* fHeta_tot;
		TH1D* fHbeta0;
		TH1D* fHeta0;
		TH1D* fHbeta_pas;
		TH1D* fHbeta_tot;
		TH1D* fHres;
		TH1D* fHdR;
		
      
		edm::EDGetTokenT<std::vector<reco::GenParticle>> genParToken_;
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
:  fHnHits(0), fHt0(0), fHbeta(0), fHeta_tot(0), fHbeta0(0), fHeta0(0), fHbeta_pas(0), fHbeta_tot(0), fHres(0), fHdR(0),
	genParToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticlesLabel"))),
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
std::vector<double>
HSCPRecHits::doFit(std::vector<TVector3> POS,std::vector<int> BX)
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
		sy+= BX[i];
		sxy+= POS[i].Mag()*BX[i];
		sxx+= POS[i].Mag()*POS[i].Mag();
		syy+= BX[i]*BX[i];
		
	}
	ssxy = sxy-sx*sy/n;
	ssxx = sxx-sx*sx/n;
	ssyy = syy-sy*sy/n;
	b = ssxy/ssxx;
	a = (1/n)*(sy-b*sx);
	s = TMath::Sqrt((ssyy - b*ssxy)/(n-2));
	aStdErr = s * TMath::Sqrt((1/n)+(sx*sx/n*n*ssxx));
	aStdErr = s/TMath::Sqrt(ssxx);
	
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
   
  Handle<RPCRecHitCollection> rechits;
  iEvent.getByToken(recHitToken_,rechits);
  
  edm::ESHandle<RPCGeometry> rpcGeo;
  iSetup.get<MuonGeometryRecord>().get(rpcGeo);
   
  edm::Handle<vector<reco::GenParticle>> genParHandle;
	iEvent.getByToken(genParToken_,genParHandle);

  vector<TVector3> vPos; //array of positions.
	vector<int> vBx;
  
  const double c = 29.979;
	double betaGen=0;
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
      }
    }
  }

  for(RPCRecHitCollection::const_iterator rechit_it = rechits->begin(); rechit_it != rechits->end() ; rechit_it++)
  {
   	int bx = rechit_it->BunchX();
   	RPCDetId idRoll(rechit_it->rpcId());
   	LocalPoint lPos = rechit_it->localPosition();
   	const RPCRoll* roll = rpcGeo->roll(idRoll);
   	const BoundPlane& rollSurface = roll->surface();
   	GlobalPoint gPos = rollSurface.toGlobal(lPos);
   	TVector3 pos(gPos.x(),gPos.y(),gPos.z());
   	
    double dr = hscpDir.DeltaR(pos);
   	fHdR->Fill(dr);
   	
    if(dr < 0.15)
		{
			vPos.push_back(pos);
			vBx.push_back(bx);
		}
  }
	vector<double> params;
	if (vPos.size() == 0) 
	{
		if(betaGen < 0.15) fHeta0->Fill(hscpDir.Eta());
		fHbeta0->Fill(betaGen);
	}
	if (vPos.size() >= 1 ) 
	{
		fHbeta_tot->Fill(betaGen);
		fHeta_tot->Fill(hscpDir.Eta());
	}
  if (vPos.size() > 2)	
  {
   	params = doFit(vPos,vBx);
   	cout << "Fit parameters: " << params[0] << " " << params[1] << " " << params[2] << " " << params[3] << endl;
   	fHt0->Fill(params[0]);
   	double betaRPC = 1/((params[2]*c)+1);
   	if (params[3]/params[2] < 0.3 && params[2]>0 ) //
   	{
   		fHbeta->Fill(betaRPC);
   		double res = (betaGen-betaRPC)/betaGen;
   		fHres->Fill(res);
   		fHbeta_pas->Fill(betaGen);
   	}
  }
  fHnHits->Fill(vPos.size());
}


// ------------ method called once each job just before starting event loop  ------------
void 
HSCPRecHits::beginJob()
{
   Service<TFileService> fs;
	fHnHits = fs->make<TH1D>("fHnHits", "n hits in sTau direction", 20, 0., 20. );
	fHt0 = fs->make<TH1D>("fHt0","initial time",100,-25.,25.);
	fHbeta = fs->make<TH1D>("fHbeta","beta",150,-1.,2.);
	fHeta_tot = fs->make<TH1D>("fHeta_tot","eta",50,-3.15,3.15);
	fHbeta0 = fs->make<TH1D>("fHbeta0","beta for stau without rechits",50,0.,1.);
	fHeta0 = fs->make<TH1D>("fHeta0","eta for stau without rechits",50,-3.15,3.15);
	fHbeta_pas = fs->make<TH1D>("fHbeta_pas","beta generated pas",50,0.,1.);
	fHbeta_tot = fs->make<TH1D>("fHbeta_tot","beta generated tot",50,0.,1.);
	fHres = fs->make<TH1D>("fHres","Beta resolution",60,-3.,3.);
	fHdR = fs->make<TH1D>("fHdR","fHdR",50,0.,5.);
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
