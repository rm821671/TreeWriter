#ifndef TREEWRITER_HPP__
#define TREEWRITER_HPP__

// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include "TTree.h"
#include "Math/VectorUtil.h"
#include "TMath.h"

#include "TreeParticles.hpp"



//
// class declaration
//

class TreeWriter : public edm::EDAnalyzer {
public:
   explicit TreeWriter(const edm::ParameterSet&);
   ~TreeWriter();

   static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   enum PhotonMatchType {UNMATCHED = 0, 
			 MATCHED_FROM_GUDSCB,
			 MATCHED_FROM_PI0,
			 MATCHED_FROM_OTHER_SOURCES};
  
private:
   virtual void beginJob() override;
   virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
   virtual void endJob() override;

   //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
   //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
   //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
   //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

   int matchToTruth(const pat::Photon &pho, 
		    const edm::Handle<edm::View<reco::GenParticle>>  &genParticles);
   int matchToTruthAlternative(const pat::Photon &pho, 
			       const edm::Handle<edm::View<reco::GenParticle>>  &genParticles);

   void findFirstNonPhotonMother(const reco::Candidate *particle,
				 int &ancestorPID, int &ancestorStatus);

   // ----------member data ---------------------------
   double dHT_cut_;
   double dPhoton_pT_cut_;

   edm::EDGetTokenT<reco::VertexCollection>    vtxToken_;
   edm::EDGetTokenT<edm::View<pat::Photon> >   photonCollectionToken_;
   edm::EDGetTokenT<pat::JetCollection>        jetCollectionToken_;
   edm::EDGetTokenT<reco::GenJetCollection>    genJetCollectionToken_;
   edm::EDGetTokenT<pat::MuonCollection>       muonCollectionToken_;
   edm::EDGetTokenT<edm::View<pat::Electron> > electronCollectionToken_;
   edm::EDGetTokenT<pat::METCollection>        metCollectionToken_;
   edm::EDGetTokenT<double>                    rhoToken_;
   edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
   // Value maps with various quantities produced upstream (for photon id)
   edm::EDGetTokenT<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMapToken_; 
   edm::EDGetTokenT<edm::ValueMap<float> > full5x5SigmaIEtaIPhiMapToken_; 
   edm::EDGetTokenT<edm::ValueMap<float> > full5x5E1x3MapToken_; 
   edm::EDGetTokenT<edm::ValueMap<float> > full5x5E2x2MapToken_; 
   edm::EDGetTokenT<edm::ValueMap<float> > full5x5E2x5MaxMapToken_; 
   edm::EDGetTokenT<edm::ValueMap<float> > full5x5E5x5MapToken_; 
   edm::EDGetTokenT<edm::ValueMap<float> > esEffSigmaRRMapToken_; 
   //
   edm::EDGetTokenT<edm::ValueMap<float> > phoChargedIsolationToken_; 
   edm::EDGetTokenT<edm::ValueMap<float> > phoNeutralHadronIsolationToken_; 
   edm::EDGetTokenT<edm::ValueMap<float> > phoPhotonIsolationToken_; 
   edm::EDGetTokenT<edm::ValueMap<float> > phoWorstChargedIsolationToken_; 

   // electron id
   edm::EDGetTokenT<edm::ValueMap<bool> > electronVetoIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > electronLooseIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > electronMediumIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > electronTightIdMapToken_;
   // photon id
   edm::EDGetTokenT<edm::ValueMap<bool> > photonLooseIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > photonMediumIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > photonTightIdMapToken_;

   // met filters to apply
   const std::vector<std::string> metFilterNames_;

   const std::string pileupHistogramName_;

   // === TREE DATA ===
   TTree *eventTree_;

   Bool_t  isRealData_; // whether data or MC
   Int_t   nPV_;   // number of reconsrtucted primary vertices
   Int_t   true_nPV_;   // true number of reconsrtucted primary vertices
   Float_t pu_weight; // pileup weight
   Int_t   nGoodVertices_;
   Float_t rho_;   // the rho variable

   ULong64_t evtNo_;
   UInt_t    runNo_;
   UInt_t    lumNo_;

   // physics Objects
   std::vector<tree::Photon>   vPhotons_;
   std::vector<tree::Jet>      vJets_;
   std::vector<tree::Particle> vGenJets_;
   std::vector<tree::Electron> vElectrons_;
   std::vector<tree::Muon>     vMuons_;
   tree::MET                   met_;
   std::vector<tree::Particle> vGenPhotons_;
   std::vector<tree::Particle> vGenElectrons_;
   // === ========  ===

   // histogram to store #evts after each "cut"
   TH1F* hCutFlow_;


   // === VARIABLES NOT STORED IN TREE ===
   // Variables that will be containers on which TMVA Reader works
   // The variables
   float varPhi_;
   float varR9_; 
   float varSieie_;
   float varSieip_; 
   float varE1x3overE5x5_; 
   float varE2x2overE5x5_; 
   float varE2x5overE5x5_; 
   float varSCEta_; 
   float varRawE_; 
   float varSCEtaWidth_; 
   float varSCPhiWidth_; 
   float varRho_;
   float varPhoIsoRaw_;
   float varChIsoRaw_; 
   float varWorstChRaw_;
   float varESEnOverRawE_; // for endcap MVA only
   float varESEffSigmaRR_; // for endcap MVA only
   // The spectators
   float varPt_; 
   float varEta_;

   // TMVA Reader for applying MVA
   TMVA::Reader *tmvaReader_[2];
   TString methodName_[2];

   // Pileup histogram(s)
   TH1F hPU_;

};


#endif /* TREEWRITER_HPP__ */
