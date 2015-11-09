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

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include "TTree.h"
#include "Math/VectorUtil.h"
#include "TMath.h"

#include "TreeParticles.hpp"

typedef std::vector<PileupSummaryInfo> PileupSummaryInfoCollection;
typedef std::vector<pat::PackedCandidate> PackedCandidateCollection;
typedef std::vector<tree::PFCandidate> PFCandidateCollection;

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
   double dR_leadingJet_gen_reco_cut_;

   edm::EDGetTokenT<reco::VertexCollection>    vtxToken_;
   edm::EDGetTokenT<PackedCandidateCollection>      packedCandidateToken_;
//   packedCandidateToken_   (consumes<pat::PackedCandidate>(iConfig.getParameter<edm::InputTag>("packedPFCandidates")))
   
   edm::EDGetTokenT<edm::View<pat::Photon> >   photonCollectionToken_;
   edm::EDGetTokenT<pat::JetCollection>        jetCollectionToken_;
   edm::EDGetTokenT<reco::GenJetCollection>    genJetCollectionToken_;
   edm::EDGetTokenT<pat::MuonCollection>       muonCollectionToken_;
   edm::EDGetTokenT<edm::View<pat::Electron> > electronCollectionToken_;
   edm::EDGetTokenT<pat::METCollection>        metCollectionToken_;
   edm::EDGetTokenT<double>                    rhoToken_;
   edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
   edm::EDGetTokenT<PileupSummaryInfoCollection>  pileUpSummaryToken_;
   edm::EDGetTokenT<LHEEventProduct>           LHEEventToken_;

   // electron id
   edm::EDGetTokenT<edm::ValueMap<bool> > electronVetoIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > electronLooseIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > electronMediumIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > electronTightIdMapToken_;
   // photon id
   edm::EDGetTokenT<edm::ValueMap<bool> > photonLooseIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > photonMediumIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > photonTightIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<float>> photonMvaValuesMapToken_;

   // met filters to apply
   const std::vector<std::string> metFilterNames_;

   // from photon ID value map producer
   edm::EDGetTokenT<edm::ValueMap<float> > phoWorstChargedIsolationToken_;

   const std::string pileupHistogramName_;
   edm::EDGetTokenT<bool> HBHENoiseFilterResult_;

   // === TREE DATA ===
   TTree *eventTree_;

   Int_t   nGoodVertices_;
   Int_t   nChargedPfCandidates_;
   
   Float_t rho_;   // the rho variable

   Float_t pu_weight_; // pileup weight
   Char_t mc_weight_; // True for positive event weights

   Float_t dummyFloat_=0.;
   Float_t dR_recoGenJet_=-1;

   Int_t genLeptonsFromW_=0.; // used e.g. to categrize ttbar decays
   Float_t genHt_;

   ULong64_t evtNo_;
   UInt_t    runNo_;
   UInt_t    lumNo_;

   // Trigger decisions
   std::vector<std::string> triggerNames_;
   std::map<std::string, Bool_t > triggerDecision_;
   std::map<std::string, int > triggerIndex_;

   // physics Objects

   PFCandidateCollection vUnpackedPFCandidates_;
   
   PackedCandidateCollection vPackedPFCandidates_;
   
   std::vector<tree::Photon>   vPhotons_;
   std::vector<tree::Jet>      vJets_;
   std::vector<tree::Particle> vGenJets_;
   std::vector<tree::Electron> vElectrons_;
   std::vector<tree::Muon>     vMuons_;
   tree::MET                   met_;
   std::vector<tree::Particle> vGenPhotons_;

   std::vector<tree::GenParticle> vGenParticles_;
   // === ========  ===

   // histogram to store #evts after each "cut"
   TH1F* hCutFlow_;

   // Pileup histogram(s)
   TH1F hPU_;

};


#endif /* TREEWRITER_HPP__ */
