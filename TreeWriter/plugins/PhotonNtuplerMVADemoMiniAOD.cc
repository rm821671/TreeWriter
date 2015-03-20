// -*- C++ -*-
//
// Package:    ElectronWork/ElectronNtupler
// Class:      PhotonNtuplerMVADemoMiniAOD
// 
/**\class PhotonNtuplerMVADemoMiniAOD PhotonNtuplerMVADemoMiniAOD.cc ElectronWork/ElectornNtupler/plugins/PhotonNtuplerMVADemoMiniAOD.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Ilya Kravchenko
//         Created:  Thu, 10 Jul 2014 09:54:13 GMT
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include "TTree.h"
#include "Math/VectorUtil.h"
#include "TMath.h"


//
// class declaration
//

class PhotonNtuplerMVADemoMiniAOD : public edm::EDAnalyzer {
public:
  explicit PhotonNtuplerMVADemoMiniAOD(const edm::ParameterSet&);
  ~PhotonNtuplerMVADemoMiniAOD();

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
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<edm::View<pat::Photon> > photonCollectionToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
  // Value maps with various quantities produced upstream
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

  TTree *photonTree_;

  Int_t nPV_;        // number of reconsrtucted primary vertices
  Float_t rho_;      // the rho variable


  // all photon variables
  Int_t nPhotons_;

  std::vector<Float_t> pt_;
  std::vector<Float_t> eta_;
  std::vector<Float_t> phi_;

  // Variables for cut based ID
  std::vector<Float_t> full5x5_sigmaIetaIeta_;
  std::vector<Float_t> hOverE_;
  std::vector<Int_t> hasPixelSeed_;
  std::vector<Int_t> passElectronVeto_;

  std::vector<Float_t> isoChargedHadrons_;
  std::vector<Float_t> isoNeutralHadrons_;
  std::vector<Float_t> isoPhotons_;

  std::vector<Float_t> isoChargedHadronsWithEA_;
  std::vector<Float_t> isoNeutralHadronsWithEA_;
  std::vector<Float_t> isoPhotonsWithEA_;

  // Extra variables for MVA ID (excluding already mentioned above)
  std::vector<Float_t> scRawEnergy_;
  std::vector<Float_t> isoWorstChargedHadrons_;
  std::vector<Float_t> r9_;
  std::vector<Float_t> full5x5_sigmaIetaIphi_;
  std::vector<Float_t> full5x5_e1x3_;
  std::vector<Float_t> full5x5_e2x2_;
  std::vector<Float_t> full5x5_e2x5Max_;
  std::vector<Float_t> full5x5_e5x5_;
  std::vector<Float_t> sigma_eta_;
  std::vector<Float_t> sigma_phi_;
  std::vector<Float_t> esEffSigmaRR_;
  std::vector<Float_t> esEnergy_;

  std::vector<Float_t> mvaValue_;
  
  std::vector<Int_t> isTrue_;
  std::vector<Int_t> isTrueAlternative_;

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

};

//
// constants, enums and typedefs
//

// Effective areas for photons from Savvas's slides
// for phys14 PU20bx25, described here:
// https://indico.cern.ch/event/367861/contribution/3/material/slides/0.pdf
namespace EffectiveAreas {
  const int nEtaBins = 7;
  const float etaBinLimits[nEtaBins+1] = {
    0.0, 1.0, 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};

  const float areaPhotons[nEtaBins] = {
    0.0894, 0.0750, 0.0423, 0.0561, 0.0882, 0.1144, 0.1684
  };
  const float areaNeutralHadrons[nEtaBins] = {
    0.049, 0.0108, 0.0019, 0.0037, 0.0062, 0.0130, 0.1699
  };
  const float areaChargedHadrons[nEtaBins] = {
    0.0089, 0.0062, 0.0086, 0.0041, 0.0113, 0.0085, 0.0039
  };
}
//


//
// static data member definitions
//

//
// constructors and destructor
//
PhotonNtuplerMVADemoMiniAOD::PhotonNtuplerMVADemoMiniAOD(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  photonCollectionToken_(consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photons"))),
  rhoToken_(consumes<double> (iConfig.getParameter<edm::InputTag>("rho"))),
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
  // Cluster shapes
  full5x5SigmaIEtaIEtaMapToken_(consumes <edm::ValueMap<float> >
				(iConfig.getParameter<edm::InputTag>("full5x5SigmaIEtaIEtaMap"))),
  full5x5SigmaIEtaIPhiMapToken_(consumes <edm::ValueMap<float> >
				(iConfig.getParameter<edm::InputTag>("full5x5SigmaIEtaIPhiMap"))),
  full5x5E1x3MapToken_(consumes <edm::ValueMap<float> >
		       (iConfig.getParameter<edm::InputTag>("full5x5E1x3Map"))),
  full5x5E2x2MapToken_(consumes <edm::ValueMap<float> >
		       (iConfig.getParameter<edm::InputTag>("full5x5E2x2Map"))),
  full5x5E2x5MaxMapToken_(consumes <edm::ValueMap<float> >
			  (iConfig.getParameter<edm::InputTag>("full5x5E2x5MaxMap"))),
  full5x5E5x5MapToken_(consumes <edm::ValueMap<float> >
		       (iConfig.getParameter<edm::InputTag>("full5x5E5x5Map"))),
  esEffSigmaRRMapToken_(consumes <edm::ValueMap<float> >
			(iConfig.getParameter<edm::InputTag>("esEffSigmaRRMap"))),
  // Isolations
  phoChargedIsolationToken_(consumes <edm::ValueMap<float> >
			    (iConfig.getParameter<edm::InputTag>("phoChargedIsolation"))),
  phoNeutralHadronIsolationToken_(consumes <edm::ValueMap<float> >
				  (iConfig.getParameter<edm::InputTag>("phoNeutralHadronIsolation"))),
  phoPhotonIsolationToken_(consumes <edm::ValueMap<float> >
			   (iConfig.getParameter<edm::InputTag>("phoPhotonIsolation"))),
  phoWorstChargedIsolationToken_(consumes <edm::ValueMap<float> >
				 (iConfig.getParameter<edm::InputTag>("phoWorstChargedIsolation")))
{

  edm::Service<TFileService> fs;
  photonTree_ = fs->make<TTree> ("PhotonTree", "Photon data");
  
  photonTree_->Branch("nPV"        ,  &nPV_     , "nPV/I");
  photonTree_->Branch("rho"        ,  &rho_ , "rho/F");
  photonTree_->Branch("nPho",  &nPhotons_ , "nPho/I");

  photonTree_->Branch("pt"           ,  &pt_    );
  photonTree_->Branch("scRawEnergy"  ,  &scRawEnergy_    );
  photonTree_->Branch("esEnergy"     , &esEnergy_);
  photonTree_->Branch("eta"          ,  &eta_ );
  photonTree_->Branch("phi"          ,  &phi_ );

  photonTree_->Branch("hasPixelSeed"           ,  &hasPixelSeed_);
  photonTree_->Branch("passElectronVeto"       ,  &passElectronVeto_);
  photonTree_->Branch("hOverE"                 ,  &hOverE_);

  // Cluster shapes
  photonTree_->Branch("r9"  , &r9_);
  photonTree_->Branch("full5x5_sigmaIetaIeta"  , &full5x5_sigmaIetaIeta_);
  photonTree_->Branch("full5x5_sigmaIetaIphi"  , &full5x5_sigmaIetaIphi_);
  photonTree_->Branch("full5x5_e1x3"   , &full5x5_e1x3_);
  photonTree_->Branch("full5x5_e2x2"   , &full5x5_e2x2_);
  photonTree_->Branch("full5x5_e2x5Max", &full5x5_e2x5Max_);
  photonTree_->Branch("full5x5_e5x5"   , &full5x5_e5x5_);
  photonTree_->Branch("sigma_eta"      , &sigma_eta_);
  photonTree_->Branch("sigma_phi"      , &sigma_phi_);
  photonTree_->Branch("esEffSigmaRR"  , &esEffSigmaRR_);

  // Isolations
  photonTree_->Branch("isoChargedHadrons"      , &isoChargedHadrons_);
  photonTree_->Branch("isoNeutralHadrons"      , &isoNeutralHadrons_);
  photonTree_->Branch("isoPhotons"             , &isoPhotons_);
  photonTree_->Branch("isoWorstChargedHadrons" , &isoWorstChargedHadrons_);

  photonTree_->Branch("isoChargedHadronsWithEA"      , &isoChargedHadronsWithEA_);
  photonTree_->Branch("isoNeutralHadronsWithEA"      , &isoNeutralHadronsWithEA_);
  photonTree_->Branch("isoPhotonsWithEA"             , &isoPhotonsWithEA_);

  photonTree_->Branch("mvaValue"           , &mvaValue_);

  photonTree_->Branch("isTrue"             , &isTrue_);
  photonTree_->Branch("isTrueAlternative"  , &isTrueAlternative_);
 
  //
  // Create and configure barrel MVA
  //
  tmvaReader_[0] = new TMVA::Reader( "!Color:!Silent:Error" );  
  tmvaReader_[0]->SetVerbose(kFALSE);
  // Add all the vars, we take the string with variable name from the weights file (the Expression field)
  tmvaReader_[0]->AddVariable("recoPhi"   , &varPhi_);
  tmvaReader_[0]->AddVariable("r9"        , &varR9_);
  tmvaReader_[0]->AddVariable("sieie_2012", &varSieie_);
  tmvaReader_[0]->AddVariable("sieip_2012", &varSieip_);
  tmvaReader_[0]->AddVariable("e1x3_2012/e5x5_2012"        , &varE1x3overE5x5_);
  tmvaReader_[0]->AddVariable("e2x2_2012/e5x5_2012"        , &varE2x2overE5x5_);
  tmvaReader_[0]->AddVariable("e2x5_2012/e5x5_2012"        , &varE2x5overE5x5_);
  tmvaReader_[0]->AddVariable("recoSCEta" , &varSCEta_);
  tmvaReader_[0]->AddVariable("rawE"      , &varRawE_);
  tmvaReader_[0]->AddVariable("scEtaWidth", &varSCEtaWidth_);
  tmvaReader_[0]->AddVariable("scPhiWidth", &varSCPhiWidth_);
  tmvaReader_[0]->AddVariable("rho"       , &varRho_);
  tmvaReader_[0]->AddVariable("phoIsoRaw" , &varPhoIsoRaw_);
  tmvaReader_[0]->AddVariable("chIsoRaw"  , &varChIsoRaw_);
  tmvaReader_[0]->AddVariable("chWorstRaw", &varWorstChRaw_);
  // Add spectators
  tmvaReader_[0]->AddSpectator("recoPt" , &varPt_);
  tmvaReader_[0]->AddSpectator("recoEta", &varEta_);

  //
  // Create and configure endcap MVA
  //
  tmvaReader_[1] = new TMVA::Reader( "!Color:!Silent:Error" );  
  tmvaReader_[1]->SetVerbose(kFALSE);
  // Add all the vars, we take the string with variable name from the weights file (the Expression field)
  tmvaReader_[1]->AddVariable("recoPhi"   , &varPhi_);
  tmvaReader_[1]->AddVariable("r9"        , &varR9_);
  tmvaReader_[1]->AddVariable("sieie_2012", &varSieie_);
  tmvaReader_[1]->AddVariable("sieip_2012", &varSieip_);
  tmvaReader_[1]->AddVariable("e1x3_2012/e5x5_2012"        , &varE1x3overE5x5_);
  tmvaReader_[1]->AddVariable("e2x2_2012/e5x5_2012"        , &varE2x2overE5x5_);
  tmvaReader_[1]->AddVariable("e2x5_2012/e5x5_2012"        , &varE2x5overE5x5_);
  tmvaReader_[1]->AddVariable("recoSCEta" , &varSCEta_);
  tmvaReader_[1]->AddVariable("rawE"      , &varRawE_);
  tmvaReader_[1]->AddVariable("scEtaWidth", &varSCEtaWidth_);
  tmvaReader_[1]->AddVariable("scPhiWidth", &varSCPhiWidth_);
  tmvaReader_[1]->AddVariable("esEn/rawE" , &varESEnOverRawE_);
  tmvaReader_[1]->AddVariable("esRR"      , &varESEffSigmaRR_);
  tmvaReader_[1]->AddVariable("rho"       , &varRho_);
  tmvaReader_[1]->AddVariable("phoIsoRaw" , &varPhoIsoRaw_);
  tmvaReader_[1]->AddVariable("chIsoRaw"  , &varChIsoRaw_);
  tmvaReader_[1]->AddVariable("chWorstRaw", &varWorstChRaw_);
  // Add spectators
  tmvaReader_[1]->AddSpectator("recoPt" , &varPt_);
  tmvaReader_[1]->AddSpectator("recoEta", &varEta_);

  //
  // Book the MVA method for each category
  //
  std::string cmssw_base_src = getenv("CMSSW_BASE");
  cmssw_base_src += "/src/";
  //
  TString localFileName1 = "EgammaAnalysis/PhotonTools/data/PHYS14/photon_general_MVA_phys14_pu20bx25_EB_V1.weights.xml";
  TString weightsFileName1 = TString(cmssw_base_src) + localFileName1;
  methodName_[0] = "BDTG photons barrel";
  tmvaReader_[0]->BookMVA(methodName_[0], weightsFileName1);
  //
  TString localFileName2 = "EgammaAnalysis/PhotonTools/data/PHYS14/photon_general_MVA_phys14_pu20bx25_EE_V1.weights.xml";
  TString weightsFileName2 = TString(cmssw_base_src) + localFileName2;
  methodName_[1] = "BDTG photons endcap";
  tmvaReader_[1]->BookMVA(methodName_[1], weightsFileName2);

}


PhotonNtuplerMVADemoMiniAOD::~PhotonNtuplerMVADemoMiniAOD()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  delete tmvaReader_[0];
  delete tmvaReader_[1];

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PhotonNtuplerMVADemoMiniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
  
  // // An object needed for isolation calculations
  // GEDPhoIDTools *GEDIdTool = new GEDPhoIDTools(iEvent);

  // Get photon collection
  edm::Handle<edm::View<pat::Photon> > collection;
  iEvent.getByToken(photonCollectionToken_, collection);

  // Get PV
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  //const reco::Vertex &pv = vertices->front();
  nPV_    = vertices->size();

  VertexCollection::const_iterator firstGoodVertex = vertices->end();
  int firstGoodVertexIdx = 0;
  for (VertexCollection::const_iterator vtx = vertices->begin(); 
       vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
    // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
    // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
    if (  /*!vtx->isFake() &&*/ 
        !(vtx->chi2()==0 && vtx->ndof()==0) 
        &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
        && fabs(vtx->position().Z())<=24.0) {
      firstGoodVertex = vtx;
      break;
    }
  }

  if ( firstGoodVertex==vertices->end() )
    return; // skip event if there are no good PVs

  // Get rho
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho_ = *rhoH;

  // Get generator level info
  // Pruned particles are the one containing "important" stuff
  Handle<edm::View<reco::GenParticle> > prunedGenParticles;
  iEvent.getByToken(prunedGenToken_,prunedGenParticles);

  // Get the full5x5 maps
  edm::Handle<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMap;
  iEvent.getByToken(full5x5SigmaIEtaIEtaMapToken_, full5x5SigmaIEtaIEtaMap);
  edm::Handle<edm::ValueMap<float> > full5x5SigmaIEtaIPhiMap;
  iEvent.getByToken(full5x5SigmaIEtaIPhiMapToken_, full5x5SigmaIEtaIPhiMap);

  edm::Handle<edm::ValueMap<float> > full5x5E1x3Map;
  iEvent.getByToken(full5x5E1x3MapToken_, full5x5E1x3Map);

  edm::Handle<edm::ValueMap<float> > full5x5E2x2Map;
  iEvent.getByToken(full5x5E2x2MapToken_, full5x5E2x2Map);

  edm::Handle<edm::ValueMap<float> > full5x5E2x5MaxMap;
  iEvent.getByToken(full5x5E2x5MaxMapToken_, full5x5E2x5MaxMap);

  edm::Handle<edm::ValueMap<float> > full5x5E5x5Map;
  iEvent.getByToken(full5x5E5x5MapToken_, full5x5E5x5Map);

  edm::Handle<edm::ValueMap<float> > esEffSigmaRRMap;
  iEvent.getByToken(esEffSigmaRRMapToken_, esEffSigmaRRMap);

  // Get the isolation maps
  edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
  iEvent.getByToken(phoChargedIsolationToken_, phoChargedIsolationMap);
  edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
  iEvent.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
  edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
  iEvent.getByToken(phoPhotonIsolationToken_, phoPhotonIsolationMap);
  edm::Handle<edm::ValueMap<float> > phoWorstChargedIsolationMap;
  iEvent.getByToken(phoWorstChargedIsolationToken_, phoWorstChargedIsolationMap);


  // Clear vectors
  nPhotons_ = 0;
  pt_.clear();
  eta_.clear();
  phi_.clear();
  full5x5_sigmaIetaIeta_.clear();
  hOverE_.clear();
  hasPixelSeed_.clear();
  passElectronVeto_.clear();
  isoChargedHadrons_.clear();
  isoNeutralHadrons_.clear();
  isoPhotons_.clear();
  isoChargedHadronsWithEA_.clear();
  isoNeutralHadronsWithEA_.clear();
  isoPhotonsWithEA_.clear();
  scRawEnergy_.clear();
  isoWorstChargedHadrons_.clear();
  r9_.clear();
  full5x5_sigmaIetaIphi_.clear();
  full5x5_e1x3_.clear();
  full5x5_e2x2_.clear();
  full5x5_e2x5Max_.clear();
  full5x5_e5x5_.clear();
  sigma_eta_.clear();
  sigma_phi_.clear();
  esEffSigmaRR_.clear();
  esEnergy_.clear();
  mvaValue_.clear();
  isTrue_.clear();
  isTrueAlternative_.clear();

  // Loop over photons
  
  // const auto& pho_refs = collection->refVector();
  
  // for( const auto& pho : pho_refs ) {

  for( View<pat::Photon>::const_iterator pho = collection->begin();
       pho != collection->end(); pho++){
    
    // Kinematics
    if( pho->pt() < 15 ) 
      continue;
    
    nPhotons_++;
    pt_          .push_back( pho->pt() );
    scRawEnergy_ .push_back( pho->superCluster()->rawEnergy() );
    esEnergy_    .push_back( pho->superCluster()->preshowerEnergy() );
    eta_         .push_back( pho->superCluster()->eta() );
    phi_         .push_back( pho->superCluster()->phi() );

    //const edm::Ptr<pat::Photon> phoPtr( pho );
    const edm::Ptr<pat::Photon> phoPtr( collection, pho - collection->begin() );

    hOverE_                .push_back( pho->hadTowOverEm() );
    hasPixelSeed_          .push_back( (Int_t)pho->hasPixelSeed() );

    passElectronVeto_      .push_back( pho->passElectronVeto() );

    sigma_eta_             .push_back( pho->superCluster()->etaWidth() );
    sigma_phi_             .push_back( pho->superCluster()->phiWidth() );
    r9_                    .push_back( pho->r9() );
    full5x5_sigmaIetaIeta_ .push_back( (*full5x5SigmaIEtaIEtaMap)[ phoPtr ] );
    full5x5_sigmaIetaIphi_ .push_back( (*full5x5SigmaIEtaIPhiMap)[ phoPtr ] );

    full5x5_e1x3_    .push_back( (*full5x5E1x3Map)[ phoPtr ] );
    full5x5_e2x2_    .push_back( (*full5x5E2x2Map)[ phoPtr ] );
    full5x5_e2x5Max_ .push_back( (*full5x5E2x5MaxMap)[ phoPtr ] );
    full5x5_e5x5_    .push_back( (*full5x5E5x5Map)[ phoPtr ] );
    esEffSigmaRR_    .push_back( (*esEffSigmaRRMap)[ phoPtr ] );

    isoChargedHadrons_ .push_back( (*phoChargedIsolationMap)[phoPtr] );
    isoNeutralHadrons_ .push_back( (*phoNeutralHadronIsolationMap)[phoPtr] );
    isoPhotons_        .push_back( (*phoPhotonIsolationMap)[phoPtr] );
    isoWorstChargedHadrons_ .push_back( (*phoWorstChargedIsolationMap)[phoPtr] );

    // Compute isolation with effective area correction for PU
    // Find eta bin first. If eta>2.5, the last eta bin is used.
    int etaBin = 0; 
    while ( etaBin < EffectiveAreas::nEtaBins-1 
	    && abs( pho->superCluster()->eta() ) > EffectiveAreas::etaBinLimits[etaBin+1] )
      { ++etaBin; };
    isoPhotonsWithEA_        .push_back( std::max( (float)0.0, (*phoPhotonIsolationMap)       [phoPtr] 
						   - rho_ * EffectiveAreas::areaPhotons[etaBin] ) );
    isoNeutralHadronsWithEA_ .push_back( std::max( (float)0.0, (*phoNeutralHadronIsolationMap)[phoPtr] 
						   - rho_ * EffectiveAreas::areaNeutralHadrons[etaBin] ) );
    isoChargedHadronsWithEA_ .push_back( std::max( (float)0.0, (*phoChargedIsolationMap)      [phoPtr] 
						   - rho_ * EffectiveAreas::areaChargedHadrons[etaBin] ) );

    // Prepare variables and find the MVA value
    varPhi_          = pho->phi();
    varR9_           = pho->r9() ;
    varSieie_        = (*full5x5SigmaIEtaIEtaMap)[ phoPtr ];
    varSieip_        = (*full5x5SigmaIEtaIPhiMap)[ phoPtr ];
    float e5x5 = (*full5x5E5x5Map)[ phoPtr ];
    // Protect from e5x5 being zero since in miniAOD not the full info is stored
    // for the poor quality photons.
    if( e5x5 != 0 ){
      varE1x3overE5x5_ = (*full5x5E1x3Map)[ phoPtr ] / e5x5;
      varE2x2overE5x5_ = (*full5x5E2x2Map)[ phoPtr ] / e5x5;
      varE2x5overE5x5_ = (*full5x5E2x5MaxMap)[ phoPtr ]/ e5x5;
    }else{
      varE1x3overE5x5_ = 0;
      varE2x2overE5x5_ = 0;
      varE2x5overE5x5_ = 0;
    }
    varSCEta_        = pho->superCluster()->eta(); 
    varRawE_         = pho->superCluster()->rawEnergy(); 
    varSCEtaWidth_   = pho->superCluster()->etaWidth(); 
    varSCPhiWidth_   = pho->superCluster()->phiWidth(); 
    varESEnOverRawE_ = pho->superCluster()->preshowerEnergy() / pho->superCluster()->rawEnergy();
    varESEffSigmaRR_ = (*esEffSigmaRRMap)[ phoPtr ];
    varRho_          = rho_; 
    varPhoIsoRaw_    = (*phoPhotonIsolationMap)[phoPtr];  
    varChIsoRaw_     = (*phoChargedIsolationMap)[phoPtr];
    varWorstChRaw_   = (*phoWorstChargedIsolationMap)[phoPtr];
    // Declare spectator vars
    varPt_ = pho->pt(); 
    varEta_ = pho->eta();

    // DEBUG
    const bool debug = false;
    if( debug && fabs( pho->superCluster()->eta())<1.479 ){
      printf("Printout of barrel electron variables for MVA:\n");
      printf("  varPhi_           %f\n", varPhi_         );
      printf("  varR9_            %f\n", varR9_          ); 
      printf("  varSieie_         %f\n", varSieie_       );
      printf("  varSieip_         %f\n", varSieip_       ); 
      printf("  varE1x3overE5x5_  %f\n", varE1x3overE5x5_); 
      printf("  varE2x2overE5x5_  %f\n", varE2x2overE5x5_); 
      printf("  varE2x5overE5x5_  %f\n", varE2x5overE5x5_); 
      printf("  varSCEta_         %f\n", varSCEta_       ); 
      printf("  varRawE_          %f\n", varRawE_        ); 
      printf("  varSCEtaWidth_    %f\n", varSCEtaWidth_  ); 
      printf("  varSCPhiWidth_    %f\n", varSCPhiWidth_  ); 
      printf("  varRho_           %f\n", varRho_         );
      printf("  varPhoIsoRaw_     %f\n", varPhoIsoRaw_   );
      printf("  varChIsoRaw_      %f\n", varChIsoRaw_    ); 
      printf("  varWorstChRaw_    %f\n", varWorstChRaw_  );
      printf("  varESEnOverRawE_  %f\n", varESEnOverRawE_); // for endcap MVA only
      printf("  varESEffSigmaRR_  %f\n", varESEffSigmaRR_); // for endcap MVA only
      // The spectators
      printf("  varPt_    %f\n", varPt_          ); 
      printf("  varEta_  %f\n", varEta_         );
    }

    //
    // Compute the MVA value for this photon. The MVA value here is stored
    // in a TTree, but one can also cut on the MVA value at this point.
    //
    if( e5x5 != 0 ){
      if( abs( pho->superCluster()->eta() ) < 1.479 )
	mvaValue_ .push_back( tmvaReader_[0]->EvaluateMVA(methodName_[0]) );
      else
	mvaValue_ .push_back( tmvaReader_[1]->EvaluateMVA(methodName_[1]) );
    }else{
      // e5x5 zero means that this photon's info hasn't been stored fully in
      // miniAOD since it is a poor quality photon. We can't run MVA on it.
      mvaValue_.push_back( -999. );
    }

    // MC match
    isTrue_.push_back( matchToTruth(*pho, prunedGenParticles) );
    isTrueAlternative_.push_back( matchToTruthAlternative(*pho, prunedGenParticles) );

  }
   
  // Save the info
  photonTree_->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
PhotonNtuplerMVADemoMiniAOD::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhotonNtuplerMVADemoMiniAOD::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
  void 
  PhotonNtuplerMVADemoMiniAOD::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void 
  PhotonNtuplerMVADemoMiniAOD::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
  void 
  PhotonNtuplerMVADemoMiniAOD::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void 
  PhotonNtuplerMVADemoMiniAOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhotonNtuplerMVADemoMiniAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

int PhotonNtuplerMVADemoMiniAOD::matchToTruth(const pat::Photon &pho, 
					      const edm::Handle<edm::View<reco::GenParticle>>  
					      &genParticles)
{
  // 
  // Explicit loop and geometric matching method 
  //

  // Find the closest status 1 gen photon to the reco photon
  double dR = 999;
  const reco::Candidate *closestPhoton = 0;
  for(size_t i=0; i<genParticles->size();i++){
    const reco::Candidate *particle = &(*genParticles)[i];
    // Drop everything that is not photon or not status 1
    if( abs(particle->pdgId()) != 22 || particle->status() != 1 )
      continue;
    //
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( pho.p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestPhoton = particle;
    }
  }
  // See if the closest photon (if it exists) is close enough.
  // If not, no match found.
  if( !(closestPhoton != 0 && dR < 0.1) ) {
    return UNMATCHED;
  }

  // Find ID of the parent of the found generator level photon match
  int ancestorPID = -999; 
  int ancestorStatus = -999;
  findFirstNonPhotonMother(closestPhoton, ancestorPID, ancestorStatus);

  // Allowed parens: quarks pdgId 1-5, or a gluon 21
  std::vector<int> allowedParents { -1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -21, 21 };
  if( !(std::find(allowedParents.begin(), 
		  allowedParents.end(), ancestorPID)
	!= allowedParents.end()) ){
    // So it is not from g, u, d, s, c, b. Check if it is from pi0 or not. 
    if( abs(ancestorPID) == 111 )
      return MATCHED_FROM_PI0;
    else
      return MATCHED_FROM_OTHER_SOURCES;
  }
  return MATCHED_FROM_GUDSCB;
   
}

void PhotonNtuplerMVADemoMiniAOD::findFirstNonPhotonMother(const reco::Candidate *particle,
							   int &ancestorPID, int &ancestorStatus){
  
  if( particle == 0 ){
    printf("PhotonNtuplerMVADemoMiniAOD: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-photon parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 22 ){
    findFirstNonPhotonMother(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }
  
  return;
}

int PhotonNtuplerMVADemoMiniAOD::matchToTruthAlternative(const pat::Photon &pho, 
							 const edm::Handle<edm::View<reco::GenParticle>>  
							 &genParticles)
{


  // 
  // Explicit loop and geometric matching method 
  //
  
  int isMatched = UNMATCHED;
  
  for(size_t i=0; i<genParticles->size();i++){
    const reco::Candidate *particle = &(*genParticles)[i];
    int pid = particle->pdgId();
    int ancestorPID = -999; 
    int ancestorStatus = -999;
    findFirstNonPhotonMother(particle, ancestorPID, ancestorStatus);
    if( pid ==22 && TMath::Abs( ancestorPID ) <= 22 ){
      double dr = ROOT::Math::VectorUtil::DeltaR( pho.p4(), particle->p4() );
      float dpt = fabs( (pho.pt() - particle->pt() )/particle->pt());
      if (dr < 0.2 && dpt < 0.2){
	isMatched = MATCHED_FROM_GUDSCB;
	if( ancestorPID == 22 ){
	  printf("Ancestor of a photon is a photon!\n");
	}
      }
    }
  }
    
  return isMatched; 
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonNtuplerMVADemoMiniAOD);
