// -*- C++ -*-
//
// Package:    TreeWriter/TreeWriter
// Class:      TreeWriter
//

//
// Original Author:  Johannes Lange (adapted parts from Ilya Kravchenko)
//

#include "TreeWriter.hpp"

// jet ID
static bool isLooseJet(const pat::Jet& jet)
{
   bool pass = true
      && jet.neutralHadronEnergyFraction() < .99
      && jet.neutralEmEnergyFraction()     < .99
      && jet.nConstituents()               > 1
      && jet.muonEnergyFraction()          < .8
      && jet.chargedEmEnergyFraction()     < .9
      ;
   // additional forward requirements:
   if (fabs(jet.eta()) > 2.4){
      pass = pass
         && jet.chargedHadronEnergyFraction() > 0
         && jet.chargedMultiplicity()         > 0
         && jet.chargedEmEnergyFraction()     < .99
         ;
   }
   return pass;
}

// compute HT using RECO objects to "reproduce" the trigger requirements
static double computeHT(const std::vector<tree::Jet>& jets)
{
   double HT=0;
   double pt=0;
   for (const tree::Jet& jet: jets){
      pt=jet.p.Pt();
      if (fabs(jet.p.Eta())<3.0 && pt>40) HT+=pt;
   }
   return HT;
}

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
TreeWriter::TreeWriter(const edm::ParameterSet& iConfig)
   : dHT_cut_(iConfig.getUntrackedParameter<double>("HT_cut"))
   , dPhoton_pT_cut_(iConfig.getUntrackedParameter<double>("photon_pT_cut"))
   , vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices")))
   , photonCollectionToken_  (consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photons")))
   , jetCollectionToken_     (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets")))
   , genJetCollectionToken_  (consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets")))
   , muonCollectionToken_    (consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons")))
   , electronCollectionToken_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons")))
   , metCollectionToken_     (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets")))
   , rhoToken_               (consumes<double> (iConfig.getParameter<edm::InputTag>("rho")))
   , prunedGenToken_         (consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedGenParticles")))
   // electron id
   , electronVetoIdMapToken_  (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronVetoIdMap"   )))
   , electronLooseIdMapToken_ (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronLooseIdMap"  )))
   , electronMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronMediumIdMap" )))
   , electronTightIdMapToken_ (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronTightIdMap"  )))
   // photon id
   , photonLooseIdMapToken_  (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonLooseIdMap"  )))
   , photonMediumIdMapToken_ (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonMediumIdMap" )))
   , photonTightIdMapToken_  (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonTightIdMap"  )))
   , photonMvaValuesMapToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("photonMvaValuesMap")))
   // met filters to apply
   , metFilterNames_(iConfig.getUntrackedParameter<std::vector<std::string>>("metFilterNames"))
   , pileupHistogramName_(iConfig.getUntrackedParameter<std::string>("pileupHistogramName"))
{

   edm::Service<TFileService> fs;
   eventTree_ = fs->make<TTree> ("eventTree", "event data");

   eventTree_->Branch("photons"  , &vPhotons_);
   eventTree_->Branch("jets"     , &vJets_);
   eventTree_->Branch("genJets"  , &vGenJets_);
   eventTree_->Branch("electrons", &vElectrons_);
   eventTree_->Branch("muons"    , &vMuons_);
   eventTree_->Branch("met"      , &met_);
   eventTree_->Branch("genPhotons"   , &vGenPhotons_);
   eventTree_->Branch("genElectrons" , &vGenElectrons_);

   eventTree_->Branch("isRealData"    , &isRealData_    , "isRealData/O");
   eventTree_->Branch("nPV"           , &nPV_           , "nPV/I");
   eventTree_->Branch("true_nPV"      , &true_nPV_      , "true_nPV/I");
   eventTree_->Branch("nGoodVertices" , &nGoodVertices_ , "nGoodVertices/I");
   eventTree_->Branch("rho"           , &rho_           , "rho/F");

   eventTree_->Branch("pu_weight"     , &pu_weight_     , "pu_weight/F");
   eventTree_->Branch("mc_weight"     , &mc_weight_     , "mc_weight/F");

   eventTree_->Branch("dummyFloat" , &dummyFloat_ , "dummyFloat/F");

   eventTree_->Branch("evtNo", &evtNo_, "evtNo/l");
   eventTree_->Branch("runNo", &runNo_, "runNo/i");
   eventTree_->Branch("lumNo", &lumNo_, "lumNo/i");

   // get pileup histogram(s)
   std::string cmssw_base_src = getenv("CMSSW_BASE");
   cmssw_base_src += "/src/";
   TFile puFile(TString(cmssw_base_src+"/TreeWriter/PUreweighting/data/puWeights.root"));
   if (puFile.IsZombie() ){
      edm::LogError("File not found") << "create puWeights.root! (see README)";
      std::exit(84);
   } else {
      hPU_=*( (TH1F*)puFile.Get( pileupHistogramName_.c_str() ) );
   }
   puFile.Close();

   // create cut-flow histogram
   std::vector<TString> vCutBinNames{{"initial","METfilters","nGoodVertices","jets","HT","photons","final"}};
   hCutFlow_ = fs->make<TH1F>("hCutFlow","hCutFlow",vCutBinNames.size(),0,vCutBinNames.size());
   for (uint i=0;i<vCutBinNames.size();i++) hCutFlow_->GetXaxis()->SetBinLabel(i+1,vCutBinNames.at(i));
}


TreeWriter::~TreeWriter(){}

//
// member functions
//

// ------------ method called for each event  ------------
void
TreeWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   hCutFlow_->Fill("initial",1);
   isRealData_=iEvent.isRealData();

   // MET Filters
   edm::Handle<edm::TriggerResults> metFilterBits;
   edm::InputTag metFilterTag("TriggerResults","","PAT");
   iEvent.getByLabel(metFilterTag, metFilterBits);
   // go through the filters and check if they were passed
   const edm::TriggerNames &allFilterNames = iEvent.triggerNames(*metFilterBits);
   for (std::string const &name: metFilterNames_){
      const int index=allFilterNames.triggerIndex(name);
      if (!metFilterBits->accept(index)) return; // not passed
   }
   hCutFlow_->Fill("METfilters",1);

   // Get photon collection
   edm::Handle<edm::View<pat::Photon> > photonColl;
   iEvent.getByToken(photonCollectionToken_, photonColl);

   // Get jet collection
   edm::Handle<pat::JetCollection> jetColl;
   iEvent.getByToken(jetCollectionToken_, jetColl);
   edm::Handle<reco::GenJetCollection> genJetColl;
   iEvent.getByToken(genJetCollectionToken_, genJetColl);
   edm::Handle<pat::MuonCollection> muonColl;
   iEvent.getByToken(muonCollectionToken_, muonColl);
   edm::Handle<edm::View<pat::Electron> > electronColl;
   iEvent.getByToken(electronCollectionToken_, electronColl);
   edm::Handle<pat::METCollection> metColl;
   iEvent.getByToken(metCollectionToken_, metColl);

   // Get PV
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   if (vertices->empty()) return; // skip the event if no PV found
   //const reco::Vertex &pv = vertices->front();
   nPV_    = vertices->size();

   reco::VertexCollection::const_iterator firstGoodVertex = vertices->end();
   nGoodVertices_=0;
   for (reco::VertexCollection::const_iterator vtx = vertices->begin();
        vtx != vertices->end(); ++vtx) {
      // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
      // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
      if (  /*!vtx->isFake() &&*/
         !(vtx->chi2()==0 && vtx->ndof()==0)
         &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
         && fabs(vtx->position().Z())<=24.0)
      {
         nGoodVertices_++;
         // first one?
         if (nGoodVertices_==1) firstGoodVertex = vtx;
      }
   }

   if (nGoodVertices_==0) return; // skip event if there are no good PVs
   hCutFlow_->Fill("nGoodVertices",1);

   // Get rho
   edm::Handle< double > rhoH;
   iEvent.getByToken(rhoToken_,rhoH);
   rho_ = *rhoH;

   // Get generator level info
   // Pruned particles are the one containing "important" stuff
   edm::Handle<edm::View<reco::GenParticle> > prunedGenParticles;
   if (!isRealData_){
      iEvent.getByToken(prunedGenToken_,prunedGenParticles);
   }

   // Jets
   vJets_.clear();
   tree::Jet trJet;
   for (const pat::Jet& jet : *jetColl){
      trJet.p.SetPtEtaPhi(jet.pt(),jet.eta(),jet.phi());
      trJet.bDiscriminator=jet.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
      trJet.someTestFloat=jet.chargedEmEnergyFraction();
      trJet.isLoose=isLooseJet(jet);
      vJets_.push_back(trJet);
   } // jet loop

   vGenJets_.clear();
   tree::Particle trGJet;
   for (const reco::GenJet& jet: *genJetColl){
      trGJet.p.SetPtEtaPhi(jet.pt(),jet.eta(),jet.phi());
      vGenJets_.push_back(trGJet);
   }

   if (vJets_.empty()) return;
   hCutFlow_->Fill("jets",1);

   double const HT=computeHT(vJets_);
   if (HT<dHT_cut_) return;
   hCutFlow_->Fill("HT",1);

   edm::Handle<edm::ValueMap<bool> > loose_id_dec;
   edm::Handle<edm::ValueMap<bool> > medium_id_dec;
   edm::Handle<edm::ValueMap<bool> > tight_id_dec;
   edm::Handle<edm::ValueMap<float>> mva_value;
   iEvent.getByToken(photonLooseIdMapToken_  ,loose_id_dec);
   iEvent.getByToken(photonMediumIdMapToken_ ,medium_id_dec);
   iEvent.getByToken(photonTightIdMapToken_  ,tight_id_dec);
   iEvent.getByToken(photonMvaValuesMapToken_,mva_value);
   // photon loop
   vPhotons_.clear();
   tree::Photon trPho;
   for(edm::View<pat::Photon>::const_iterator pho = photonColl->begin(); pho != photonColl->end(); pho++){
      // Kinematics
      if( pho->pt() < 15 )
         continue;

      trPho.p.SetPtEtaPhi(pho->pt(),pho->superCluster()->eta(),pho->superCluster()->phi());

      const edm::Ptr<pat::Photon> phoPtr( photonColl, pho - photonColl->begin() );

      trPho.sigmaIetaIeta=pho->full5x5_sigmaIetaIeta(); // from reco::Photon
      trPho.hOverE=pho->hadTowOverEm() ;
      trPho.hasPixelSeed=(Int_t)pho->hasPixelSeed() ;
      trPho.passElectronVeto= pho->passElectronVeto() ;
      trPho.r9  = pho->r9();

      trPho.isoChargedHadronsEA=pho->chargedHadronIso();
      trPho.isoNeutralHadronsEA=pho->neutralHadronIso();
      trPho.isoPhotonsEA       =pho->photonIso();

      // Compute isolation with effective area correction for PU
      // Find eta bin first. If eta>2.5, the last eta bin is used.
      int etaBin = 0;
      while ( etaBin < EffectiveAreas::nEtaBins-1
              && abs( pho->superCluster()->eta() ) > EffectiveAreas::etaBinLimits[etaBin+1] ){
         ++etaBin;
      };
      trPho.isoPhotonsEA        = std::max(float(0.0), trPho.isoPhotonsEA
                                           - rho_ * EffectiveAreas::areaPhotons[etaBin] );
      trPho.isoNeutralHadronsEA = std::max(float(0.0), trPho.isoNeutralHadronsEA
                                           - rho_ * EffectiveAreas::areaNeutralHadrons[etaBin] );
      trPho.isoChargedHadronsEA = std::max(float(0.0), trPho.isoChargedHadronsEA
                                           - rho_ * EffectiveAreas::areaChargedHadrons[etaBin] );

      trPho.mvaValue=(*mva_value)[phoPtr];

      // MC match
      if (!isRealData_){
         trPho.isTrue=matchToTruth(*pho, prunedGenParticles);
         trPho.isTrueAlternative=matchToTruthAlternative(*pho, prunedGenParticles);
      }else{
         trPho.isTrue=           UNMATCHED;
         trPho.isTrueAlternative=UNMATCHED;
      }

      // check photon working points
      trPho.isLoose = (*loose_id_dec) [phoPtr];
      trPho.isMedium= (*medium_id_dec)[phoPtr];
      trPho.isTight = (*tight_id_dec) [phoPtr];

      // write the photon to collection
      vPhotons_.push_back(trPho);
   } // photon loop

   if (vPhotons_.empty() || vPhotons_.at(0).p.Pt()<dPhoton_pT_cut_) return;
   hCutFlow_->Fill("photons",1);

   // Muons
   vMuons_.clear();
   tree::Muon trMuon;
   for (const pat::Muon &mu : *muonColl) {
      if (!mu.isLooseMuon()) continue;
      trMuon.p.SetPtEtaPhi(mu.pt(),mu.eta(),mu.phi());
      trMuon.isTight=mu.isTightMuon(*firstGoodVertex);
      // trMuon.someTestFloat=mu.isLooseMuon();
      vMuons_.push_back(trMuon);
   } // muon loop

   // Electrons
   // Get the electron ID data from the event stream
   edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
   edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
   edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
   edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
   iEvent.getByToken(electronVetoIdMapToken_  ,veto_id_decisions);
   iEvent.getByToken(electronLooseIdMapToken_ ,loose_id_decisions);
   iEvent.getByToken(electronMediumIdMapToken_,medium_id_decisions);
   iEvent.getByToken(electronTightIdMapToken_ ,tight_id_decisions);

   vElectrons_.clear();
   tree::Electron trEl;
   for(edm::View<pat::Electron>::const_iterator el = electronColl->begin();el != electronColl->end(); el++){
      const edm::Ptr<pat::Electron> elPtr(electronColl, el - electronColl->begin() );
      if (!(*veto_id_decisions)[elPtr]) continue; // take only 'veto' electrons
      trEl.isLoose =(*loose_id_decisions) [elPtr];
      trEl.isMedium=(*medium_id_decisions)[elPtr];
      trEl.isTight =(*tight_id_decisions) [elPtr];
      trEl.p.SetPtEtaPhi(el->pt(),el->superCluster()->eta(),el->superCluster()->phi());
      vElectrons_.push_back(trEl);
   }

   // MET
   const pat::MET &met = metColl->front();
   pat::MET::LorentzVector metRaw=met.shiftedP4(pat::MET::NoShift, pat::MET::Raw);
   double metPt=met.pt();
   met_.p.SetPtEtaPhi(metPt,met.eta(),met.phi());
   met_.p_raw.SetPtEtaPhi(metRaw.pt(),metRaw.eta(),metRaw.phi());

   // jet resolution shift is set to 0 for 74X
   met_.uncertainty=0;
   // loop over all up-shifts save for last one (=NoShift)
   for (uint iShift=0; iShift<(pat::MET::METUncertaintySize-1); iShift+=2){
      // up and down shifts
      const double u=fabs(met.shiftedPt(pat::MET::METUncertainty(iShift))  -metPt);
      const double d=fabs(met.shiftedPt(pat::MET::METUncertainty(iShift+1))-metPt);
      // average
      const double a=.5*(u+d);
      // add deviations in quadrature
      met_.uncertainty+=a*a;
   }
   met_.uncertainty=TMath::Sqrt(met_.uncertainty);

   // Generated Particles
   vGenPhotons_  .clear();
   vGenElectrons_.clear();
   tree::Particle trP;
   if (!isRealData_){
      for (const reco::GenParticle &genP: *prunedGenParticles){
         if (genP.status() != 1) continue; // only final state particles
         if (genP.pt() < 30)     continue;
         if (abs(genP.pdgId()) == 11){ // electron
            trP.p.SetPtEtaPhi(genP.pt(),genP.eta(),genP.phi());
            vGenElectrons_.push_back(trP);
         } else if (genP.pdgId() == 22){ // photon
            trP.p.SetPtEtaPhi(genP.pt(),genP.eta(),genP.phi());
            vGenPhotons_.push_back(trP);
         }
      }
   }

   // PileUp weights
   if (!isRealData_){
      edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
      iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      float Tnpv = -1;
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
         int BX = PVI->getBunchCrossing();
         if(BX == 0) {
            Tnpv = PVI->getTrueNumInteractions();
            continue;
         }
      }
      true_nPV_=Tnpv;
      pu_weight_=hPU_.GetBinContent(hPU_.FindBin(Tnpv));
   }else{ // real data
      true_nPV_=-1;
      pu_weight_=1.;
   }

   // generator weights
   mc_weight_=0.;
   if (!isRealData_){
      edm::Handle<GenEventInfoProduct> GenEventInfoHandle;
      iEvent.getByLabel("generator", GenEventInfoHandle);
      mc_weight_=GenEventInfoHandle->weight();
   }
   
   hCutFlow_->Fill("final",1);
   // store event identity
   evtNo_=iEvent.id().event();
   runNo_=iEvent.run();
   lumNo_=iEvent.luminosityBlock();
   // write the event
   eventTree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void
TreeWriter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
TreeWriter::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
/*
  void
  TreeWriter::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void
  TreeWriter::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
  void
  TreeWriter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void
  TreeWriter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TreeWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
   //The following says we do not know what parameters are allowed so do no validation
   // Please change this to state exactly what you do use, even if it is no parameters
   edm::ParameterSetDescription desc;
   desc.setUnknown();
   descriptions.addDefault(desc);
}

int TreeWriter::matchToTruth(const pat::Photon &pho,
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

void TreeWriter::findFirstNonPhotonMother(const reco::Candidate *particle,
                                                           int &ancestorPID, int &ancestorStatus){

   if( particle == 0 ){
      printf("TreeWriter: ERROR! null candidate pointer, this should never happen\n");
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

int TreeWriter::matchToTruthAlternative(const pat::Photon &pho,
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
DEFINE_FWK_MODULE(TreeWriter);
