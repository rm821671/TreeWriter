#include <TLorentzVector.h>
#include <TVector3.h>

#ifndef TREEPARTICLES_H
#define TREEPARTICLES_H


namespace tree {

   class Particle {
   public:
      TVector3 p;
      float    someTestFloat=0.;
   };

   class Photon : public Particle {
   public:
      // Variables for cut based ID
      Float_t full5x5_sigmaIetaIeta;
      Float_t hOverE;
      Int_t hasPixelSeed;
      Int_t passElectronVeto;

      Float_t isoChargedHadrons;
      Float_t isoNeutralHadrons;
      Float_t isoPhotons;

      Float_t isoChargedHadronsWithEA;
      Float_t isoNeutralHadronsWithEA;
      Float_t isoPhotonsWithEA;

      // Extra variables for MVA ID (excluding already mentioned above)
      Float_t scRawEnergy;
      Float_t isoWorstChargedHadrons;
      Float_t r9;
      Float_t full5x5_sigmaIetaIphi;
      Float_t full5x5_e1x3;
      Float_t full5x5_e2x2;
      Float_t full5x5_e2x5Max;
      Float_t full5x5_e5x5;
      Float_t sigma_eta;
      Float_t sigma_phi;
      Float_t esEffSigmaRR;
      Float_t esEnergy;

      Float_t mvaValue;
  
      Int_t isTrue;
      Int_t isTrueAlternative;

      Bool_t  isLoose;
      Bool_t  isMedium;
      Bool_t  isTight;

   };

   class Jet : public Particle{
   public:
      float bDiscriminator;
   };
   // class Jet : public Particle{
   // public:
   //    float bCSV;
   //    float chargedHadronEnergy,
   // 	 neutralHadronEnergy,
   // 	 photonEnergy,
   // 	 electronEnergy,
   // 	 muonEnergy,
   // 	 HFHadronEnergy,
   // 	 HFEMEnergy,
   // 	 chargedEmEnergy,
   // 	 chargedMuEnergy,
   // 	 neutralEmEnergy;
   // };

   bool EtGreater(const tree::Particle, const tree::Particle);

} // end namespace definition

#endif /* TREEPARTICLES_H */
