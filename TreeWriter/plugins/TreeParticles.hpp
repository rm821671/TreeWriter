#ifndef TREEPARTICLES_H
#define TREEPARTICLES_H

#include <TLorentzVector.h>
#include <TVector3.h>

namespace tree 
{
   struct Particle
   {
      TVector3 p;
      float    someTestFloat=0.; // only for quick peeking
   };

   struct Photon : public Particle
   {
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

   struct Jet : public Particle
   {
      bool isLoose;
      float bDiscriminator;
   };
   // struct Jet : public Particle{
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

   struct Muon: public Particle
   {
      bool isTight;
   };

   struct Electron: public Particle
   {
      bool isLoose;
      bool isMedium;
      bool isTight;
   };

   struct MET : public Particle
   {
      TVector3 p_raw;
      Float_t  uncertainty;
   };

   bool EtGreater(const tree::Particle, const tree::Particle);

} // end namespace definition

#endif /* TREEPARTICLES_H */
