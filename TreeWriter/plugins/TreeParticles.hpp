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
      Float_t hOverE;
      Int_t hasPixelSeed;
      Int_t passElectronVeto;

      // Extra variables for MVA ID (excluding already mentioned above)
      Float_t scRawEnergy;
      Float_t r9;
      Float_t sigma_eta;
      Float_t sigma_phi;
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
