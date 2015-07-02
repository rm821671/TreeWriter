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
      Float_t sigmaIetaIeta; // full 5x5
      Float_t hOverE;
      Int_t hasPixelSeed;
      Int_t passElectronVeto;
      Float_t r9;

      Float_t isoChargedHadronsEA;
      Float_t isoNeutralHadronsEA;
      Float_t isoPhotonsEA;
  
      Int_t isTrue;
      Int_t isTrueAlternative;

      // IDs
      Bool_t  isLoose;
      Bool_t  isMedium;
      Bool_t  isTight;
      Float_t mvaValue;
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

   inline bool EtGreater(const tree::Particle p1, const tree::Particle p2) {
      return p1.p.Pt() > p2.p.Pt();
   }

} // end namespace definition
#endif /* TREEPARTICLES_H */
