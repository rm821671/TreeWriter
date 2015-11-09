#ifndef TREEPARTICLES_H
#define TREEPARTICLES_H

#include <TLorentzVector.h>
#include <TVector3.h>

namespace tree
{
   struct Particle
   {
      TVector3 p;
   };

   struct GenParticle: public Particle
   {
      Int_t pdgId=0;
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
      Float_t isoWorstChargedHadrons;

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
      bool hasPhotonMatch;
      bool hasElectronMatch;
      float bDiscriminator;
   };

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
   
   struct PFCandidate : public Particle
   {
      Int_t pdgId=0;      
      Int_t charge;
      
      // 0: no match
      // 1: electron/positron id-charge match
      // 2: pi+/pi- id-charge match
      //Int_t IdChargeMatch;
      
      // primary vertex information
      //
      // enum PVAssociationQuality {
      // 	NotReconstructedPrimary=0,
      // 	OtherDeltaZ=1,
      // 	CompatibilityBTag=4,
      // 	CompatibilityDz=5,
      // 	UsedInFitLoose=6,
      // 	UsedInFitTight=7	
      // };
      Int_t pvAssociationQuality;
      
      
   };

   inline bool EtGreater(const tree::Particle p1, const tree::Particle p2) {
      return p1.p.Pt() > p2.p.Pt();
   }

} // end namespace definition
#endif /* TREEPARTICLES_H */
