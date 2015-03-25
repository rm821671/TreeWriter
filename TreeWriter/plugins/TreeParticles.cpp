#include "TreeParticles.hpp"



bool tree::EtGreater(const tree::Particle p1, const tree::Particle p2) {
   return p1.p.Pt() > p2.p.Pt();
}

