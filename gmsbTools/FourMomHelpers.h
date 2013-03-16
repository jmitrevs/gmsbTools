#ifndef GMSBTOOLS_FOURMOMHELPERS_H  
#define GMSBTOOLS_FOURMOMHELPERS_H 

#include <cmath>

namespace FourMomHelpers {

  inline float deltaPhi( float phiA, float phiB )
  {
    return  -remainderf( -phiA + phiB, 2*M_PI );
  }

  inline float deltaEta( float etaA, float etaB )
  {
    return (etaA - etaB);
  }

  inline bool isInDeltaR( float etaA, float phiA, float etaB, float phiB,
			  float dR )
  {
    using std::abs;
    using std::sqrt;
    const double dPhi = abs(deltaPhi(phiA, phiB));
    if ( dPhi > dR ) return false;                         // <==
    const double dEta = abs(deltaEta(etaA, etaB));
    if ( dEta > dR ) return false;                         // <==
    const double deltaR2 = dEta*dEta + dPhi*dPhi;
    if ( deltaR2 > dR*dR ) return false;                   // <==
    return true;
  }


  inline float deltaR( float etaA, float phiA, float etaB, float phiB )
  {
    return hypotf(deltaPhi(phiA, phiB), deltaEta(etaA, etaB));
  }			
}
#endif
