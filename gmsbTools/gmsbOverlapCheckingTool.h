#ifndef GMSBTOOLS_GMSBOVERLAPCHECKINGTOOL_H  
#define GMSBTOOLS_GMSBOVERLAPCHECKINGTOOL_H 

/*****************************************************************************
Name    : gmsbOverlapCheckingTool.h
Package : gmsbTools
Author  : Ketevi A. Assamagan (adopted by Jovan Mitrevski)
Created : November 2007
Purpose : User tools for checking overlaps at deltaR, TrackParticle/Cluster and Hit/Cell levels
*****************************************************************************/

#include "AthenaBaseComps/AthAlgTool.h"

#include "GaudiKernel/ToolHandle.h"

#include "VxVertex/VxContainer.h"
#include "Particle/TrackParticleContainer.h"
#include "CaloEvent/CaloClusterContainer.h"
#include "TrkSegment/SegmentCollection.h"

#include "muonEvent/MuonContainer.h"
#include "egammaEvent/ElectronContainer.h"
#include "egammaEvent/PhotonContainer.h"
#include "tauEvent/TauJetContainer.h"
#include "JetEvent/JetCollection.h"
#include "MissingETEvent/MissingET.h"

#include "NavFourMom/IParticleContainer.h"
#include "NavFourMom/INavigable4MomentumCollection.h"

#include "FourMomUtils/P4Helpers.h"

#include <string>
#include <map>
#include <vector>

/** Interface ID */  
static const InterfaceID IID_gmsbOverlapCheckingTool("gmsbOverlapCheckingTool", 1, 0);

class gmsbOverlapCheckingTool : public AthAlgTool {

public:

  /** Standard Constructor */
  gmsbOverlapCheckingTool(const std::string& type, const std::string& name,
	                          const IInterface* parent);

  /** AlgTool and IAlgTool interface methods */
  static const InterfaceID& interfaceID() { return IID_gmsbOverlapCheckingTool; };

  /** Overriding initialize, finalize, and execute */
  virtual StatusCode initialize();
  virtual StatusCode finalize();

  /** overlaps */
  bool overlap(const I4Momentum *object1, const I4Momentum *object2) const ;

protected:

   /** Standard destructor */
   virtual ~gmsbOverlapCheckingTool();

private:

  /** deltaR overlap */
  double m_deltaR;
  double m_deltaRWithJets;

};

/** check for oeverlap in deltaR and return as well the deltaR value */
inline bool gmsbOverlapCheckingTool::overlap(const I4Momentum *object1,
                                                     const I4Momentum *object2) const 
{

  const Jet * jet1 = dynamic_cast<const Jet*> (object1);
  const Jet * jet2 = dynamic_cast<const Jet*> (object2);
  if ( jet1 || jet2 ) {
    return P4Helpers::isInDeltaR( *object1, *object2, m_deltaRWithJets);
  } else {
    return P4Helpers::isInDeltaR( *object1, *object2, m_deltaR);
  }
}

#endif // GMSBTOOLS_GMSBOVERLAPCHECKINGTOOL_H 




