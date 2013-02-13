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

#include <string>
#include <map>
#include <vector>

#include "gmsbTools/FourMomHelpers.h"

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
  bool overlap(float etaA, float phiA, float etaB, float phiB, bool jet) const ;

protected:

   /** Standard destructor */
   virtual ~gmsbOverlapCheckingTool();

private:


  /** deltaR overlap */
  double m_deltaR;
  double m_deltaRWithJets;

};

/** check for oeverlap in deltaR and return as well the deltaR value */
inline bool gmsbOverlapCheckingTool::overlap(float etaA, float phiA, float etaB, float phiB, bool jet) const 
{
  return FourMomHelpers::isInDeltaR(etaA, phiA, etaB, phiB, jet ? m_deltaRWithJets : m_deltaR);
}

#endif // GMSBTOOLS_GMSBOVERLAPCHECKINGTOOL_H 




