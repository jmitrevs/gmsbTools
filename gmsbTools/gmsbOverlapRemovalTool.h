#ifndef GMSBTOOLS_GMSBOVERLAPREMOVALTOOL_H  
#define GMSBTOOLS_GMSBOVERLAPREMOVALTOOL_H 

/*****************************************************************************
Name    : gmsbOverlapRemovalTool.h
Package : offline/PhysicsAnalysis/AnalysisCommon/gmsbTools
Author  : Ketevi A. Assamagan
Created : November 2007
Purpose : User tools for analyis overlap removal on ESD/AOD/DPD in Athena
          - do overlap removal given a set of containers
          - Return lists of non-overlapping Particles, leptons, etc
	  - Call overlap checking tools down do cell and hit level 
*****************************************************************************/

#include "AthenaBaseComps/AthAlgTool.h"

#include "GaudiKernel/ToolHandle.h"

#include "gmsbTools/gmsbSelectionTool.h"
#include "gmsbTools/gmsbOverlapCheckingTool.h"


#include <string>
#include <map>
#include <vector>

/** Interface ID */  
static const InterfaceID IID_gmsbOverlapRemovalTool("gmsbOverlapRemovalTool", 1, 0);

class ElectronD3PDObject;
class MuonD3PDObject;
class JetD3PDObject;
class PhotonD3PDObject;

class genericParticle {
public:
  enum particleType {photon, electron, muon, jet}
  float pt() const {return m_pt;};
  float eta() const {return m_eta;};
  float phi() const {return m_phi;};
  particleType type() const {return m_isJet;} 
  float cl_eta() const {return m_cl_eta;};
  float cl_phi() const {return m_cl_phi;};
  genericParticle(float pt, float eta, float phi, particleType type) :
    m_pt(pt), m_eta(eta), m_phi(phi), m_type(type), 
    m_etaClus(-999), m_phiClus(-999) {};
  genericParticle(float pt, float eta, float phi, particleType type,
		  float etaClus, float phiClus) :
    m_pt(pt), m_eta(eta), m_phi(phi), m_type(type),
    m_cl_eta(etaClus), m_cl_phi(phiClus) {};
private:
  float m_pt;
  float m_eta;
  float m_phi;
  praticleType m_type;
  float m_cl_eta; // cluster values for electrons/photons
  float m_cl_phi; // cluster values for electrons/photons
}

class gmsbOverlapRemovalTool : public AthAlgTool {

public:

  /** Standard Constructor */
  gmsbOverlapRemovalTool(const std::string& type, const std::string& name,
	                         const IInterface* parent);

  /** AlgTool and IAlgTool interface methods */
  static const InterfaceID& interfaceID() { return IID_gmsbOverlapRemovalTool; };

  /** Overriding initialize, finalize, and execute */
  virtual StatusCode initialize();
  virtual StatusCode execute();
  virtual StatusCode finalize();

  /** access to containers after preparation */
  /// NOTE: These are factories:  The user needs to delete the created object
  const PhotonD3PDObject               * finalStatePhotons();  // including converted photons
  const ElectronD3PDObject             * finalStateElectrons();
  const MuonD3PDObject       * finalStateMuons();
  const JetD3PDObject          * finalStateJets();
  // const JetD3PDObject          * finalStateBJets();
  // const JetD3PDObject          * finalStateLightJets();

  /** summary of pre-selections and overlap removal - will be called at the end of the job
      in the finalize of this tool - the first number is reconstrued and the second is the pre-selected */
  void summarize();
  const std::pair<unsigned int, unsigned int>& electronSummary() const;
  const std::pair<unsigned int, unsigned int>& photonSummary() const;
  const std::pair<unsigned int, unsigned int>& muonSummary() const;
  const std::pair<unsigned int, unsigned int>& jetSummary() const;
  //const std::pair<unsigned int, unsigned int>& bJetSummary() const;
  //const std::pair<unsigned int, unsigned int>& lightJetSummary() const;

protected:

   /** Standard destructor */
   virtual ~gmsbOverlapRemovalTool();

private:

  /** container preparation */
  StatusCode electronPreparation( std::string key );
  StatusCode photonPreparation( std::string key );
  StatusCode muonPreparation( std::string key );
  StatusCode jetPreparation( std::string key );

  /** for debugging purposes - called if MSG_Level = DEBUG */
  void print();


private:

  /** a handle on selection  and on overlap checking */
  ToolHandle <gmsbSelectionTool> m_userSelectionTool;
  ToolHandle <gmsbOverlapCheckingTool> m_userOverlapCheckingTool;

  /** should contain the StoreGate keys to be passed in job options */ 
  std::vector<std::string> m_inputContainerKeys;
 
  /** number of various particles <before selection, after selection> 
      used in the summarize() method print summary information */
  std::pair<unsigned int, unsigned int> m_numElectrons;
  std::pair<unsigned int, unsigned int> m_numPhotons;
  std::pair<unsigned int, unsigned int> m_numMuons;
  std::pair<unsigned int, unsigned int> m_numJets;
  // std::pair<unsigned int, unsigned int> m_numBJets;
  // std::pair<unsigned int, unsigned int> m_numLightJets;

  /** output collection prefix and keys 
      the output collection key are built form the inputCollectionKeys with the prefix appended
      the use can set the prefix in the job options */
  //std::string m_outputObjectKey;
  //std::string m_outputLeptonKey;
  std::string m_outputElectronKey;
  std::string m_outputPhotonKey;
  std::string m_outputMuonKey;
  std::string m_outputJetKey;
  //std::string m_outputBJetKey;
  //std::string m_outputLightJetKey;

  /** is ATLFAST data */
  bool m_isAtlfast;

  /** remove overlay in same container */
  bool m_removeOverlapInSameContainer;

  /** doing debugging */
  bool m_debug;

  std::vector<genericParticle> m_allParticles;

};

inline const std::pair<unsigned int, unsigned int>& gmsbOverlapRemovalTool::electronSummary() const
{
  return m_numElectrons;
}

inline const std::pair<unsigned int, unsigned int>& gmsbOverlapRemovalTool::photonSummary() const
{
  return m_numPhotons;
}

inline const std::pair<unsigned int, unsigned int>& gmsbOverlapRemovalTool::muonSummary() const
{
  return m_numMuons;
}

inline const std::pair<unsigned int, unsigned int>& gmsbOverlapRemovalTool::jetSummary() const
{
  return m_numJets;
}

// inline const std::pair<unsigned int, unsigned int>& gmsbOverlapRemovalTool::bJetSummary() const
// {
//   return m_numBJets;
// }

// inline const std::pair<unsigned int, unsigned int>& gmsbOverlapRemovalTool::lightJetSummary() const
// {
//   return m_numLightJets;
// }


#endif // GMSBTOOLS_GMSBOVERLAPREMOVALTOOL_H 




