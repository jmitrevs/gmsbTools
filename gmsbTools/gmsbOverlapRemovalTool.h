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

#include <string>
#include <map>
#include <vector>

/** Interface ID */  
static const InterfaceID IID_gmsbOverlapRemovalTool("gmsbOverlapRemovalTool", 1, 0);

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
  const INavigable4MomentumCollection * finalStateObjects();
  const PhotonContainer               * finalStatePhotons();  // including converted photons
  const ElectronContainer             * finalStateElectrons();
  const Analysis::MuonContainer       * finalStateMuons();
  const INavigable4MomentumCollection * finalStateLeptons();  // Electrons or Muons
  const Analysis::TauJetContainer     * finalStateTauJets();
  const JetCollection          * finalStateJets();
  const JetCollection          * finalStateBJets();
  const JetCollection          * finalStateLightJets();
  const Rec::TrackParticleContainer   * finalStateTrackParticles();
  const CaloClusterContainer          * finalStateCaloClusters();

  /** summary of pre-selections and overlap removal - will be called at the end of the job
      in the finalize of this tool - the first number is reconstrued and the second is the pre-selected */
  void summarize();
  const std::pair<unsigned int, unsigned int>& electronSummary() const;
  const std::pair<unsigned int, unsigned int>& photonSummary() const;
  const std::pair<unsigned int, unsigned int>& muonSummary() const;
  const std::pair<unsigned int, unsigned int>& tauJetSummary() const;
  const std::pair<unsigned int, unsigned int>& jetSummary() const;
  const std::pair<unsigned int, unsigned int>& bJetSummary() const;
  const std::pair<unsigned int, unsigned int>& lightJetSummary() const;
  const std::pair<unsigned int, unsigned int>& trackParticleSummary() const;
  const std::pair<unsigned int, unsigned int>& caloClusterSummary() const;

  /** check if execute() is already called for this tool in this job for this event */
  bool isExecuted();

protected:

   /** Standard destructor */
   virtual ~gmsbOverlapRemovalTool();

private:

  /** container preparation */
  StatusCode prepareContainers();
  StatusCode electronPreparation( std::string key );
  StatusCode photonPreparation( std::string key );
  StatusCode muonPreparation( std::string key );
  StatusCode tauJetPreparation( std::string key );
  StatusCode jetPreparation( std::string key );
  StatusCode trackParticlePreparation( std::string key );
  StatusCode caloClusterPreparation( std::string key );
  StatusCode lockContainers();

  /** for debugging purposes - called if MSG_Level = DEBUG */
  void print();

private:

  INavigable4MomentumCollection * allParticles();
  INavigable4MomentumCollection * allLeptons();
  PhotonContainer               * allPhotons();  // including converted photons
  ElectronContainer             * allElectrons();
  Analysis::MuonContainer       * allMuons();
  Analysis::TauJetContainer     * allTauJets();
  JetCollection          * allJets();
  JetCollection          * allBJets();
  JetCollection          * allLightJets();
  Rec::TrackParticleContainer   * allTrackParticles();
  CaloClusterContainer          * allCaloClusters();

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
  std::pair<unsigned int, unsigned int> m_numTauJets;
  std::pair<unsigned int, unsigned int> m_numJets;
  std::pair<unsigned int, unsigned int> m_numBJets;
  std::pair<unsigned int, unsigned int> m_numLightJets;
  std::pair<unsigned int, unsigned int> m_numTrackParticles;
  std::pair<unsigned int, unsigned int> m_numCaloClusters; 

  /** output collection prefix and keys 
      the output collection key are built form the inputCollectionKeys with the prefix appended
      the use can set the prefix in the job options */
  std::string m_outputObjectKey;
  std::string m_outputLeptonKey;
  std::string m_outputElectronKey;
  std::string m_outputPhotonKey;
  std::string m_outputMuonKey;
  std::string m_outputTauJetKey;
  std::string m_outputJetKey;
  std::string m_outputBJetKey;
  std::string m_outputLightJetKey;
  std::string m_outputTrackParticleKey;
  std::string m_outputCaloClusterKey;

  /** is ATLFAST data */
  bool m_isAtlfast;

  /** remove overlay in same container */
  bool m_removeOverlapInSameContainer;

  /** doing debugging */
  bool m_debug;

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

inline const std::pair<unsigned int, unsigned int>& gmsbOverlapRemovalTool::tauJetSummary() const
{
  return m_numTauJets;
}

inline const std::pair<unsigned int, unsigned int>& gmsbOverlapRemovalTool::jetSummary() const
{
  return m_numJets;
}

inline const std::pair<unsigned int, unsigned int>& gmsbOverlapRemovalTool::bJetSummary() const
{
  return m_numBJets;
}

inline const std::pair<unsigned int, unsigned int>& gmsbOverlapRemovalTool::lightJetSummary() const
{
  return m_numLightJets;
}

inline const std::pair<unsigned int, unsigned int>& gmsbOverlapRemovalTool::trackParticleSummary() const
{
  return m_numTrackParticles;
}

inline const std::pair<unsigned int, unsigned int>& gmsbOverlapRemovalTool::caloClusterSummary() const
{
  return m_numCaloClusters;
}

#endif // GMSBTOOLS_GMSBOVERLAPREMOVALTOOL_H 




