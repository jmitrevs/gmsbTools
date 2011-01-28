#ifndef GMSBTOOLS_GMSBPREPARATIONTOOL_H  
#define GMSBTOOLS_GMSBPREPARATIONTOOL_H 

/*****************************************************************************
Name    : gmsbPreparationTool.h
Package : offline/PhysicsAnalysis/AnalysisCommon/gmsbTools
Author  : Ketevi A. Assamagan
Created : November 2007
Purpose : User tools for analyis preparation on ESD/AOD/DPD in Athena
          - selections
          - write out contianer of selected particles to StoreGate
*****************************************************************************/

#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ToolHandle.h"

#include "gmsbTools/gmsbSelectionTool.h"

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

class MsgStream;

/** Interface ID */  
static const InterfaceID IID_gmsbPreparationTool("gmsbPreparationTool", 1, 0);

class gmsbPreparationTool : public AthAlgTool {

public:

  /** Standard Constructor */
  gmsbPreparationTool(const std::string& type, const std::string& name,
	                      const IInterface* parent);

  /** AlgTool and IAlgTool interface methods */
  static const InterfaceID& interfaceID() { return IID_gmsbPreparationTool; };

  /** Overriding initialize, finalize, and execute */
  virtual StatusCode initialize();
  virtual StatusCode execute();
  virtual StatusCode finalize();

  /** access to containers after preparation */
  const ElectronContainer             * selectedElectrons();
  const PhotonContainer               * selectedPhotons();  
  const Analysis::MuonContainer       * selectedMuons();
  const Analysis::TauJetContainer     * selectedTauJets();
  const JetCollection                 * selectedJets();
  const Rec::TrackParticleContainer   * selectedTrackParticles();
  const CaloClusterContainer          * selectedCaloClusters();

  /** summary of pre-selections and overlap removal - will be called at the end of the job
      in the finalize of this tool - the first number is reconstrued and the second is the pre-selected */
  void summarize();
  const std::pair<unsigned int, unsigned int>& electronSummary() const;
  const std::pair<unsigned int, unsigned int>& photonSummary() const;
  const std::pair<unsigned int, unsigned int>& muonSummary() const;
  const std::pair<unsigned int, unsigned int>& tauJetSummary() const;
  const std::pair<unsigned int, unsigned int>& jetSummary() const;
  const std::pair<unsigned int, unsigned int>& trackParticleSummary() const;
  const std::pair<unsigned int, unsigned int>& caloClusterSummary() const;

protected:

   /** Standard destructor */
   virtual ~gmsbPreparationTool();

private:

  /** container preparation */
  StatusCode electronPreparation( std::string key, int runNum );
  StatusCode photonPreparation( std::string key, int runNum );
  StatusCode muonPreparation( std::string key );
  StatusCode tauJetPreparation( std::string key );
  StatusCode jetPreparation( std::string key );
  StatusCode trackParticlePreparation( std::string key );
  StatusCode caloClusterPreparation( std::string key );

  /** for debugging purposes - called if MSG_Level = DEBUG */
  void print();

private:

  /** a handle on selection */
  ToolHandle <gmsbSelectionTool> m_userSelectionTool;

  /** should contain the StoreGate keys to be passed in job options */ 
  std::vector<std::string> m_inputContainerKeys;
 
  /** number of various particles <before selection, after selection> 
      used in the summarize() method print summary information */
  std::pair<unsigned int, unsigned int> m_numElectrons;
  std::pair<unsigned int, unsigned int> m_numPhotons;
  std::pair<unsigned int, unsigned int> m_numMuons;
  std::pair<unsigned int, unsigned int> m_numTauJets;
  std::pair<unsigned int, unsigned int> m_numJets;
  std::pair<unsigned int, unsigned int> m_numTrackParticles;
  std::pair<unsigned int, unsigned int> m_numCaloClusters; 

  /** output collection prefix and keys 
      the output collection key are built form the inputCollectionKeys with the prefix appended
      the use can set the prefix in the job options */
  std::vector<std::string> m_outputContainerKeys;
  std::string m_outputElectronKey;
  std::string m_outputPhotonKey;
  std::string m_outputMuonKey;
  std::string m_outputTauJetKey;
  std::string m_outputJetKey;
  std::string m_outputTrackParticleKey;
  std::string m_outputCaloClusterKey;

  /** is ATLFAST data */
  bool m_isAtlfast;

  /** doing debugging */
  bool m_debug;

  /** on first event */
  bool m_first;

  /** run number to use for OQ: -1 means use event run number */
  int m_OQRunNum;

};

inline const std::pair<unsigned int, unsigned int>& gmsbPreparationTool::electronSummary() const
{
  return m_numElectrons;
}

inline const std::pair<unsigned int, unsigned int>& gmsbPreparationTool::photonSummary() const
{
  return m_numPhotons;
}

inline const std::pair<unsigned int, unsigned int>& gmsbPreparationTool::muonSummary() const
{
  return m_numMuons;
}

inline const std::pair<unsigned int, unsigned int>& gmsbPreparationTool::tauJetSummary() const
{
  return m_numTauJets;
}

inline const std::pair<unsigned int, unsigned int>& gmsbPreparationTool::jetSummary() const
{
  return m_numJets;
}

inline const std::pair<unsigned int, unsigned int>& gmsbPreparationTool::trackParticleSummary() const
{
  return m_numTrackParticles;
}

inline const std::pair<unsigned int, unsigned int>& gmsbPreparationTool::caloClusterSummary() const
{
  return m_numCaloClusters;
}

#endif // GMSBTOOLS_GMSBPREPARATIONTOOL_H 




