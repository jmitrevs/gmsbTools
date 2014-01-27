#ifndef GMSBTOOLS_GMSBSELECTIONTOOL_H  
#define GMSBTOOLS_GMSBSELECTIONTOOL_H 

/*****************************************************************************
Name    : gmsbSelectionTool.h
Package : offline/PhysicsAnalysis/AnalysisCommon/gmsbTools
Author  : Ketevi A. Assamagan
Created : November 2007
Purpose : User tools for analyis preparation on ESD/AOD/DPD in Athena - selections
          - Take a list of input containers
          - Call the selections tools to see is containee passed selection
          - if passed selection put containee in a new container  
*****************************************************************************/

#include "AthenaBaseComps/AthAlgTool.h"

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

#include "egammaAnalysisUtils/checkOQ.h"
#include "egammaAnalysisUtils/EnergyRescaler.h"

#include "AthenaKernel/IUserDataSvc.h"
#include "GaudiKernel/ToolHandle.h"

#include <string>
#include <map>
#include <vector>

class IPAUcaloIsolationTool;


/** Interface ID */  
static const InterfaceID IID_gmsbSelectionTool("gmsbSelectionTool", 1, 0);

class IMCTruthClassifier;

class gmsbSelectionTool : public AthAlgTool {

public:

  /** Standard Constructor */
  gmsbSelectionTool(const std::string& type, const std::string& name,
	                    const IInterface* parent);

  /** AlgTool and IAlgTool interface methods */
  static const InterfaceID& interfaceID() { return IID_gmsbSelectionTool; };

  /** Overriding initialize, finalize, and execute */
  virtual StatusCode initialize();
  virtual StatusCode finalize();

  /** pre-selections */
  bool isSelected( const Analysis::Electron * electron, 
		   unsigned int runNum = 0,
		   unsigned int nPV = 0) const;
  bool isSelected( const Analysis::Photon * photon, 
		   unsigned int runNum = 0,
		   unsigned int nPV = 0) const;
  bool isSelected( const Analysis::Muon * muon ) const;
  bool isSelected( const Analysis::TauJet * tauJet ) const;
  bool isSelected( const Jet* jet ) const;
  bool isSelected( const Rec::TrackParticle * trackParticle ) const;
  bool isSelected( const CaloCluster* caloCluster ) const;

  bool isBJet( const Jet * jet ) const;

protected:

   /** Standard destructor */
   virtual ~gmsbSelectionTool();

private:

  /** this is Atlfast data */
  bool m_isAtlfast;

  bool m_simple; 		// if true, don't smear or decorate event


  /** MC */
  bool m_isMC;
  bool m_doTruth;		// require electron and photon to be true
  bool m_smearMC;
  bool m_MCHasConstantTerm;
  // int m_randomSeed; // use SUSY prescription
  int m_elScaleShift;
  int m_elSmearShift;
  int m_phoScaleShift;
  int m_phoSmearShift;
  double m_mcEtconeScale;
  bool m_useAltIsoCorrection;

  /** Electron selection */
  double m_electronPt;
  double m_electronEta;
  int    m_electronID;
  bool   m_doOldElectronIsolation;
  bool   m_doNewElectronIsolation;
  bool   m_doElectronTrackIsolation;
  bool   m_authorEgammaOnly;
  bool   m_doElectronEtaWindCut;/// apply (or not) eta cut in bad region window   
  bool   m_doEDElectronIsolation;
 
  double m_electronEtaWindMin;
  double m_electronEtaWindMax;
  double m_electronEtcone20ovEt;  // for old
  double m_electronEtcone20corrected; // for new
  double m_electronPtcone20ovEt;  // for track

  /** Photon selection */
  double m_photonPt;
  double m_photonEta;
  int    m_photonID;
  unsigned int m_photonIsEM;
  bool   m_doOldPhotonIsolation;
  bool   m_doNewPhotonIsolation;
  bool   m_doEDPhotonIsolation;
  bool   m_doPhotonEtaWindCut;/// apply (or not) eta cut in bad region window  
  double m_photonEtaWindMin;
  double m_photonEtaWindMax;
  double m_photonEtcone20ovEt;  // for old
  double m_photonEtcone20corrected;  // for new

  /** Muon selection */
  double m_muonPt;
  double m_muonEta;
  // double m_matchChi2Max;
  bool m_sel_combined;
  bool m_sel_seg_tag;
  bool m_do_iso_cut;
  bool m_do_flat_iso_cut;
  float m_flat_isolation_cut;
  float m_isolation_cut;
  //float m_ms_pt_limit;
  //float m_ms_p_diff_limit;
  std::string m_muonResSyst;  // Valid values: {"MSLOW", "MSUP", "IDLOW", "IDUP"} 

  /** TauJet selection */
  double m_tauJetPt;
  double m_tauJetEta;
  double m_tauJetLikelihood;
  double m_tauElTauLikelihood;

  /** Jet selection */
  double m_jetPt;
  double m_jetEta;
  double m_bJetLikelihood;
  double m_rejNegEJets;		// reject jets with neg E

  /** caloCluster selection */
  double m_caloClusterE;

  /** TrackParticle selection */
  double m_trackParticlePt;


  // user data
  ServiceHandle<IUserDataSvc> m_userdatasvc;

  ToolHandle<IMCTruthClassifier> m_MCTruthClassifier;

  /** @brief Tool handle for corrected isolation */
  ToolHandle<IPAUcaloIsolationTool> m_PAUcaloIsolationTool;

  // the OQ utility
  egammaOQ m_OQ;

  mutable eg2011::EnergyRescaler m_eRescale;

};

#endif // GMSBTOOLS_GMSBSELECTIONTOOL_H 




