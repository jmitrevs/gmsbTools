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


//#include "egammaAnalysisUtils/checkOQ.h"
//#include "egammaAnalysisUtils/EnergyRescaler.h"

//#include "MuonMomentumCorrections/SmearingClass.h"

#include "GaudiKernel/ToolHandle.h"

#include <string>
#include <map>
#include <vector>


/** Interface ID */  
static const InterfaceID IID_gmsbSelectionTool("gmsbSelectionTool", 1, 0);

class ElectronD3PDObject;
class MuonD3PDObject;
class JetD3PDObject;
class PhotonD3PDObject;

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
  // note: the references are not const since an update can be made to the calibrations
  bool isSelected( ElectronD3PDObject& electron, std::size_t idx ) const; 
  bool isSelected( PhotonD3PDObject& photon, std::size_t idx ) const; 
  bool isSelected( MuonD3PDObject& muon, std::size_t idx ) const;
  bool isSelected( JetD3PDObject& jet, std::size_t idx ) const;

  //  bool isBJet( const Jet * jet ) const;

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
  bool   m_doNewElectronIsolation;
  bool   m_doElectronTrackIsolation;
  bool   m_authorEgammaOnly;
  bool   m_doElectronEtaWindCut;/// apply (or not) eta cut in bad region window   
  bool   m_doEDElectronIsolation;
 
  double m_electronEtaWindMin;
  double m_electronEtaWindMax;
  double m_electronEtcone20corrected; // for new
  double m_electronPtcone20ovEt;  // for track

  /** Photon selection */
  double m_photonPt;
  double m_photonEta;
  int    m_photonID;
  unsigned int m_photonIsEM;
  bool   m_doPhotonTrackIsolation;
  bool   m_doEDPhotonIsolation;
  bool   m_doPhotonEtaWindCut;/// apply (or not) eta cut in bad region window  
  double m_photonEtaWindMin;
  double m_photonEtaWindMax;
  double m_photonPtcone20ovEt; // for track
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

  /** Jet selection */
  double m_jetPt;
  double m_jetEta;
  double m_bJetLikelihood;
  double m_rejNegEJets;		// reject jets with neg E



  // // the OQ utility
  // egammaOQ m_OQ;

  // mutable eg2011::EnergyRescaler m_eRescale;

  // mutable MuonSmear::SmearingClass m_muonSmear;

};

#endif // GMSBTOOLS_GMSBSELECTIONTOOL_H 




