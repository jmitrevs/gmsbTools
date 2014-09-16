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

#include "egammaEvent/egammaPIDdefs.h"

//#include "egammaAnalysisUtils/checkOQ.h"
#include "egammaAnalysisUtils/EnergyRescalerUpgrade.h"
#include "egammaAnalysisUtils/FudgeMCTool.h"

#include "MuonMomentumCorrections/SmearingClass.h"

#include "gmsbTools/gmsbSystError.h"

#include "GaudiKernel/ToolHandle.h"

#include <string>
#include <map>
#include <vector>

#include <TLorentzVector.h>

/** Interface ID */  
static const InterfaceID IID_gmsbSelectionTool("gmsbSelectionTool", 1, 0);

class ElectronD3PDObject;
class MuonD3PDObject;
class JetD3PDObject;
class PhotonD3PDObject;
class JetCalibrationTool;

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

  float calibrate(const ElectronD3PDObject& electron, std::size_t idx) const; // return is energy
  float calibrate(const PhotonD3PDObject& photon, std::size_t idx) const; // return is energy
  float calibrate(const MuonD3PDObject& muon, std::size_t idx) const; // return is pt
  TLorentzVector calibrate(const JetD3PDObject& jet, std::size_t idx, 
			   float rhoKt4LC, float mu, int nPV2) const; // return is TLV

  /** pre-selections */
  // note: the references are not const since an update can be made to the calibrations
  bool isSelected( ElectronD3PDObject& electron, std::size_t idx, int nPV) const; 
  bool isSelected( PhotonD3PDObject& photon, std::size_t idx ) const; 
  bool isSelected( MuonD3PDObject& muon, std::size_t idx, int nPV) const;
  bool isSelected( JetD3PDObject& jet, std::size_t idx, 
		   float rhoKt4LC = -999.0, float mu = -999.0, int nPV2 = -999) const;

  //  bool isBJet( const Jet * jet ) const;

  enum isoType {NONE = 0, EDIso, TrackIso, FlatTrackIso, LooseIso, MediumIso, TightIso, LooserIso};

protected:

   /** Standard destructor */
   virtual ~gmsbSelectionTool();

  void fudgeID(PhotonD3PDObject& photon, std::size_t idx) const;

private:

  /** this is Atlfast data */
  bool m_isAtlfast;

  bool m_simple; 		// if true, don't smear or decorate event

  int m_whichsyste;  // really SystErr::Syste


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
  float m_mcEtconeScale;
  bool m_useAltIsoCorrection;

  /** Electron selection */
  float m_electronPt;
  float m_electronEta;
  int    m_electronID;
  bool   m_authorEgammaOnly;
  bool   m_doElectronEtaWindCut;/// apply (or not) eta cut in bad region window   
   float m_electronEtaWindMin;
  float m_electronEtaWindMax;
  int   m_doElectronIsolation; // enum value from above
  //bool   m_doEDElectronIsolation;
  //bool   m_doElectronTrackIsolation;
  float m_electronEtcone20corrected; // for new
  float m_electronPtcone20ovEt;  // for track

  std::string m_rescalerData; // name of the rescalar root file

  /** Photon selection */
  float m_photonPt;
  float m_photonEta;
  int    m_photonID;
  unsigned int m_photonIsEM;
  bool   m_doPhotonTrackIsolation;
  bool   m_doEDPhotonIsolation;
  bool   m_doPhotonEtaWindCut;/// apply (or not) eta cut in bad region window  
  float m_photonEtaWindMin;
  float m_photonEtaWindMax;
  float m_photonPtcone20ovEt; // for track
  float m_photonEtcone20corrected;  // for new
  bool m_doFF;
  int m_FFset;

  /** Muon selection */
  float m_muonPt;
  float m_muonEta;
  // float m_matchChi2Max;
  bool m_sel_combined;
  bool m_sel_seg_tag;
  int m_doMuonIsolation;
  float m_flat_isolation_cut;
  float m_isolation_cut;
  //float m_ms_pt_limit;
  //float m_ms_p_diff_limit;
  std::string m_muonResSyst;  // Valid values: {"MSLOW", "MSUP", "IDLOW", "IDUP"} 

  /** Jet selection */
  float m_jetPt;
  float m_jetEta;
  float m_bJetLikelihood;
  float m_rejNegEJets;		// reject jets with neg E

  std::string m_fsJESConfigFile;
  std::string m_af2JESConfigFile;

  // // the OQ utility
  // egammaOQ m_OQ;

  // the photon utilities
  mutable FudgeMCTool  m_ft;

  mutable egRescaler::EnergyRescalerUpgrade m_eRescale;

  mutable MuonSmear::SmearingClass m_muonSmear;

  JetCalibrationTool *m_jetCalibrator;

};


#endif // GMSBTOOLS_GMSBSELECTIONTOOL_H 




