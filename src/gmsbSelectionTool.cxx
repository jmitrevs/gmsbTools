/*****************************************************************************
Name    : gmsbSelectionTool.cxx
Package : offline/PhysicsAnalysis/AnalysisCommon/gmsbTools
Author  : Ketevi A. Assamagan
Created : November 2007
Purpose : User Analysis Selections - see gmsbSelectionTool.h for details
*****************************************************************************/

//#include "GaudiKernel/GaudiException.h"
//#include "GaudiKernel/Property.h"

// Accessing data:
#include "CLHEP/Units/SystemOfUnits.h"
#include "gmsbD3PDObjects/ElectronD3PDObject.h"
#include "gmsbD3PDObjects/MuonD3PDObject.h"
#include "gmsbD3PDObjects/JetD3PDObject.h"
#include "gmsbD3PDObjects/PhotonD3PDObject.h"

// User Tools
#include "gmsbTools/gmsbSelectionTool.h"

//#include "egammaAnalysisUtils/CaloIsoCorrection.h"

#include "ApplyJetCalibration/ApplyJetCalibration.h"

#include "egammaAnalysisUtils/PhotonIDTool.h"

#include "PathResolver/PathResolver.h"

#include "MCTruthClassifier/MCTruthClassifierDefs.h"

#include <sstream>
#include <iomanip>
#include <iostream>


//------------------------------------------------------------------------------
gmsbSelectionTool::gmsbSelectionTool( const std::string& type,
				      const std::string& name, 
				      const IInterface* parent )
  : AthAlgTool( type, name, parent ),
    m_muonSmear("Data12","staco","q_pT","Rel17.2Sum13","")
{
  declareInterface<gmsbSelectionTool>( this );

  declareProperty("Atlfast",          m_isAtlfast=false);

  declareProperty("WhichSyste", m_whichsyste = SystErr::NONE);

  declareProperty("IsMC", m_isMC=false);
  declareProperty("SmearMC", m_smearMC = true);
  declareProperty("MCHasConstantTerm", m_MCHasConstantTerm = false);
  //  declareProperty("RandomSeed", m_randomSeed = 0); // use SUSY prescription
  declareProperty("ElScaleShift", m_elScaleShift = 0);
  declareProperty("ElSmearShift", m_elSmearShift = 0);
  declareProperty("PhoScaleShift", m_phoScaleShift = 0);
  declareProperty("PhoSmearShift", m_phoSmearShift = 0);
  declareProperty("MCEtconeScale", m_mcEtconeScale = 1.5);
  declareProperty("MCUseAltIsoCorrection", m_useAltIsoCorrection = false);


  declareProperty("DoTruth", m_doTruth = false);


  /** Electron selection */
  declareProperty("ElectronPt",       m_electronPt=20*GeV);
  declareProperty("ElectronEta",      m_electronEta=2.47);
  declareProperty("ElectronID", m_electronID = egammaPID::ElectronIDMediumPP);
  declareProperty("AuthorEgammaOnly", m_authorEgammaOnly=true);
  declareProperty("DoElectronEtaWindowCut", m_doElectronEtaWindCut = false);
  declareProperty("ElectronEtaWindowMin", m_electronEtaWindMin = 1.37);
  declareProperty("ElectronEtaWindowMax", m_electronEtaWindMax = 1.52);
  declareProperty("ElectronEtcone20corrected", m_electronEtcone20corrected=5*GeV);
  //declareProperty("DoElectronTrackIsolation", m_doElectronTrackIsolation = false);
  declareProperty("ElectronPtcone20ovEt", m_electronPtcone20ovEt=0.1);
  //declareProperty("DoEDElectronIsolation", m_doEDElectronIsolation = false);
  declareProperty("DoElectronIsolation", m_doElectronIsolation = gmsbSelectionTool::NONE); // isoType
  declareProperty("Simple", m_simple=false); // don't smear or decorate object
                                             // (useful for selecting already selected)

  declareProperty("RescalerData", m_rescalerData = "EnergyRescalerData.root");

  /** Photon selection */
  declareProperty("PhotonPt",   m_photonPt=20*GeV);
  declareProperty("PhotonEta",  m_photonEta=2.37);
  declareProperty("PhotonID", m_photonID = egammaPID::PhotonIDTightAR);
  declareProperty("PhotonIsEM", m_photonIsEM = 0);
  declareProperty("DoPhotonEtaWindowCut", m_doPhotonEtaWindCut = false);
  declareProperty("PhotonEtaWindowMin", m_photonEtaWindMin = 1.37);
  declareProperty("PhotonEtaWindowMax", m_photonEtaWindMax = 1.52);
  declareProperty("PhotonPtcone20ovEt", m_photonPtcone20ovEt=0.1);
  declareProperty("DoPhotonTrackIsolation", m_doPhotonTrackIsolation = false);
  declareProperty("PhotonEtcone20corrected", m_photonEtcone20corrected=5*GeV);
  declareProperty("DoEDPhotonIsolation", m_doEDPhotonIsolation = false);
  declareProperty("ApplyFF", m_doFF = true); // only applied on fullsim MC
  declareProperty("FFSet", m_FFset = 14); //

  /** Muon selection */
  declareProperty("MuonPt", m_muonPt = 6.0*GeV);
  declareProperty("MuonEta", m_muonEta = 2.5);
  //declareProperty("MuonMatchChi2Max", m_matchChi2Max = 150.0);
  declareProperty("SelectCombined", m_sel_combined=true);
  declareProperty("SelectSegmentTag", m_sel_seg_tag=true);
  declareProperty("DoMuonIsolation", m_doMuonIsolation = gmsbSelectionTool::NONE); // isoType
  declareProperty("MuonFlatIsoCut",m_flat_isolation_cut = 1.8*GeV);
  declareProperty("MuonIsoCut",m_isolation_cut = 0.1);
  //declareProperty("MuonSpecPtLimit", m_ms_pt_limit = 50.*GeV);
  //declareProperty("MuonSegMomDiffLimit", m_ms_p_diff_limit = -0.4);
  declareProperty("MuonResSyst", m_muonResSyst = "");


  /** Jet selection */
  declareProperty("JetPt",          m_jetPt=20*GeV);
  declareProperty("JetEta",         m_jetEta=100);
  declareProperty("rejectNegativeEnergyJets", m_rejNegEJets = false); 
  declareProperty("BJetLikelihood", m_bJetLikelihood=6.0);

  declareProperty("FullSimJESConfigFile", 
		  m_fsJESConfigFile = "JES_Full2012dataset_Preliminary_Jan13.config");
  declareProperty("AF2JESConfigFile", 
		  m_af2JESConfigFile = "JES_Full2012dataset_Preliminary_AFII_Jan13.config");

}

//------------------------------------------------------------------------------
StatusCode gmsbSelectionTool::initialize() {

  ATH_MSG_DEBUG("in initialize()");

  if (!m_simple) {
    std::string rescalerData = PathResolver::find_file(m_rescalerData, "DATAPATH");

    if (rescalerData == "") {
      ATH_MSG_ERROR("rescaler data file " << m_rescalerData << " not found. Exiting");
      return StatusCode::FAILURE;
    }

    m_eRescale.Init(rescalerData, "2012", "es2012");
    if (m_doFF && m_isMC && !m_isAtlfast) {
      m_ft.SetPreselection(m_FFset);
    }

    m_jetCalibrator = new JetCalibrationTool("AntiKt4LCTopo", 
					     m_isAtlfast ? m_af2JESConfigFile : m_fsJESConfigFile,
					     !m_isMC);

  } else{
    m_jetCalibrator = 0;
  }
  // m_eRescale.useDefaultCalibConstants("2011");
  // // m_eRescale.SetRandomSeed(m_randomSeed);

  // m_muonSmear.UseScale(1);
  // m_muonSmear.UseImprovedCombine();

  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
StatusCode gmsbSelectionTool::finalize() {

  ATH_MSG_DEBUG("in finalize()");
 
  delete m_jetCalibrator;

  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
gmsbSelectionTool::~gmsbSelectionTool()
{}


float gmsbSelectionTool::calibrate(const ElectronD3PDObject& electron, std::size_t idx) const
{
  const float cleta = electron.cl_eta(idx);
  const float uncorrectedE = electron.cl_E(idx);

  float energy = uncorrectedE;

  if ( m_authorEgammaOnly && !electron.author(idx, egammaParameters::AuthorElectron)) {
    return energy;
  }
  if (m_isMC) {

    /// electron energy scale uncertainty
    switch(m_whichsyste) {
    case SystErr::EGZEEUP:
      energy = m_eRescale.applyEnergyCorrection(cleta, uncorrectedE, 
						egRescaler::EnergyRescalerUpgrade::Electron, 
						egRescaler::EnergyRescalerUpgrade::ZeeAllUp);
      break;

    case SystErr::EGZEEDOWN:
      energy = m_eRescale.applyEnergyCorrection(cleta, uncorrectedE, 
						egRescaler::EnergyRescalerUpgrade::Electron, 
						egRescaler::EnergyRescalerUpgrade::ZeeAllDown);
      break;

    case SystErr::EGMATUP:
      energy = m_eRescale.applyEnergyCorrection(cleta, uncorrectedE, 
						egRescaler::EnergyRescalerUpgrade::Electron, 
						egRescaler::EnergyRescalerUpgrade::R12StatUp);
      break;

    case SystErr::EGMATDOWN:
      energy = m_eRescale.applyEnergyCorrection(cleta, uncorrectedE, 
						egRescaler::EnergyRescalerUpgrade::Electron, 
						egRescaler::EnergyRescalerUpgrade::R12StatDown);
      break;

    case SystErr::EGPSUP:
      energy = m_eRescale.applyEnergyCorrection(cleta, uncorrectedE, 
						egRescaler::EnergyRescalerUpgrade::Electron,
						egRescaler::EnergyRescalerUpgrade::PSStatUp);
      break;

    case SystErr::EGPSDOWN:
      energy = m_eRescale.applyEnergyCorrection(cleta, uncorrectedE, 
						egRescaler::EnergyRescalerUpgrade::Electron, 
						egRescaler::EnergyRescalerUpgrade::PSStatDown);
      break;

    case SystErr::EGLOWUP:
      energy = m_eRescale.applyEnergyCorrection(cleta, uncorrectedE, 
						egRescaler::EnergyRescalerUpgrade::Electron, 
						egRescaler::EnergyRescalerUpgrade::LowPtUp);
      break;

    case SystErr::EGLOWDOWN: 
      energy = m_eRescale.applyEnergyCorrection(cleta, uncorrectedE, 
						egRescaler::EnergyRescalerUpgrade::Electron,
						egRescaler::EnergyRescalerUpgrade::LowPtDown);
      break;
    default:
      break;
    }

    if (m_smearMC) {
      /// electron energy smearing
      int seed = int(1.e+5*fabs(electron.cl_phi(idx)));
      if(!seed) ++seed;
      m_eRescale.SetRandomSeed(seed); 
      
      double smearcorr = 1;
      if (m_whichsyste == SystErr::EGRESUP) {
	smearcorr = m_eRescale.getSmearingCorrection(cleta, uncorrectedE, 
						     egRescaler::EnergyRescalerUpgrade::ERR_UP);
      } else if (m_whichsyste == SystErr::EGRESDOWN) {
	smearcorr = m_eRescale.getSmearingCorrection(cleta, uncorrectedE, 
						     egRescaler::EnergyRescalerUpgrade::ERR_DOWN);
      } else {
	smearcorr = m_eRescale.getSmearingCorrection(cleta, uncorrectedE, 
						     egRescaler::EnergyRescalerUpgrade::NOMINAL);
      }
      energy *= smearcorr;
    }
    /// Atlfast specific calibration corrections
    if (m_isAtlfast) {
      energy *= m_eRescale.applyAFtoG4(cleta);
    }
  } else { /// Residual energy scale corrections to be applied to data
    energy = m_eRescale.applyEnergyCorrection(cleta, uncorrectedE, 
					      egRescaler::EnergyRescalerUpgrade::Electron, 
					      egRescaler::EnergyRescalerUpgrade::Nominal);
  }
  return energy;
}


bool gmsbSelectionTool::isSelected( ElectronD3PDObject& electron, std::size_t idx, int nPV ) const
{
  ATH_MSG_DEBUG("in electron isSelected(), with electron = " << idx);

  bool select = true;

  if ( m_authorEgammaOnly ) select = select && electron.author(idx, egammaParameters::AuthorElectron);

  if (!select) return false;

  if (m_isMC && m_doTruth) {
    if (electron.type(idx) != MCTruthPartClassifier::IsoElectron) return false;
  }
  

  const float eta2 = electron.etas2(idx);
  const float absClusEta = fabsf(eta2);


  const float eta = electron.tracketa(idx);
  const float phi = electron.trackphi(idx);

  const float uncorrectedE = electron.cl_E(idx);
  //const float uncorrectedEt = uncorrectedE/cosh(eta);

  float energy = uncorrectedE;

  if (!m_simple) {

    energy = calibrate(electron, idx);

    ATH_MSG_DEBUG("Original electron E = " << uncorrectedE << ", corrected E = " << energy);

    electron.E(idx) = energy;
    electron.eta(idx) = eta;
    electron.phi(idx) = phi;
    electron.pt(idx) = energy/cosh(eta); // massless
  }

  const float pt = electron.pt(idx);

  select = electron.passID(idx, static_cast<egammaPID::egammaIDQuality>(m_electronID)); 

  select = select && pt > m_electronPt && absClusEta <m_electronEta;
  ATH_MSG_DEBUG("select is now " << select);

  // check if electron is in bad eta region

  if(m_doElectronEtaWindCut) {
    const bool isCrack = absClusEta > m_electronEtaWindMin && absClusEta < m_electronEtaWindMax; 
    select = select && !isCrack;
  }

  ATH_MSG_DEBUG("after crack, select is now " << select);

  // check OQ
  bool badOQ = (electron.OQ(idx) & egammaPID::BADCLUSELECTRON); // 0 == good
  // if (m_isMC) {
  //   badOQ = badOQ || 
  //     m_OQ.checkOQClusterElectron(runNum, electron->cluster()->eta(), electron->cluster()->phi())==3;
  // }
  select = select && !badOQ;

  switch (m_doElectronIsolation) {
  case TrackIso:
    {
      const float isol = electron.ptcone20_zpv05(idx) / pt;
      select = select && isol < m_electronPtcone20ovEt;
    }
    break;
  case EDIso:
    {
      const float isol = electron.topoEtcone20_corrected(idx);
      select = select && isol < m_electronEtcone20corrected;
    }
    break;
  case LooserIso:
    {
      float id_isocut = 0.16;
      id_isocut *= pt;
      select = select && electron.ptcone30(idx) < id_isocut;
    }
    break;
  case LooseIso:
    {
      float id_isocut = 0.16;
      const float z0cut = 0.4;
      id_isocut *= pt;
      select = select && electron.ptcone30(idx) < id_isocut
	&& fabs(electron.trackIPEstimate_z0_unbiasedpvunbiased(idx)*sin(electron.tracktheta(idx))) <= z0cut;
    }
    break;
  case MediumIso:
    {
      float id_isocut = 0.16;
      const float d0sigcut = 5.;
      const float z0cut = 0.4;
      id_isocut *= pt;
      select = select && electron.ptcone30(idx) < id_isocut
	&& fabs(electron.trackIPEstimate_z0_unbiasedpvunbiased(idx)*sin(electron.tracktheta(idx))) <= z0cut;

      if (electron.trackIPEstimate_sigd0_unbiasedpvunbiased(idx) != 0) {
	select = select && fabs(electron.trackIPEstimate_d0_unbiasedpvunbiased(idx)/electron.trackIPEstimate_sigd0_unbiasedpvunbiased(idx)) <= d0sigcut;
      }
    }
    break;
  case TightIso:
    {
      float id_isocut = 0.16;
      float calo_isocut = 0.18;
      const float d0sigcut = 5.;
      const float z0cut = 0.4;
      const float corr = m_isMC ? 17.94 : 20.15;
      const float calo_iso = electron.topoEtcone30_corrected(idx) - corr * nPV;

      id_isocut *= pt;
      calo_isocut *= pt;     

      select = select && electron.ptcone30(idx) < id_isocut
	&& fabs(electron.trackIPEstimate_z0_unbiasedpvunbiased(idx)*sin(electron.tracktheta(idx))) <= z0cut
	&& calo_iso < calo_isocut;

      if (electron.trackIPEstimate_sigd0_unbiasedpvunbiased(idx) != 0) {
	select = select && fabs(electron.trackIPEstimate_d0_unbiasedpvunbiased(idx)/electron.trackIPEstimate_sigd0_unbiasedpvunbiased(idx)) <= d0sigcut;
      }
    }
    break;
  case NONE:
    break; // do nothing
  default:
    ATH_MSG_ERROR("Using unimplemented electron isolation: " << m_doElectronIsolation);
    return StatusCode::FAILURE;
  }

  ATH_MSG_DEBUG("after iso, select is now " << select);

  return select;
}

float gmsbSelectionTool::calibrate(const PhotonD3PDObject& photon, std::size_t idx) const
{

  const float uncorrectedE = photon.cl_E(idx);
  float energy = uncorrectedE;

  const float cleta = photon.cl_eta(idx);

  egRescaler::EnergyRescalerUpgrade::ParticleType particleType = 
    (photon.isConv(idx)) ? 
    egRescaler::EnergyRescalerUpgrade::Converted : 
    egRescaler::EnergyRescalerUpgrade::Unconverted;
  
  if (m_isMC) {
      
    /// electron energy scale uncertainty
    switch(m_whichsyste) {
    case SystErr::EGZEEUP:
      energy = m_eRescale.applyEnergyCorrection(cleta, uncorrectedE, 
						particleType, 
						egRescaler::EnergyRescalerUpgrade::ZeeAllUp);
      break;
      
    case SystErr::EGZEEDOWN:
      energy = m_eRescale.applyEnergyCorrection(cleta, uncorrectedE, 
						particleType, 
						egRescaler::EnergyRescalerUpgrade::ZeeAllDown);
      break;
      
    case SystErr::EGMATUP:
      energy = m_eRescale.applyEnergyCorrection(cleta, uncorrectedE, 
						particleType, 
						egRescaler::EnergyRescalerUpgrade::R12StatUp);
      break;
      
    case SystErr::EGMATDOWN:
      energy = m_eRescale.applyEnergyCorrection(cleta, uncorrectedE, 
						particleType, 
						egRescaler::EnergyRescalerUpgrade::R12StatDown);
      break;
      
    case SystErr::EGPSUP:
      energy = m_eRescale.applyEnergyCorrection(cleta, uncorrectedE, 
						particleType,
						egRescaler::EnergyRescalerUpgrade::PSStatUp);
      break;
      
    case SystErr::EGPSDOWN:
      energy = m_eRescale.applyEnergyCorrection(cleta, uncorrectedE, 
						particleType, 
						egRescaler::EnergyRescalerUpgrade::PSStatDown);
      break;
      
    case SystErr::EGLOWUP:
      energy = m_eRescale.applyEnergyCorrection(cleta, uncorrectedE, 
						particleType, 
						egRescaler::EnergyRescalerUpgrade::LowPtUp);
      break;
      
    case SystErr::EGLOWDOWN: 
      energy = m_eRescale.applyEnergyCorrection(cleta, uncorrectedE, 
						particleType,
						egRescaler::EnergyRescalerUpgrade::LowPtDown);
      break;
    default:
      break;
    }
    
    if (m_smearMC) {
      /// electron energy smearing
      int seed = int(1.e+5*fabs(photon.cl_phi(idx)));
      if(!seed) ++seed;
      m_eRescale.SetRandomSeed(seed); 
      
      double smearcorr = 1;
      if (m_whichsyste == SystErr::EGRESUP) {
	smearcorr = m_eRescale.getSmearingCorrection(cleta, uncorrectedE, 
						     egRescaler::EnergyRescalerUpgrade::ERR_UP);
      } else if (m_whichsyste == SystErr::EGRESDOWN) {
	smearcorr = m_eRescale.getSmearingCorrection(cleta, uncorrectedE, 
						     egRescaler::EnergyRescalerUpgrade::ERR_DOWN);
      } else {
	smearcorr = m_eRescale.getSmearingCorrection(cleta, uncorrectedE, 
						     egRescaler::EnergyRescalerUpgrade::NOMINAL);
      }
      energy *= smearcorr;
    }
    /// Atlfast specific calibration corrections
    if (m_isAtlfast) {
      energy *= m_eRescale.applyAFtoG4(cleta);
    }
  } else { /// Residual energy scale corrections to be applied to data
    energy = m_eRescale.applyEnergyCorrection(cleta, uncorrectedE, 
					      particleType, 
					      egRescaler::EnergyRescalerUpgrade::Nominal);
  }
  return energy;
}

bool gmsbSelectionTool::isSelected( PhotonD3PDObject& photon, std::size_t idx ) const
{
  ATH_MSG_DEBUG("in photon isSelected()");

  bool select = false;

    
  if (m_isMC && m_doTruth) {
    if (photon.type(idx) != MCTruthPartClassifier::IsoPhoton) return false;
  }


  const float eta2 = photon.etas2(idx);
  const float absClusEta = fabs(eta2);

  const float uncorrectedE = photon.cl_E(idx);


  float energy = uncorrectedE;

  if (!m_simple) {

    energy = calibrate(photon, idx);
    ATH_MSG_DEBUG("Original photon E = " << photon.E(idx) << ", corrected E = " << energy);

    const float scale = energy / photon.E(idx);

    ATH_MSG_DEBUG("original pT = " <<  photon.pt(idx) << ", final pT = " << photon.pt(idx) * scale);

    photon.E(idx) = energy;
    photon.pt(idx) *= scale; 

  }

  const float pt = photon.pt(idx);

  select = pt > m_photonPt && absClusEta < m_photonEta; 

  ATH_MSG_DEBUG("isEM before fudgind = " << std::hex << photon.isEM(idx) << std::dec << " and passID = " << photon.passID(idx, static_cast<egammaPID::egammaIDQuality>(m_photonID)));

  if (m_doFF && !m_simple) {
    fudgeID(photon, idx);
  } 

    ATH_MSG_DEBUG("isEM after fudging = " << std::hex << photon.isEM(idx) << ", mask = " << egammaPID::PhotonTightAR << std::dec << " and passID = " << photon.passID(idx, static_cast<egammaPID::egammaIDQuality>(m_photonID)));

  select = select && photon.passID(idx, static_cast<egammaPID::egammaIDQuality>(m_photonID));

  // also do isEM selection for special requirements (usually m_photonIsEM == 0, so this does nothing)
  select = select && ((photon.isEM(idx) & m_photonIsEM) == 0);

  ATH_MSG_DEBUG("after pt, eta (= " << eta2 << "), and isEM cut, select = " << select);

  // check if photon is in bad eta region

  if(m_doPhotonEtaWindCut) {
    const bool isCrack = absClusEta > m_photonEtaWindMin && absClusEta < m_photonEtaWindMax; 
    select = select && !isCrack;
  }

  ATH_MSG_DEBUG("after crack cut, select = " << select);

  // check OQ
  bool badOQ = (photon.OQ(idx) & egammaPID::BADCLUSPHOTON); // 0 == good
  // if (m_isMC) {
  //   badOQ = badOQ || 
  //     m_OQ.checkOQClusterPhoton(runNum, photon->cluster()->eta(), photon->cluster()->phi(), photon->conversion())==3;
  // }
  select = select && !badOQ;

  ATH_MSG_DEBUG("after OTX, select = " << select);

  if ( m_doPhotonTrackIsolation ) {
    const float isol = photon.ptcone20_zpv05(idx) / pt;
    select = select && isol < m_photonPtcone20ovEt;
  }

  if ( m_doEDPhotonIsolation ) {
    const float isol = photon.topoEtcone20_corrected(idx);
    select = select && isol < m_photonEtcone20corrected;
  }

  ATH_MSG_DEBUG("photon select = " << select);

  return select;
}

float gmsbSelectionTool::calibrate(const MuonD3PDObject& muon, std::size_t idx) const
{
  float pt = muon.pt(idx);
  if (m_isMC && m_smearMC) {
    int seed = int(fabs(muon.phi(idx)*1.e+5));
    if (!seed) seed = 1;
    ATH_MSG_INFO("seed = " << seed);

    m_muonSmear.SetSeed(seed);
    
    // float charge = muon->charge();
    const float eta = muon.eta(idx);
    const float ptcb = pt;
    const float ptid = (muon.id_qoverp_exPV(idx) != 0.) ? 
      fabs(sin(muon.id_theta_exPV(idx))/muon.id_qoverp_exPV(idx)) : 0.;
    const float ptms = (muon.me_qoverp_exPV(idx) != 0.) ? 
      fabs(sin(muon.me_theta_exPV(idx))/muon.me_qoverp_exPV(idx)) : 0.;
    
    if (muon.isCombinedMuon(idx)) {
      m_muonSmear.Event(ptms,ptid,ptcb,eta,muon.charge(idx),muon.phi(idx)); 
      ATH_MSG_INFO("m_muonSmear.Event(" << ptms << ", " << ptid << " ," << ptcb << ", " 
		   << eta << ", " << muon.charge(idx) << ");");
    }
    else if (muon.isSegmentTaggedMuon(idx))
      m_muonSmear.Event(ptid,eta,"ID",muon.charge(idx),muon.phi(idx)); 
    else
      m_muonSmear.Event(ptms,eta,"MS",muon.charge(idx),muon.phi(idx)); 
    
    std::string THESTRING = "";
    bool doSyst = false;
    switch(m_whichsyste) {
    case SystErr::MMSLOW: THESTRING = "MSLOW";  doSyst = true; break;
    case SystErr::MMSUP:   THESTRING = "MSUP";  doSyst = true;  break;
    case SystErr::MIDLOW: THESTRING = "IDLOW"; doSyst = true; break;
    case SystErr::MIDUP:   THESTRING = "IDUP"; doSyst = true;  break;
    case SystErr::MSCALELOW: THESTRING = "SCALELOW"; doSyst = true; break;
    case SystErr::MSCALEUP:   THESTRING = "SCALEUP"; doSyst = true;  break;
    default: break;
    }
    if (!doSyst) {
      if (muon.isCombinedMuon(idx))
	pt = m_muonSmear.pTCB();
      else if (muon.isSegmentTaggedMuon(idx))
	pt = m_muonSmear.pTID();
      else
	pt = m_muonSmear.pTMS();
    } else {
      double pTMS_smeared = pt;
      double pTID_smeared = pt;
      double pTCB_smeared = pt;
      
      // Valid values for "THESTRING":{"IDLOW", IDUP", "MSLOW", MSUP","SCALELOW", "SCALEUP"} 
      m_muonSmear.PTVar(pTMS_smeared, pTID_smeared, pTCB_smeared, THESTRING);
      
      if (muon.isCombinedMuon(idx)) pt = pTCB_smeared;
      else if (muon.isSegmentTaggedMuon(idx)) pt = pTID_smeared;
      else pt = pTMS_smeared;    
    }
  }
  return pt;
}

bool gmsbSelectionTool::isSelected( MuonD3PDObject& muon, std::size_t idx, int nPV ) const
{
  ATH_MSG_INFO("in muon isSelected(), with muon = " << idx << ", pt = " << muon.pt(idx) << ", eta = " << muon.eta(idx));

  // do ID cut
  bool select = ((m_sel_combined && muon.isCombinedMuon(idx)) ||
		 (m_sel_seg_tag && muon.isSegmentTaggedMuon(idx)));
  
  select = select && muon.loose(idx);

  if (!select) return false;

  ATH_MSG_INFO("Pass type and loose; isCombined = " << muon.isCombinedMuon(idx));

  //float pt = (muon->isCombinedMuon()) ? muon->pt() : muon->inDetTrackParticle()->pt(); ;
  float pt = muon.pt(idx);

  // ATH_MSG_DEBUG("Here 2");
  if (!m_simple) {
    pt = calibrate(muon, idx);
    muon.pt(idx) = pt;
  }
    
  // select must be true before in order to get here, so can overwrite
  select = pt >m_muonPt && fabs(muon.eta(idx))<m_muonEta;

  ATH_MSG_INFO("calibrated pt = " << pt);

  // do iso cut
  
  switch (m_doMuonIsolation) {
  case FlatTrackIso:
    select = select && muon.ptcone20(idx)<m_flat_isolation_cut;
    break;
  case TrackIso:
    select = select && muon.ptcone20(idx)/pt<m_isolation_cut;
    break;
  case LooserIso:
    {
      float id_isocut = 0.12;
      id_isocut *= pt;

      select = select && muon.ptcone30_trkelstyle(idx) < id_isocut;
    }
    break;
  case LooseIso:
    {
      float id_isocut = 0.12;
      const float z0cut = 0.4;
      id_isocut *= pt;

      select = select && muon.ptcone30_trkelstyle(idx) < id_isocut
	&& fabs(muon.trackIPEstimate_z0_unbiasedpvunbiased(idx)*sin(muon.tracktheta(idx))) <= z0cut;
    }
    break;
  case MediumIso:
    {
      float id_isocut = 0.12;
      const float d0sigcut = 3.;
      const float z0cut = 0.4;
      id_isocut *= pt;
      select = select && muon.ptcone30_trkelstyle(idx) < id_isocut
	&& fabs(muon.trackIPEstimate_z0_unbiasedpvunbiased(idx)*sin(muon.tracktheta(idx))) <= z0cut;

      if (muon.trackIPEstimate_sigd0_unbiasedpvunbiased(idx) != 0) {
	select = select && fabs(muon.trackIPEstimate_d0_unbiasedpvunbiased(idx)/muon.trackIPEstimate_sigd0_unbiasedpvunbiased(idx)) <= d0sigcut;
      }
    }
    break;
  case TightIso:
    {
      float id_isocut = 0.12;
      float calo_isocut = 0.12;
      const float d0sigcut = 3.;
      const float z0cut = 0.4;
      const float corr = m_isMC ? 69.2 : 64.8;
      const float corr2 = m_isMC ? 0.76 : 0.98;
      const float calo_iso = muon.etcone30(idx) - corr * nPV - corr2 * nPV * nPV;

      id_isocut *= pt;
      calo_isocut *= pt;     

      select = select && muon.ptcone30_trkelstyle(idx) < id_isocut
	&& fabs(muon.trackIPEstimate_z0_unbiasedpvunbiased(idx)*sin(muon.tracktheta(idx))) <= z0cut
	&& calo_iso < calo_isocut;

      if (muon.trackIPEstimate_sigd0_unbiasedpvunbiased(idx) != 0) {
	select = select && fabs(muon.trackIPEstimate_d0_unbiasedpvunbiased(idx)/muon.trackIPEstimate_sigd0_unbiasedpvunbiased(idx)) <= d0sigcut;
      }
    }
    break;
  case NONE:
    break; // do nothing
  default:
    ATH_MSG_ERROR("Using unimplemented muon isolation: " << m_doMuonIsolation);
    return StatusCode::FAILURE;
  }

  if (!select) return false;

  ATH_MSG_INFO("pass iso");

  // do track cuts
  //if (muon.expectBLayerHit(idx) && muon.nBLHits(idx) == 0) return false;
  if (muon.nPixHits(idx) + muon.nPixelDeadSensors(idx) < 1) return false;
  if (muon.nSCTHits(idx) + muon.nSCTDeadSensors(idx) < 5) return false;
  if (muon.nPixHoles(idx) + muon.nSCTHoles(idx) >= 3) return false;

  ATH_MSG_INFO("pass silicon req");

  const int nTRTOutliers = muon.nTRTOutliers(idx);
  const int nTRTTotal    = nTRTOutliers + muon.nTRTHits(idx);
  const float trackEta   = fabs(log(tan(muon.tracktheta(idx)/2.0)));
  if (trackEta < 1.9 && trackEta > 0.1) {
    select = (nTRTTotal > 5 &&  nTRTOutliers < 0.9 * nTRTTotal);
  } 
  
  ATH_MSG_INFO("Muon with pt = " << pt << ", select = " << select);
  
  return select;
}

TLorentzVector gmsbSelectionTool::calibrate(const JetD3PDObject& jet, std::size_t idx, 
					    float rhoKt4LC, float mu, int nPV2) const
{
  return m_jetCalibrator->ApplyJetAreaOffsetEtaJES(jet.constscale_E(idx),
						   jet.constscale_eta(idx),
						   jet.constscale_phi(idx),
						   jet.constscale_m(idx),
						   jet.ActiveAreaPx(idx),  
						   jet.ActiveAreaPy(idx), 
						   jet.ActiveAreaPz(idx),  
						   jet.ActiveAreaE(idx), 
						   rhoKt4LC,
						   mu,
						   nPV2);
}


bool gmsbSelectionTool::isSelected( JetD3PDObject& jet, std::size_t idx, 
				    float rhoKt4LC, float mu, int nPV2) const
{

  if (!m_simple) {

    // ATH_MSG_WARNING("calling isSelected with rhoKt4LC = " << rhoKt4LC << ", mu = " 
    // 		    <<  mu << ", nPV2 = " << nPV2);

    TLorentzVector jlv = calibrate(jet, idx, rhoKt4LC, mu, nPV2);

    jet.E(idx) = jlv.E();
    jet.eta(idx) = jlv.Eta();
    jet.phi(idx) = jlv.Phi();
    jet.pt(idx) = jlv.Pt();
  }

  if (m_rejNegEJets && jet.E(idx) < 0) {
    return false;
  }

  return (jet.pt(idx) > m_jetPt && fabs(jet.eta(idx)) < m_jetEta);

}


// bool gmsbSelectionTool::isBJet( const Jet * jet ) const {

//   /** first check that it is a selected jet */
//   if ( !jet ) return false;

//   bool is_bjet = this->isSelected( jet );
//   return ( is_bjet && jet->getFlavourTagWeight()>m_bJetLikelihood );
// }


void gmsbSelectionTool::fudgeID(PhotonD3PDObject& photon, 
				std::size_t idx) const
{

  // calculate some veraibles

  const float eta2 = photon.etas2(idx);
  float et = 0;

  if (fabs(eta2)<999.)
    et = cosh(eta2)!=0. ? photon.cl_E(idx)/cosh(eta2) : 0.;

  float rhad1 = fabsf(et)>0. ? photon.Ethad1(idx)/et : 0.;
  float rhad  = fabsf(et)>0. ? photon.Ethad(idx)/et : 0.;


  float reta = fabsf(photon.E277(idx))>0. ? photon.E237(idx)/photon.E277(idx) : 0.;
  float rphi = fabsf(photon.E237(idx))>0. ? photon.E233(idx)/photon.E237(idx) : 0.;

  const float emax2 = photon.Emax2(idx);
  const float emin = photon.Emins1(idx);
  const float emaxs1 = photon.emaxs1(idx);

  float deltae = emax2 - emin;  
  float eratio = (emaxs1+emax2)==0. ? 0 : (emaxs1 - emax2)/(emaxs1+emax2);

  ATH_MSG_DEBUG("fudging with et = " << et << ", eta2 = " << eta2 
		<< ", reta = " << reta << ", rphi = " << rphi 
		<< ", deltae = " << deltae <<  ", eratio = " << eratio 
		<< ", emaxs1 = " << emaxs1 << ", emax2 = " << emax2
		<< ", rhad1 = " << rhad1 << ", rhad = " << rhad
		<< ", wstot = " << photon.wstot(idx)
		<< ", ws3 = " << photon.ws3(idx)
		<< ", isConv = " << photon.isConv(idx)); 

  if (m_isMC && !m_isAtlfast) {
    m_ft.FudgeShowers(photon.pt(idx),
		      eta2,
		      rhad1,
		      rhad,
		      photon.E277(idx),
		      reta,
		      rphi,
		      photon.weta2(idx),
		      photon.f1(idx),
		      photon.fside(idx),
		      photon.wstot(idx),
		      photon.ws3(idx),
		      deltae,
		      eratio,
		      photon.isConv(idx));
  }
  ATH_MSG_DEBUG("after fudge"
		<< " reta = " << reta << ", rphi = " << rphi 
		<< ", deltae = " << deltae <<  ", eratio = " << eratio 
		<< ", emaxs1 = " << emaxs1 << ", emax2 = " << emax2
		<< ", rhad1 = " << rhad1 << ", rhad = " << rhad
		<< ", wstot = " << photon.wstot(idx)
		<< ", ws3 = " << photon.ws3(idx)
		<< ", isConv = " << photon.isConv(idx)); 

  PhotonIDTool myTool(et,
		      eta2,
		      rhad1,
		      rhad,
		      photon.E277(idx),
		      reta,
		      rphi,
		      photon.weta2(idx),
		      photon.f1(idx),
		      photon.fside(idx),
		      photon.wstot(idx),
		      photon.ws3(idx),
		      deltae,
		      eratio,
		      photon.isConv(idx));

  const unsigned int fudgedIsEM = myTool.isEM(4,2012);

  const unsigned int origIsEM = photon.isEM(idx);

  const unsigned int newIsEM = (origIsEM & egammaPID::AMBIGUITYRESOLVE_PHOTON) | fudgedIsEM;

  photon.isEM(idx) = newIsEM;

}

