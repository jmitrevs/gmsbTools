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
#include "CLHEP/Units/PhysicalConstants.h"
#include "gmsbD3PDObjects/ElectronD3PDObject.h"
#include "gmsbD3PDObjects/MuonD3PDObject.h"
#include "gmsbD3PDObjects/JetD3PDObject.h"
#include "gmsbD3PDObjects/PhotonD3PDObject.h"

// User Tools
#include "gmsbTools/gmsbSelectionTool.h"

//#include "egammaAnalysisUtils/CaloIsoCorrection.h"

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
  : AthAlgTool( type, name, parent )
    //    m_muonSmear("Data11","staco","q_pT","Rel17","")
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
  declareProperty("ElectronPt",       m_electronPt=25*GeV);
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
  declareProperty("PhotonPt",   m_photonPt=25*GeV);
  declareProperty("PhotonEta",  m_photonEta=2.37);
  declareProperty("PhotonID", m_photonID = egammaPID::PhotonIDTightAR);
  declareProperty("PhotonIsEM", m_photonIsEM = 0);
  declareProperty("DoPhotonEtaWindowCut", m_doPhotonEtaWindCut = false);
  declareProperty("PhotonEtaWindowMin", m_photonEtaWindMin = 1.37);
  declareProperty("PhotonEtaWindowMax", m_photonEtaWindMax = 1.52);
  declareProperty("PhotonPtcone20ovEt", m_photonPtcone20ovEt=0.1);
  declareProperty("DoPhotonTrackIsolation", m_doPhotonTrackIsolation = false);
  declareProperty("PhotonEtcone20corrected", m_photonEtcone20corrected=5*GeV);
  declareProperty("DoEDPhotonIsolation", m_doEDPhotonIsolation = true);
  declareProperty("ApplyFF", m_doFF = true); // only applied on fullsim MC
  declareProperty("FFSet", m_FFset = 14); //

  /** Muon selection */
  declareProperty("MuonPt", m_muonPt = 6.0*GeV);
  declareProperty("MuonEta", m_muonEta = 2.5);
  //declareProperty("MuonMatchChi2Max", m_matchChi2Max = 150.0);
  declareProperty("SelectCombined", m_sel_combined=true);
  declareProperty("SelectSegmentTag", m_sel_seg_tag=true);
  declareProperty("DoMuonIsoCut",m_do_iso_cut = false);
  declareProperty("DoMuonFlatIsoCut",m_do_flat_iso_cut = true);
  declareProperty("MuonFlatIsoCut",m_flat_isolation_cut = 1.8*GeV);
  declareProperty("MuonIsoCut",m_isolation_cut = 0.1);
  //declareProperty("MuonSpecPtLimit", m_ms_pt_limit = 50.*GeV);
  //declareProperty("MuonSegMomDiffLimit", m_ms_p_diff_limit = -0.4);
  declareProperty("MuonResSyst", m_muonResSyst = "");


  /** Jet selection */
  declareProperty("JetPt",          m_jetPt=20*GeV);
  declareProperty("JetEta",         m_jetEta=100);
  declareProperty("rejectNegativeEnergyJets", m_rejNegEJets = true); 
  declareProperty("BJetLikelihood", m_bJetLikelihood=6.0);

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
 
  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
gmsbSelectionTool::~gmsbSelectionTool()
{}

// pt = -1 means use real one
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

  const float cleta = electron.cl_eta(idx);

  const float eta = electron.tracketa(idx);
  const float phi = electron.trackphi(idx);

  const float uncorrectedE = electron.cl_E(idx);
  //const float uncorrectedEt = uncorrectedE/cosh(eta);

  float energy = uncorrectedE;

  if (!m_simple) {
    
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
	&& fabs(electron.trackIPEstimate_z0_unbiasedpvunbiased(idx)*sin(electron.tracktheta(idx))) <= z0cut
	&& fabs(electron.trackIPEstimate_d0_unbiasedpvunbiased(idx)/electron.trackIPEstimate_sigd0_unbiasedpvunbiased(idx)) <= d0sigcut;
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
	&& fabs(electron.trackIPEstimate_d0_unbiasedpvunbiased(idx)/electron.trackIPEstimate_sigd0_unbiasedpvunbiased(idx)) <= d0sigcut
	&& calo_iso < calo_isocut;
    }
    break;
  default:
    break;
  }

  ATH_MSG_DEBUG("after iso, select is now " << select);

  return select;
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

  const float cleta = photon.cl_eta(idx);

  egRescaler::EnergyRescalerUpgrade::ParticleType particleType = 
    (photon.convFlag(idx)) ? 
    egRescaler::EnergyRescalerUpgrade::Converted : 
    egRescaler::EnergyRescalerUpgrade::Unconverted;

  float energy = uncorrectedE;

  if (!m_simple) {
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
    
    ATH_MSG_DEBUG("Original photon E = " << photon.E(idx) << ", corrected E = " << energy);

    const float scale = energy / photon.E(idx);

    photon.E(idx) = energy;
    photon.pt(idx) *= scale; 

  }

  const float pt = photon.pt(idx);

  select = pt > m_photonPt && absClusEta < m_photonEta; 

  if (m_doFF && m_isMC && !m_isAtlfast && !m_simple) {
    fudgeID(photon, idx);
  } 

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

bool gmsbSelectionTool::isSelected( MuonD3PDObject& muon, std::size_t idx, int nPV ) const
{
  ATH_MSG_DEBUG("in muon isSelected(), with muon = " << idx);

  // do ID cut
  bool select = ((m_sel_combined && muon.isCombinedMuon(idx)) ||
		 (m_sel_seg_tag && muon.isSegmentTaggedMuon(idx)));
  
  select = select && muon.loose(idx);

  if (!select) return false;

  //float pt = (muon->isCombinedMuon()) ? muon->pt() : muon->inDetTrackParticle()->pt(); ;
  float pt = muon.pt(idx);

  // ATH_MSG_DEBUG("Here 2");
  // if (!m_simple) {
  //   if (m_isMC && m_smearMC) {
  //     ATH_MSG_DEBUG("Here 2a");
  //     int seed = int(1.e+5*fabs(muon->phi()));
  //     if (!seed) seed = 1;
  //     m_muonSmear.SetSeed(seed);
  //     ATH_MSG_DEBUG("Here 2b");
  //     ATH_MSG_DEBUG(" args = " << muon->muonExtrapolatedTrackParticle() << ", "
  // 		    << muon->inDetTrackParticle() << ", "
  // 		    << muon->pt() << ", "
  // 		    <<  muon->eta());
      
  //     // float charge = muon->charge();
  //     float eta = muon->eta();
  //     float ptcb = muon->pt();
  //     float ptms = muon->muonExtrapolatedTrackParticle() ? muon->muonExtrapolatedTrackParticle()->pt() : 1;
  //     float ptid = muon->inDetTrackParticle() ? muon->inDetTrackParticle()->pt() : 1;
      
  //     m_muonSmear.Event(ptms,ptid,ptcb,eta);
      
  //     if (m_muonResSyst == "") {
  // 	if (muon->isCombinedMuon()) {
  // 	  pt = m_muonSmear.pTCB();
  // 	} else {
  // 	  pt = m_muonSmear.pTID();
  // 	}
  //     } else {
  // 	float pTMS_smeared = 0.;
  // 	float pTID_smeared = 0.;
  // 	float pTCB_smeared = 0.;
	
  // 	// Valid values for "THESTRING": {"MSLOW", "MSUP", "IDLOW", "IDUP"} 
  // 	m_muonSmear.PTVar(pTMS_smeared, pTID_smeared, pTCB_smeared, m_muonResSyst);
	
  // 	if (muon->isCombinedMuon()) {
  // 	  pt = pTCB_smeared;
  // 	} else {
  // 	  pt = pTID_smeared;
  // 	}
  //     }
    

  //     // let's cosnt-cast the four-mom
  //     if (pt != 0) {
  // 	Analysis::Muon* volmu = const_cast<Analysis::Muon*>(muon);
  // 	if (!volmu) {
  // 	  ATH_MSG_ERROR("Const-cast for muon did not work");
  // 	  return false;
  // 	}
  // 	volmu->setIPt(1.0/pt);
  //     }
  //   }
  // }
    
  // select must be true before in order to get here, so can overwrite
  select = pt >m_muonPt && fabs(muon.eta(idx))<m_muonEta;

  // do iso cut
  if (m_do_iso_cut) {
    if (m_do_flat_iso_cut) {
      select = select && muon.ptcone20(idx)<m_flat_isolation_cut;
  } else {
      select = select && muon.ptcone20(idx)/pt<m_isolation_cut;
    }
  }

  if (!select) return false;

  // do track cuts
  if (muon.expectBLayerHit(idx) && muon.nBLHits(idx) == 0) return false;
  if (muon.nPixHits(idx) + muon.nPixelDeadSensors(idx) < 1) return false;
  if (muon.nSCTHits(idx) + muon.nSCTDeadSensors(idx) < 5) return false;
  if (muon.nPixHoles(idx) + muon.nSCTHoles(idx) >= 3) return false;
  const int nTRTOutliers = muon.nTRTOutliers(idx);
  const int nTRTTotal    = nTRTOutliers + muon.nTRTHits(idx);
  const float trackEta   = fabs(log(tan(muon.tracktheta(idx)/2.0)));
  if (trackEta < 1.9 && trackEta > 0.1) {
    select = (nTRTTotal > 5 &&  nTRTOutliers < 0.9 * nTRTTotal);
  } else if (nTRTTotal > 5) {
    select = nTRTOutliers < 0.9*nTRTTotal;
  }
  
  ATH_MSG_DEBUG("Muon with pt = " << pt << ", select = " << select);
  
  return select;
}

bool gmsbSelectionTool::isSelected( JetD3PDObject& jet, std::size_t idx ) const
{

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
		    photon.convFlag(idx));

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
		      photon.convFlag(idx));

  const unsigned int fudgedIsEM = myTool.isEM(4,2012);

  const unsigned int origIsEM = photon.isEM(idx);

  const unsigned int newIsEM = (origIsEM & egammaPID::AMBIGUITYRESOLVE_PHOTON) | fudgedIsEM;

  photon.isEM(idx) = newIsEM;

}

