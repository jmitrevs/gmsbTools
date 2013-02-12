/*****************************************************************************
Name    : gmsbSelectionTool.cxx
Package : offline/PhysicsAnalysis/AnalysisCommon/gmsbTools
Author  : Ketevi A. Assamagan
Created : November 2007
Purpose : User Analysis Selections - see gmsbSelectionTool.h for details
*****************************************************************************/

//#include "GaudiKernel/GaudiException.h"
//#include "GaudiKernel/Property.h"

#include "egammaEvent/EMShower.h"

// Accessing data:
#include "CLHEP/Units/PhysicalConstants.h"
#include "gmsbD3PDObjects/ElectronD3PDObject.h"
#include "gmsbD3PDObjects/MuonD3PDObject.h"
#include "gmsbD3PDObjects/JetD3PDObject.h"
#include "gmsbD3PDObjects/PhotonD3PDObject.h"

// User Tools
#include "gmsbTools/gmsbSelectionTool.h"

//#include "egammaAnalysisUtils/CaloIsoCorrection.h"

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

  declareProperty("IsMC", m_isMC=false);
  declareProperty("SmearMC", m_smearMC = true);
  declareProperty("MCHasConstantTerm", m_MCHasConstantTerm = false);
  //  declareProperty("RandomSeed", m_randomSeed = 0); // use SUSY prescription
  declareProperty("ElScaleShift", m_elScaleShift = eg2011::EnergyRescaler::NOMINAL);
  declareProperty("ElSmearShift", m_elSmearShift = eg2011::EnergyRescaler::NOMINAL);
  declareProperty("PhoScaleShift", m_phoScaleShift = eg2011::EnergyRescaler::NOMINAL);
  declareProperty("PhoSmearShift", m_phoSmearShift = eg2011::EnergyRescaler::NOMINAL);
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
  declareProperty("DoElectronTrackIsolation", m_doElectronTrackIsolation = false);
  declareProperty("ElectronPtcone20ovEt", m_electronPtcone20ovEt=0.1);
  declareProperty("DoEDElectronIsolation", m_doEDElectronIsolation = false);
  declareProperty("Simple", m_simple=false); // don't smear or decorate object
                                             // (useful for selecting already selected)

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

  /** Muon selection */
  declareProperty("MuonPt", m_muonPt = 10.0*GeV);
  declareProperty("MuonEta", m_muonEta = 2.4);
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
bool gmsbSelectionTool::isSelected( ElectronD3PDObject& electron, std::size_t idx ) const
{
  ATH_MSG_DEBUG("in electron isSelected(), with electron = " << idx);

  bool select = true;

  if ( m_authorEgammaOnly ) select = select && electron.author(idx, egammaParameters::AuthorElectron);

  if (!select) return false;

  if (m_isMC && m_doTruth) {
    if (electron.type(idx) != MCTruthPartClassifier::IsoElectron) return false;
  }
  

  const float eta2 = electron.etas2(idx);
  const float absClusEta = fabs(eta2);

  const float eta = electron.tracketa(idx);
  const float phi = electron.trackphi(idx);

  const float uncorrectedE = electron.cl_E(idx);
  const float uncorrectedEt = uncorrectedE/cosh(eta);

  float energy = uncorrectedE;

  if (!m_simple) {
    // Energy calibration around the barrel-endcap crack region (both data/MC)
    // energy *= m_eRescale.applyMCCalibrationMeV(electron->cluster()->eta(), uncorrectedEt,"ELECTRON");
    
    // if (m_isMC) {
    //   if (m_elScaleShift) {
    // 	energy *= m_eRescale.applyEnergyCorrectionMeV(electron->cluster()->eta(),  
    // 						     electron->cluster()->phi(),  
    // 						     uncorrectedE, 
    // 						     uncorrectedEt, 
    // 						     m_elScaleShift, 
    // 						     "ELECTRON") / 
    // 	  m_eRescale.applyEnergyCorrectionMeV(electron->cluster()->eta(),  
    // 					      electron->cluster()->phi(),  
    // 					      uncorrectedE, 
    // 					      uncorrectedEt, 
    // 					      eg2011::EnergyRescaler::NOMINAL, 
    // 					      "ELECTRON"); 
    //   }
    //   if (m_smearMC) {	  

    // 	int seed = int(1.e+5*fabs(electron->cluster()->phi()));
    // 	if (!seed) seed = 1;
    // 	m_eRescale.SetRandomSeed(seed);
    // 	energy *= m_eRescale.getSmearingCorrectionMeV(electron->cluster()->eta(),
    // 						      uncorrectedE,
    // 						      m_elSmearShift,
    // 						      m_MCHasConstantTerm,
    // 						      "2011");
    //   } 
    // } else {
    //   energy *= m_eRescale.applyEnergyCorrectionMeV(electron->cluster()->eta(),  
    // 						    electron->cluster()->phi(),  
    // 						    uncorrectedE, 
    // 						    uncorrectedEt, 
    // 						    m_elScaleShift, 
    // 						    "ELECTRON") / uncorrectedE; 
    // }

    ATH_MSG_DEBUG("Original electron E = " << uncorrectedE << ", corrected E = " << energy);

    // let's cosnt-cast the four-mom

    electron.E(idx) = energy;
    electron.eta(idx) = eta;
    electron.phi(idx) = phi;
    electron.pt(idx) = energy/cosh(eta); // massless
  }


  const float pt = electron.pt(idx);

  if ( m_isAtlfast ) {
    select = pt > m_electronPt && absClusEta < m_electronEta;
    return select;
  }

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

  if ( m_doElectronTrackIsolation ) {
    const float isol = electron.ptcone20_zpv05(idx) / pt;
    select = select && isol < m_electronPtcone20ovEt;
  }

  if ( m_doEDElectronIsolation ) {
    const float isol = electron.topoEtcone20_corrected(idx);
    select = select && isol < m_electronEtcone20corrected;
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


  const float eta2 = photon.eta2(idx);
  const float absClusEta = fabs(eta2);


  float energy = photon.E(idx);
  if (!m_simple) {
    // if (m_isMC) {
    //   if (m_phoScaleShift) {
    // 	energy *= m_eRescale.applyEnergyCorrectionMeV(photon->cluster()->eta(),  
    // 						      photon->cluster()->phi(),  
    // 						      photon->e(), 
    // 						      photon->et(), 
    // 						      m_phoScaleShift, 
    // 						      (photon->conversion()) ?  
    // 						      "CONVERTED_PHOTON" : "UNCONVERTED_PHOTON") /
    // 	  m_eRescale.applyEnergyCorrectionMeV(photon->cluster()->eta(),  
    // 					      photon->cluster()->phi(),  
    // 					      photon->e(), 
    // 					      photon->et(), 
    // 					      eg2011::EnergyRescaler::NOMINAL,
    // 					      (photon->conversion()) ?  
    // 					      "CONVERTED_PHOTON" : "UNCONVERTED_PHOTON"); 
	  
    //   }
    //   if (m_smearMC) {
    // 	m_eRescale.SetRandomSeed(int(1.e+5*fabs(photon->cluster()->phi())));
    // 	energy *= m_eRescale.getSmearingCorrectionMeV(photon->cluster()->eta(),
    // 						      photon->e(),
    // 						      m_phoSmearShift,
    // 						      m_MCHasConstantTerm);
    //   }  
    // } else { 
    //   energy = m_eRescale.applyEnergyCorrectionMeV(photon->cluster()->eta(),  
    // 						   photon->cluster()->phi(),  
    // 						   photon->e(), 
    // 						   photon->et(), 
    // 						   m_phoScaleShift, 
    // 						   (photon->conversion()) ?  
    // 						   "CONVERTED_PHOTON" : "UNCONVERTED_PHOTON"); 
      
    // }

    ATH_MSG_DEBUG("Original photon E = " << photon.E(idx) << ", corrected E = " << energy);

    const float scale = energy / photon.E(idx);

    photon.E(idx) = energy;
    photon.pt(idx) *= scale; 

  }

  const float pt = photon.pt(idx);

  if ( m_isAtlfast ) {
    select = pt >m_photonPt && absClusEta < m_photonEta;
    return select;
  }
 
  select = pt > m_photonPt && absClusEta < m_photonEta && photon.passID(idx, static_cast<egammaPID::egammaIDQuality>(m_photonID)); 

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

bool gmsbSelectionTool::isSelected( MuonD3PDObject& muon, std::size_t idx ) const
{
  ATH_MSG_DEBUG("in muon isSelected(), with muon = " << idx);

  if ( m_isAtlfast ) {
    return (muon.pt(idx) > m_muonPt && fabs(muon.eta(idx)) < m_muonEta);
  }


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
  if (muon.nPixHits(idx) + muon.nPixelDeadSensors(idx) <= 1) return false;
  if (muon.nSCTHits(idx) + muon.nSCTDeadSensors(idx) < 6) return false;
  if (muon.nPixHoles(idx) + muon.nSCTHoles(idx) >= 3) return false;
  const int nTRTOutliers = muon.nTRTOutliers(idx);
  const int nTRTTotal    = nTRTOutliers + muon.nTRTHits(idx);
  const float trackEta   = - log(tan(muon.tracktheta(idx)/2.0));;
  if (fabs(trackEta) < 1.9) {
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

