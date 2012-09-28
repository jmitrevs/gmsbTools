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

// User Tools
#include "gmsbTools/gmsbSelectionTool.h"

#include "egammaAnalysisUtils/CaloIsoCorrection.h"

#include "MCTruthClassifier/IMCTruthClassifier.h"
#include "MCTruthClassifier/MCTruthClassifierDefs.h"

#include "PhotonAnalysisUtils/IPAUcaloIsolationTool.h"

#include <sstream>
#include <iomanip>
#include <iostream>


//------------------------------------------------------------------------------
gmsbSelectionTool::gmsbSelectionTool( const std::string& type,
				      const std::string& name, 
				      const IInterface* parent )
  : AthAlgTool( type, name, parent ), 
    m_userdatasvc("UserDataSvc", name),
    m_muonSmear("Data11","staco","q_pT","Rel17","")
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

  declareProperty("PAUcaloIsolationTool", m_PAUcaloIsolationTool);

  declareProperty("MCTruthClassifier", m_MCTruthClassifier);

  declareProperty("DoTruth", m_doTruth = false);


  /** caloCluster selection */
  declareProperty("CaloClusterE", m_caloClusterE=1.0*GeV);

  /** TrackParticle Pt */
  declareProperty("TrackParticlePt", m_trackParticlePt=1.0*GeV);

  /** Electron selection */
  declareProperty("ElectronPt",       m_electronPt=25*GeV);
  declareProperty("ElectronEta",      m_electronEta=2.47);
  declareProperty("ElectronID", m_electronID = egammaPID::ElectronIDMediumPP);
  declareProperty("AuthorEgammaOnly", m_authorEgammaOnly=true);
  declareProperty("DoElectronEtaWindowCut", m_doElectronEtaWindCut = false);
  declareProperty("ElectronEtaWindowMin", m_electronEtaWindMin = 1.37);
  declareProperty("ElectronEtaWindowMax", m_electronEtaWindMax = 1.52);
  declareProperty("DoOldElectronIsolation", m_doOldElectronIsolation = false);
  declareProperty("ElectronEtcone20ovEt", m_electronEtcone20ovEt=0.15);
  declareProperty("DoNewElectronIsolation", m_doNewElectronIsolation = false);
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
  declareProperty("DoOldPhotonIsolation", m_doOldPhotonIsolation = false);
  declareProperty("PhotonEtcone20ovEt", m_photonEtcone20ovEt=0.1);
  declareProperty("DoNewPhotonIsolation", m_doNewPhotonIsolation = false);
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

  /** TauJet selection */
  declareProperty("TauJetPt",           m_tauJetPt=20*GeV);
  declareProperty("TauJetEta",          m_tauJetEta=2.5);
  declareProperty("TauJetLikelihood",   m_tauJetLikelihood=-6.0);
  declareProperty("TauElTauLikelihood", m_tauElTauLikelihood=1000.0); // not yet set - No 23 1007

  /** Jet selection */
  declareProperty("JetPt",          m_jetPt=20*GeV);
  declareProperty("JetEta",         m_jetEta=100);
  declareProperty("rejectNegativeEnergyJets", m_rejNegEJets = true); 
  declareProperty("BJetLikelihood", m_bJetLikelihood=6.0);

}

//------------------------------------------------------------------------------
StatusCode gmsbSelectionTool::initialize() {

  ATH_MSG_DEBUG("in initialize()");

  // if ( !m_userdatasvc.retrieve().isSuccess() ) {
  //   ATH_MSG_ERROR("Unable to retrieve pointer to UserDataSvc");
  //   return StatusCode::FAILURE;
  // }
  
  if (m_isMC && m_doTruth) {
    if(m_MCTruthClassifier.retrieve().isFailure()) {
      ATH_MSG_ERROR("Failed to retrieve " << m_MCTruthClassifier);
      return StatusCode::FAILURE; // why success?
    }
    else {
      ATH_MSG_DEBUG("Retrieved MCTruthClassifier " << m_MCTruthClassifier);   
    }
  }

  if (m_doEDPhotonIsolation || m_doEDElectronIsolation) {
    if(m_PAUcaloIsolationTool.retrieve().isFailure()) {
      ATH_MSG_ERROR("Failed to retrieve " << m_PAUcaloIsolationTool);
      return StatusCode::FAILURE;
    } else {
      ATH_MSG_DEBUG("Retrieved PAUcaloIsolationTool " << m_PAUcaloIsolationTool);   
    }
  }

  m_eRescale.useDefaultCalibConstants("2011");
  // m_eRescale.SetRandomSeed(m_randomSeed);

  m_muonSmear.UseScale(1);
  m_muonSmear.UseImprovedCombine();

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
bool gmsbSelectionTool::isSelected( const Analysis::Electron * electron, 
				    unsigned int /* runNum */, 
				    unsigned int nPV) const
{
  ATH_MSG_DEBUG("in electron isSelected(), with electron = " << electron);

  if ( !electron ) return false;

  bool select = true;

  if ( m_authorEgammaOnly ) select = select && electron->author(egammaParameters::AuthorElectron);

  if (!select) return false;

  if (m_isMC && m_doTruth) {
    std::pair<MCTruthPartClassifier::ParticleType, MCTruthPartClassifier::ParticleOrigin> res =
      m_MCTruthClassifier->particleTruthClassifier(electron);
    if (res.first != MCTruthPartClassifier::IsoElectron) return false;
  }
  

  const double eta2 = electron->cluster()->etaBE(2);
  const double absClusEta = fabs(eta2);

  const double eta = (electron->trackParticle()) ? electron->trackParticle()->eta() : electron->eta();
  const double phi = (electron->trackParticle()) ? electron->trackParticle()->phi() : electron->phi();

  const double uncorrectedE = electron->cluster()->e();
  const double uncorrectedEt = uncorrectedE/cosh(eta);

  double energy = uncorrectedE;

  if (!m_simple) {
    // Energy calibration around the barrel-endcap crack region (both data/MC)
    energy *= m_eRescale.applyMCCalibrationMeV(electron->cluster()->eta(), uncorrectedEt,"ELECTRON");
    
    if (m_isMC) {
      if (m_elScaleShift) {
	energy *= m_eRescale.applyEnergyCorrectionMeV(electron->cluster()->eta(),  
						     electron->cluster()->phi(),  
						     uncorrectedE, 
						     uncorrectedEt, 
						     m_elScaleShift, 
						     "ELECTRON") / 
	  m_eRescale.applyEnergyCorrectionMeV(electron->cluster()->eta(),  
					      electron->cluster()->phi(),  
					      uncorrectedE, 
					      uncorrectedEt, 
					      eg2011::EnergyRescaler::NOMINAL, 
					      "ELECTRON"); 
      }
      if (m_smearMC) {	  

	int seed = int(1.e+5*fabs(electron->cluster()->phi()));
	if (!seed) seed = 1;
	m_eRescale.SetRandomSeed(seed);
	energy *= m_eRescale.getSmearingCorrectionMeV(electron->cluster()->eta(),
						      uncorrectedE,
						      m_elSmearShift,
						      m_MCHasConstantTerm,
						      "2011");
      } 
    } else {
      energy *= m_eRescale.applyEnergyCorrectionMeV(electron->cluster()->eta(),  
						    electron->cluster()->phi(),  
						    uncorrectedE, 
						    uncorrectedEt, 
						    m_elScaleShift, 
						    "ELECTRON") / uncorrectedE; 
    }

    ATH_MSG_DEBUG("Original electron E = " << uncorrectedE << ", corrected E = " << energy);

    // let's cosnt-cast the four-mom
    Analysis::Electron* volel = const_cast<Analysis::Electron*>(electron);
    if (!volel) {
      ATH_MSG_ERROR("Const-cast did not work");
      return false;
    }
    volel->setE(energy);
    volel->setEta(eta);
    volel->setPhi(phi);
  }


  const double pt = electron->pt();

  // if (!m_simple) {
  //   // add this as user data
  //   if (m_userdatasvc->decorateElement(*electron, std::string("corrPt"), pt)
  // 	!= StatusCode::SUCCESS) {
  //     ATH_MSG_ERROR("Error in electron decoration");
  //     return false;
  //   }
  // }

  if ( m_isAtlfast ) {
    select = pt > m_electronPt && absClusEta < m_electronEta;
    return select;
  }

  select = electron->passID(static_cast<egammaPID::egammaIDQuality>(m_electronID)); 

  select = select && pt > m_electronPt && absClusEta <m_electronEta;
  ATH_MSG_DEBUG("select is now " << select);

  // check if electron is in bad eta region

  if(m_doElectronEtaWindCut) {
    const bool isCrack = absClusEta > m_electronEtaWindMin && absClusEta < m_electronEtaWindMax; 
    select = select && !isCrack;
  }

  ATH_MSG_DEBUG("after crack, select is now " << select);

  // check OQ
  bool badOQ = electron->isgoodoq(egammaPID::BADCLUSELECTRON); // 0 == good
  // if (m_isMC) {
  //   badOQ = badOQ || 
  //     m_OQ.checkOQClusterElectron(runNum, electron->cluster()->eta(), electron->cluster()->phi())==3;
  // }
  select = select && !badOQ;

  if ( m_doOldElectronIsolation ) {
    const EMShower* egdetail = electron->detail<EMShower>();
    double isol = 10000000;
    if(egdetail && pt > 0.0) {
      double etcone = egdetail->etcone20();
      if (m_isMC) {
	if (!m_useAltIsoCorrection) {
	  etcone *= m_mcEtconeScale;
	} else {
	  if (absClusEta < 0.6) {
	    etcone += 0.007 * uncorrectedEt;
	  } else if (absClusEta < 1.37) {
	    etcone += 0.009 * uncorrectedEt;
	  } else {
	    etcone += 0.008 * uncorrectedEt;
	  }
	}
      }
      isol = etcone / uncorrectedEt;
    }
    select = select && isol < m_electronEtcone20ovEt;
  }

  if ( m_doNewElectronIsolation ) {
    const EMShower* egdetail = electron->detail<EMShower>();
    float isol = 1000000;
    if(egdetail && pt > 0.0) {

      isol = CaloIsoCorrection::GetPtNPVCorrectedIsolation(nPV,
							   energy,
							   eta2,
							   egdetail->parameter(egammaParameters::etap),
							   electron->cluster()->eta(),
							   20,
							   m_isMC,
							   egdetail->etcone20(),
							   false,
							   CaloIsoCorrection::ELECTRON);
    }
    select = select && isol < m_electronEtcone20corrected;
  }

  if ( m_doElectronTrackIsolation ) {
    const EMShower* egdetail = electron->detail<EMShower>();
    double isol = 10000000;
    if(egdetail && pt > 0.0) {
      double ptcone = egdetail->ptcone20();
      isol = ptcone / uncorrectedEt;
    }
    select = select && isol < m_electronPtcone20ovEt;
  }

  if ( m_doEDElectronIsolation ) {
    const EMShower* egdetail = electron->detail<EMShower>();
    float isol = 1000000;
    if(egdetail && pt > 0.0) {

      const int removeNHardestJets = 0;  // default value for now
      const double ED_correction_40 = m_PAUcaloIsolationTool->EtConeCorrectionJetAreas(electron, .40,
										       removeNHardestJets);

      ATH_MSG_DEBUG("ED_correction_40 = " << ED_correction_40);
      
      isol = CaloIsoCorrection::GetPtEDCorrectedIsolation(ED_correction_40, 0.0,
							  energy,
							  eta2,
							  egdetail->parameter(egammaParameters::etap),
							  electron->cluster()->eta(),
							  20,
							  m_isMC,
							  egdetail->etcone20(),
							  electron->conversion(),
							  CaloIsoCorrection::ELECTRON);

      ATH_MSG_DEBUG("ptcone20 = " << egdetail->etcone20() << ", pt_corrected = " << egdetail->etcone20_ptcorrected() << ", +ED = " << isol);

      EMShower *newDetail = const_cast<EMShower *>(egdetail);
      newDetail->set_parameter(egammaParameters::etcone20_ptcorrected,isol, true) ;

    }
    select = select && isol < m_electronEtcone20corrected;
  }


  ATH_MSG_DEBUG("after iso, select is now " << select);

  return select;
}

bool gmsbSelectionTool::isSelected( const Analysis::Photon * photon, 
				    unsigned int /* runNum */,
				    unsigned int nPV) const 
{
  ATH_MSG_DEBUG("in photon isSelected()");

  bool select = false;
  if ( !photon ) return select;

  const double eta2 = photon->cluster()->etaBE(2);
  const double absClusEta = fabs(eta2);

    
  if (m_isMC && m_doTruth) {
    std::pair<MCTruthPartClassifier::ParticleType, MCTruthPartClassifier::ParticleOrigin> res =
      m_MCTruthClassifier->particleTruthClassifier(photon);
    if (res.first != MCTruthPartClassifier::IsoPhoton) return false;
  }


  double energy = photon->e();
  if (!m_simple) {
    if (m_isMC) {
      if (m_phoScaleShift) {
	energy *= m_eRescale.applyEnergyCorrectionMeV(photon->cluster()->eta(),  
						      photon->cluster()->phi(),  
						      photon->e(), 
						      photon->et(), 
						      m_phoScaleShift, 
						      (photon->conversion()) ?  
						      "CONVERTED_PHOTON" : "UNCONVERTED_PHOTON") /
	  m_eRescale.applyEnergyCorrectionMeV(photon->cluster()->eta(),  
					      photon->cluster()->phi(),  
					      photon->e(), 
					      photon->et(), 
					      eg2011::EnergyRescaler::NOMINAL,
					      (photon->conversion()) ?  
					      "CONVERTED_PHOTON" : "UNCONVERTED_PHOTON"); 
	  
      }
      if (m_smearMC) {
	m_eRescale.SetRandomSeed(int(1.e+5*fabs(photon->cluster()->phi())));
	energy *= m_eRescale.getSmearingCorrectionMeV(photon->cluster()->eta(),
						      photon->e(),
						      m_phoSmearShift,
						      m_MCHasConstantTerm);
      }  
    } else { 
      energy = m_eRescale.applyEnergyCorrectionMeV(photon->cluster()->eta(),  
						   photon->cluster()->phi(),  
						   photon->e(), 
						   photon->et(), 
						   m_phoScaleShift, 
						   (photon->conversion()) ?  
						   "CONVERTED_PHOTON" : "UNCONVERTED_PHOTON"); 
      
    }

    ATH_MSG_DEBUG("Original photon E = " << photon->e() << ", corrected E = " << energy);

    // let's cosnt-cast the four-mom
    Analysis::Photon* volpho = const_cast<Analysis::Photon*>(photon);
    if (!volpho) {
      ATH_MSG_ERROR("Const-cast did not work");
      return false;
    }
    volpho->setE(energy);

  }

  const double pt = photon->pt();

  // double pt = energy/cosh(photon->eta());

  // if (!m_simple) {
  //   // add this as a barcode
  //   if (m_userdatasvc->decorateElement(*photon, std::string("corrPt"), pt)
  // 	!= StatusCode::SUCCESS) {
  //     ATH_MSG_ERROR("Error in photon decoration");
  //     return false;
  //   }
  // }

  // ATH_MSG_DEBUG("Original pt = " << photon->pt() << ", corrected photon pt = " << pt);

  // // for test
  // // get the user data
  // double readBackPt;
  // if (m_userdatasvc->getInMemElementDecoration(*photon, std::string("corrPt"), readBackPt)
  //     != StatusCode::SUCCESS) {
  //   ATH_MSG_ERROR("Error in geting photon decoration");
  //   return StatusCode::FAILURE;
  // }

  // ATH_MSG_DEBUG("Original Written photon pt = " << pt << ", read back = " << readBackPt); 

  if ( m_isAtlfast ) {
    select = pt >m_photonPt && absClusEta < m_photonEta;
    return select;
  }
 
  select = pt > m_photonPt && absClusEta < m_photonEta && photon->passID(static_cast<egammaPID::egammaIDQuality>(m_photonID)); 

  // also do isEM selection for special requirements (usually m_photonIsEM == 0, so this does nothing)
  select = select && photon->isPhoton(m_photonIsEM);

  ATH_MSG_DEBUG("after pt, eta (= " << eta2 << "), and isEM cut, select = " << select);

  // check if photon is in bad eta region

  if(m_doPhotonEtaWindCut) {
    const bool isCrack = absClusEta > m_photonEtaWindMin && absClusEta < m_photonEtaWindMax; 
    select = select && !isCrack;
  }

  ATH_MSG_DEBUG("after crack cut, select = " << select);

  // check OQ
  bool badOQ = photon->isgoodoq(egammaPID::BADCLUSPHOTON); // 0 == good
  // if (m_isMC) {
  //   badOQ = badOQ || 
  //     m_OQ.checkOQClusterPhoton(runNum, photon->cluster()->eta(), photon->cluster()->phi(), photon->conversion())==3;
  // }
  select = select && !badOQ;

  ATH_MSG_DEBUG("after OTX, select = " << select);

  if ( m_doOldPhotonIsolation ) {
    const EMShower* egdetail = photon->detail<EMShower>();
    double isol = 10000000;
    if(egdetail && pt > 0.0) {
      double etcone = egdetail->etcone20();
      ATH_MSG_DEBUG("etcone20 = " << etcone);
      if (m_isMC) {
	if (!m_useAltIsoCorrection) {
	  etcone *= m_mcEtconeScale;
	} else {
	  if (absClusEta < 0.6) {
	    etcone += 0.007 * photon->pt();
	  } else if (absClusEta < 1.37) {
	    etcone += 0.009 * photon->pt();
	  } else {
	    etcone += 0.008 * photon->pt();
	  }
	}
      }
      isol = etcone / photon->pt();
    }
    select = select && isol < m_photonEtcone20ovEt;
  }

  if ( m_doNewPhotonIsolation ) {
    const EMShower* egdetail = photon->detail<EMShower>();
    float isol = 1000000;
    if(egdetail && pt > 0.0) {
      // ATH_MSG_DEBUG("arguments to GetPtNPVCorrectedIsolation: " << nPV << ", " <<
      // 		    energy << ", " <<
      // 		    eta2 << ", " <<
      // 		    egdetail->parameter(egammaParameters::etap) << ", " <<
      // 		    photon->cluster()->eta() << ", " <<
      // 		    20 << ", " <<
      // 		    m_isMC << ", " <<
      // 		    egdetail->etcone20() << ", " <<
      // 		    photon->conversion() << ", " <<
      // 		    CaloIsoCorrection::PHOTON);

      isol = CaloIsoCorrection::GetPtNPVCorrectedIsolation(nPV,
							   energy,
							   eta2,
							   egdetail->parameter(egammaParameters::etap),
							   photon->cluster()->eta(),
							   20,
							   m_isMC,
							   egdetail->etcone20(),
							   photon->conversion(),
							   CaloIsoCorrection::PHOTON);
    }
    select = select && isol < m_photonEtcone20corrected;
  }

  if ( m_doEDPhotonIsolation ) {
    const EMShower* egdetail = photon->detail<EMShower>();
    float isol = 1000000;
    if(egdetail && pt > 0.0) {

      const int removeNHardestJets = 0;  // default value for now
      const double ED_correction_40 = m_PAUcaloIsolationTool->EtConeCorrectionJetAreas(photon, .40,
										       removeNHardestJets);

      ATH_MSG_DEBUG("ED_correction_40 = " << ED_correction_40);
      
      isol = CaloIsoCorrection::GetPtEDCorrectedIsolation(ED_correction_40, 0.0,
							  energy,
							  eta2,
							  egdetail->parameter(egammaParameters::etap),
							  photon->cluster()->eta(),
							  20,
							  m_isMC,
							  egdetail->etcone20(),
							  photon->conversion(),
							  CaloIsoCorrection::PHOTON);

      ATH_MSG_DEBUG("ptcone20 = " << egdetail->etcone20() << ", pt_corrected = " << egdetail->etcone20_ptcorrected() << ", +ED = " << isol);

      EMShower *newDetail = const_cast<EMShower *>(egdetail);
      newDetail->set_parameter(egammaParameters::etcone20_ptcorrected,isol, true) ;

    }
    select = select && isol < m_photonEtcone20corrected;
  }

  ATH_MSG_DEBUG("photon select = " << select);

  return select;
}

bool gmsbSelectionTool::isSelected( const Analysis::Muon * muon ) const
{
  ATH_MSG_DEBUG("in muon isSelected(), with muon = " << muon);

  if ( !muon ) return false;

  if ( m_isAtlfast ) {
    return (muon->pt()>m_muonPt && fabs(muon->eta())<m_muonEta);
  }

  ATH_MSG_DEBUG("Here");

  // do ID cut
  bool select = ((m_sel_combined && muon->isCombinedMuon()) ||
		 (m_sel_seg_tag && muon->isSegmentTaggedMuon()));
  
  select = select && muon->isLoose();

  if (!select) return false;

  //double pt = (muon->isCombinedMuon()) ? muon->pt() : muon->inDetTrackParticle()->pt(); ;
  double pt = muon->pt();

  // ATH_MSG_DEBUG("Here 2");
  if (!m_simple) {
    if (m_isMC && m_smearMC) {
      ATH_MSG_DEBUG("Here 2a");
      int seed = int(1.e+5*fabs(muon->phi()));
      if (!seed) seed = 1;
      m_muonSmear.SetSeed(seed);
      ATH_MSG_DEBUG("Here 2b");
      ATH_MSG_DEBUG(" args = " << muon->muonExtrapolatedTrackParticle() << ", "
		    << muon->inDetTrackParticle() << ", "
		    << muon->pt() << ", "
		    <<  muon->eta());
      
      // double charge = muon->charge();
      double eta = muon->eta();
      double ptcb = muon->pt();
      double ptms = muon->muonExtrapolatedTrackParticle() ? muon->muonExtrapolatedTrackParticle()->pt() : 1;
      double ptid = muon->inDetTrackParticle() ? muon->inDetTrackParticle()->pt() : 1;
      
      m_muonSmear.Event(ptms,ptid,ptcb,eta);
      
      if (m_muonResSyst == "") {
	if (muon->isCombinedMuon()) {
	  pt = m_muonSmear.pTCB();
	} else {
	  pt = m_muonSmear.pTID();
	}
      } else {
	double pTMS_smeared = 0.;
	double pTID_smeared = 0.;
	double pTCB_smeared = 0.;
	
	// Valid values for "THESTRING": {"MSLOW", "MSUP", "IDLOW", "IDUP"} 
	m_muonSmear.PTVar(pTMS_smeared, pTID_smeared, pTCB_smeared, m_muonResSyst);
	
	if (muon->isCombinedMuon()) {
	  pt = pTCB_smeared;
	} else {
	  pt = pTID_smeared;
	}
      }
    

      // let's cosnt-cast the four-mom
      if (pt != 0) {
	Analysis::Muon* volmu = const_cast<Analysis::Muon*>(muon);
	if (!volmu) {
	  ATH_MSG_ERROR("Const-cast for muon did not work");
	  return false;
	}
	volmu->setIPt(1.0/pt);
      }
    }
  }
    
  // select must be true before in order to get here, so can overwrite
  select = pt >m_muonPt && fabs(muon->eta())<m_muonEta;

  // do iso cut
  if (m_do_iso_cut) {
    if (m_do_flat_iso_cut) {
      select = select && muon->parameter(MuonParameters::ptcone20)<m_flat_isolation_cut;
    } else {
      select = select && muon->parameter(MuonParameters::ptcone20)/muon->et()<m_isolation_cut;
    }
  }

  if (!select) return false;

  // do track cuts
  const Rec::TrackParticle* id_trk_part = muon->inDetTrackParticle();
  if(!id_trk_part) return false;
  const Trk::TrackSummary* id_trk_sum = id_trk_part->trackSummary();
  if(!id_trk_sum) return false;
  if(id_trk_sum->get(Trk::expectBLayerHit) && id_trk_sum->get(Trk::numberOfBLayerHits) == 0) return false;
  if (id_trk_sum->get(Trk::numberOfPixelHits) + id_trk_sum->get(Trk::numberOfPixelDeadSensors) <= 1) return false;
  if (id_trk_sum->get(Trk::numberOfSCTHits) + id_trk_sum->get(Trk::numberOfSCTDeadSensors) < 6) return false;
  if (id_trk_sum->get(Trk::numberOfPixelHoles) + id_trk_sum->get(Trk::numberOfSCTHoles) >= 3) return false;
  const int nTRTOutliers = id_trk_sum->get(Trk::numberOfTRTOutliers);
  const int nTRTTotal    = nTRTOutliers + id_trk_sum->get(Trk::numberOfTRTHits);
  const float trackEta   = id_trk_part->eta();
  if (fabs(trackEta) < 1.9) {
    select = (nTRTTotal > 5 &&  nTRTOutliers < 0.9 * nTRTTotal);
  } else if (nTRTTotal > 5) {
    select = nTRTOutliers < 0.9*nTRTTotal;
  }
  
  ATH_MSG_DEBUG("Muon with pt = " << pt << ", select = " << select);
  
  return select;
}

bool gmsbSelectionTool::isSelected( const Jet* jet ) const
{
  if ( !jet ) return false;

  if (m_rejNegEJets && jet->e() < 0) {
    return false;
  }

  return (jet->pt() > m_jetPt && fabs(jet->eta()) < m_jetEta);

}

bool gmsbSelectionTool::isSelected( const Rec::TrackParticle * trackParticle ) const
{
  bool select = false;
  if ( !trackParticle ) return select;
  select = trackParticle->pt() > m_trackParticlePt; 

  return select;
}

bool gmsbSelectionTool::isSelected( const CaloCluster* caloCluster ) const
{
  bool select = false;
  if ( !caloCluster ) return select;
  select = caloCluster->e() > m_caloClusterE;

  return select;
}

bool gmsbSelectionTool::isSelected( const Analysis::TauJet * /* tauJet */ ) const {

  bool select = false;
  // if ( !tauJet ) return select;

  // int numTrack = tauJet->numTrack();
  // select = tauJet->pt()>m_tauJetPt &&
  //          fabs(tauJet->eta())<m_tauJetEta &&
  //          fabs(tauJet->charge())==1.0 &&
  //          (numTrack==1 || numTrack==3);

  // const Analysis::TauPID* tauId = tauJet->tauID();
  // if ( tauId ) {
  //   select = select &&
  //            tauId->discriminant( TauJetParameters::Likelihood ) > m_tauJetLikelihood &&
  //            tauId->discriminant( TauJetParameters::TauElTauLikelihood )< m_tauElTauLikelihood;
  // }

  return select;

}

bool gmsbSelectionTool::isBJet( const Jet * jet ) const {

  /** first check that it is a selected jet */
  if ( !jet ) return false;

  bool is_bjet = this->isSelected( jet );
  return ( is_bjet && jet->getFlavourTagWeight()>m_bJetLikelihood );
}

