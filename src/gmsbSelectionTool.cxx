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

#include <sstream>
#include <iomanip>
#include <iostream>


//------------------------------------------------------------------------------
gmsbSelectionTool::gmsbSelectionTool( const std::string& type,
				      const std::string& name, 
				      const IInterface* parent )
  : AthAlgTool( type, name, parent ), 
    m_userdatasvc("UserDataSvc", name)
    //    m_muonSmear("staco")
{
  declareInterface<gmsbSelectionTool>( this );

  declareProperty("Atlfast",          m_isAtlfast=false);

  declareProperty("IsMC", m_isMC=false);
  declareProperty("SmearMC", m_smearMC = false);
  declareProperty("MCHasConstantTerm", m_MCHasConstantTerm = true);
  //  declareProperty("RandomSeed", m_randomSeed = 0); // use SUSY prescription
  declareProperty("EgammaScaleShift", m_egammaScaleShift = eg2011::EnergyRescaler::NOMINAL);
  declareProperty("MCEtconeScale", m_mcEtconeScale = 1.5);
  declareProperty("MCUseAltIsoCorrection", m_useAltIsoCorrection = false);

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
  declareProperty("DoNewElectronIsolation", m_doNewElectronIsolation = true);
  declareProperty("ElectronEtcone20corrected", m_electronEtcone20corrected=5*GeV);
  declareProperty("DoElectronTrackIsolation", m_doElectronTrackIsolation = false);
  declareProperty("ElectronPtcone20ovEt", m_electronPtcone20ovEt=0.1);
  declareProperty("Simple", m_simple=false); // don't smear or decorate object
                                             // (useful for selecting already selected)

  /** Photon selection */
  declareProperty("PhotonPt",   m_photonPt=25*GeV);
  declareProperty("PhotonEta",  m_photonEta=2.37);
  declareProperty("PhotonID", m_photonID = egammaPID::PhotonIDTightAR);
  declareProperty("DoPhotonEtaWindowCut", m_doPhotonEtaWindCut = false);
  declareProperty("PhotonEtaWindowMin", m_photonEtaWindMin = 1.37);
  declareProperty("PhotonEtaWindowMax", m_photonEtaWindMax = 1.52);
  declareProperty("DoOldPhotonIsolation", m_doOldPhotonIsolation = false);
  declareProperty("PhotonEtcone20ovEt", m_photonEtcone20ovEt=0.1);
  declareProperty("DoNewPhotonIsolation", m_doNewPhotonIsolation = true);
  declareProperty("PhotonEtcone20corrected", m_photonEtcone20corrected=5*GeV);

  /** Muon selection */
  declareProperty("MuonPt", m_muonPt = 10.0*GeV);
  declareProperty("MuonEta", m_muonEta = 2.4);
  declareProperty("MuonMatchChi2Max", m_matchChi2Max = 150.0);
  declareProperty("SelectCombined", m_sel_combined=true);
  declareProperty("SelectSegmentTag", m_sel_seg_tag=true);
  declareProperty("DoMuonIsoCut",m_do_iso_cut = false);
  declareProperty("DoMuonFlatIsoCut",m_do_flat_iso_cut = true);
  declareProperty("MuonFlatIsoCut",m_flat_isolation_cut = 10.*GeV);
  declareProperty("MuonIsoCut",m_isolation_cut = 0.1);
  declareProperty("MuonSpecPtLimit", m_ms_pt_limit = 50.*GeV);
  declareProperty("MuonSegMomDiffLimit", m_ms_p_diff_limit = -0.4);

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

  if ( !m_userdatasvc.retrieve().isSuccess() ) {
    ATH_MSG_ERROR("Unable to retrieve pointer to UserDataSvc");
    return StatusCode::FAILURE;
  }
  
  if (m_isMC && m_doTruth) {
    if(m_MCTruthClassifier.retrieve().isFailure()) {
      ATH_MSG_ERROR("Failed to retrieve " << m_MCTruthClassifier);
      return StatusCode::FAILURE; // why success?
    }
    else {
      ATH_MSG_DEBUG("Retrieved MCTruthClassifier " << m_MCTruthClassifier);   
    }
  }

  m_eRescale.useDefaultCalibConstants();
  // m_eRescale.SetRandomSeed(m_randomSeed);

  // m_muonSmear.UseScale(1);

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
				    unsigned int nPV,
				    double pt) const
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

  const double uncorrectedE = electron->cluster()->e();
  const double uncorrectedEt = uncorrectedE/cosh(eta);

  double energy = uncorrectedE;

  if (!m_simple) {
    if (m_isMC) {
      if (m_smearMC) {
	m_eRescale.SetRandomSeed(int(1.e+5*fabs(electron->cluster()->phi())));
	energy *= m_eRescale.getSmearingCorrectionMeV(electron->cluster()->eta(),
						      uncorrectedE,
						      eg2011::EnergyRescaler::NOMINAL,
						      m_MCHasConstantTerm);
	double er_up=-1,er_do=-1;
	m_eRescale.getErrorMeV(electron->cluster()->eta(), energy/cosh(eta), er_up, er_do,
			       "ELECTRON");
	
	
	if (m_egammaScaleShift == eg2011::EnergyRescaler::ERR_UP) {
	  energy *= (1+er_up);
	} else if (m_egammaScaleShift == eg2011::EnergyRescaler::ERR_DOWN) {
	  energy *= (1+er_do);
	}

      } 
    } else {
      energy = m_eRescale.applyEnergyCorrectionMeV(electron->cluster()->eta(),  
						   electron->cluster()->phi(),  
						   uncorrectedE, 
						   uncorrectedEt, 
						   m_egammaScaleShift, 
						   "ELECTRON"); 
    }
  }

  if (pt == -1.0) {
    pt = energy/cosh(eta);
  }

  if (!m_simple) {
    // add this as user data
    if (m_userdatasvc->decorateElement(*electron, std::string("corrPt"), pt)
	!= StatusCode::SUCCESS) {
      ATH_MSG_ERROR("Error in electron decoration");
      return false;
    }
  }

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


  ATH_MSG_DEBUG("after iso, select is now " << select);

  return select;
}

bool gmsbSelectionTool::isSelected( const Analysis::Photon * photon, 
				    unsigned int /* runNum */,
				    unsigned int nPV,
				    double pt) const 
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
      if (m_smearMC) {
	m_eRescale.SetRandomSeed(int(1.e+5*fabs(photon->cluster()->phi())));
	energy = photon->e() * m_eRescale.getSmearingCorrectionMeV(photon->cluster()->eta(),
								   photon->e(),
								   eg2011::EnergyRescaler::NOMINAL,
								   m_MCHasConstantTerm);
	
	double er_up=-1,er_do=-1;
	m_eRescale.getErrorMeV(photon->cluster()->eta(), energy/cosh(photon->eta()), er_up, er_do,
			       (photon->conversion()) ? "CONVERTED_PHOTON" : "UNCONVERTED_PHOTON");
	
	
	if (m_egammaScaleShift == eg2011::EnergyRescaler::ERR_UP) {
	  energy *= (1+er_up);
	} else if (m_egammaScaleShift == eg2011::EnergyRescaler::ERR_DOWN) {
	  energy *= (1+er_do);
	}
      }  
    } else { 
      energy = m_eRescale.applyEnergyCorrectionMeV(photon->cluster()->eta(),  
						   photon->cluster()->phi(),  
						   photon->e(), 
						   photon->et(), 
						   m_egammaScaleShift, 
						   (photon->conversion()) ?  
						   "CONVERTED_PHOTON" : "UNCONVERTED_PHOTON"); 
      
    }
  }

  if (pt == -1.0) {
    pt = energy/cosh(photon->eta());
  }

  if (!m_simple) {
    // add this as a barcode
    if (m_userdatasvc->decorateElement(*photon, std::string("corrPt"), pt)
	!= StatusCode::SUCCESS) {
      ATH_MSG_ERROR("Error in photon decoration");
      return false;
    }
  }

  ATH_MSG_DEBUG("Original pt = " << photon->pt() << ", corrected photon pt = " << pt);

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

  // ATH_MSG_DEBUG("Here 1");

  double pt = (muon->isCombinedMuon()) ? muon->pt() : muon->inDetTrackParticle()->pt(); ;

  // ATH_MSG_DEBUG("Here 2");
  // if (m_isMC && m_smearMC) {
  //   ATH_MSG_DEBUG("Here 2a");
  //   m_muonSmear.SetSeed(int(1.e+5*fabs(muon->phi())));
  //   ATH_MSG_DEBUG("Here 2b");
  //   ATH_MSG_DEBUG(" args = " << muon->muonExtrapolatedTrackParticle() << ", "
  // 		  << muon->inDetTrackParticle() << ", "
  // 		  << muon->pt() << ", "
  // 		  <<  muon->eta());

  //   if (muon->isCombinedMuon()) {
  //     m_muonSmear.Event(muon->muonExtrapolatedTrackParticle()->pt(),
  // 			muon->inDetTrackParticle()->pt(),
  // 			muon->pt(),
  // 			muon->eta());
  //     pt = m_muonSmear.pTCB();
  //   } else {
  //     m_muonSmear.Event(muon->inDetTrackParticle()->pt(),
  // 			muon->eta(), "ID");
  //     pt = m_muonSmear.pTID();
  //   }
  // }
    
  // ATH_MSG_DEBUG("Here3");

  // select must be true before in order to get here, so can overwrite
  select = pt >m_muonPt && fabs(muon->eta())<m_muonEta;

  // do iso cut
  if (m_do_iso_cut) {
    if (m_do_flat_iso_cut) {
      select = select && muon->parameter(MuonParameters::etcone20)<m_flat_isolation_cut;
    } else {
      select = select && muon->parameter(MuonParameters::etcone20)/muon->et()<m_isolation_cut;
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
  //ATH_MSG_DEBUG("Here6, select = " << select);
  
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

