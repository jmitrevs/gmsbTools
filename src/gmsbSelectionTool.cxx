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

#include <sstream>
#include <iomanip>
#include <iostream>

//using namespace Analysis;
//using namespace Rec;
//using namespace std;

//------------------------------------------------------------------------------
gmsbSelectionTool::gmsbSelectionTool( const std::string& type,
                                                      const std::string& name, 
                                                      const IInterface* parent )
  : AthAlgTool( type, name, parent ), m_userdatasvc("UserDataSvc", name) {
  declareInterface<gmsbSelectionTool>( this );

  declareProperty("IsAtlfastData",          m_isAtlfast=false);

  declareProperty("IsMC", m_isMC=false);
  declareProperty("SmearMC", m_smearMC = false);
  declareProperty("MCHasConstantTerm", m_MCHasConstantTerm = true);
  declareProperty("RandomSeed", m_randomSeed = 0);
  declareProperty("EgammaScaleShift", m_egammaScaleShift = EnergyRescaler::NOMINAL);
  declareProperty("MCEtconeShift", m_mcEtconeShift = 314.0);

  /** caloCluster selection */
  declareProperty("CaloClusterE", m_caloClusterE=1.0*GeV);

  /** TrackParticle Pt */
  declareProperty("TrackParticlePt", m_trackParticlePt=1.0*GeV);

  /** Electron selection */
  declareProperty("ElectronPt",       m_electronPt=20*GeV);
  declareProperty("ElectronEta",      m_electronEta=2.47);
  declareProperty("ElectronIsEM", m_electronIsEM = egammaPID::ElectronMedium_WithTrackMatch);
  declareProperty("AuthorEgammaOnly", m_authorEgammaOnly=true);
  declareProperty("DoElectronEtaWindowCut", m_doElectronEtaWindCut = true);
  declareProperty("ElectronEtaWindowMin", m_electronEtaWindMin = 1.37);
  declareProperty("ElectronEtaWindowMax", m_electronEtaWindMax = 1.52);
  declareProperty("DoElectronIsolation", m_doElectronIsolation = true);
  declareProperty("ElectronEtcone20ovEt", m_electronEtcone20ovEt=0.15);

  /** Photon selection */
  declareProperty("PhotonPt",   m_photonPt=20*GeV);
  declareProperty("PhotonEta",  m_photonEta=1.81);
  declareProperty("PhotonIsEM", m_photonIsEM = egammaPID::PhotonTightAR);
  declareProperty("DoPhotonEtaWindowCut", m_doPhotonEtaWindCut = true);
  declareProperty("PhotonEtaWindowMin", m_photonEtaWindMin = 1.37);
  declareProperty("PhotonEtaWindowMax", m_photonEtaWindMax = 1.52);
  declareProperty("DoPhotonIsolation", m_doPhotonIsolation = true);
  declareProperty("PhotonEtcone20ovEt", m_photonEtcone20ovEt=0.1);

  /** Muon selection */
  declareProperty("MuonPt", m_muonPt = 20.0*GeV);
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
  declareProperty("BJetLikelihood", m_bJetLikelihood=6.0);

}

//------------------------------------------------------------------------------
StatusCode gmsbSelectionTool::initialize() {

  ATH_MSG_DEBUG("in initialize()");

  if ( !m_userdatasvc.retrieve().isSuccess() ) {
    ATH_MSG_ERROR("Unable to retrieve pointer to UserDataSvc");
    return StatusCode::FAILURE;
  }
  
  // initialize the OQ 
  m_OQ.initialize();
  m_eRescale.useDefaultCalibConstants();
  m_eRescale.SetRandomSeed(m_randomSeed);

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

bool gmsbSelectionTool::isSelected( const Analysis::Electron * electron, int runNum ) const
{
  ATH_MSG_DEBUG("in electron isSelected(), with electron = " << electron);

  if ( !electron ) return false;

  bool select = true;

  if ( m_authorEgammaOnly ) select = select && electron->author(egammaParameters::AuthorElectron);

  if (!select) return false;

  const double absClusEta = fabs(electron->cluster()->etaBE(2));

  const double eta = (electron->trackParticle()) ? electron->trackParticle()->eta() : electron->eta();

  const double uncorrectedE = electron->cluster()->e();
  const double uncorrectedEt = uncorrectedE/cosh(eta);

  double energy = uncorrectedE;
  
  if (m_isMC) {
    if (m_smearMC) {
      energy *= m_eRescale.getSmearingCorrection(electron->cluster()->eta(),
						 uncorrectedE,
						 EnergyRescaler::NOMINAL,
						 m_MCHasConstantTerm);
      double er_up=-1,er_do=-1;
      m_eRescale.getError(electron->cluster()->eta(), energy/cosh(eta), er_up, er_do,
			  "ELECTRON");


      if (m_egammaScaleShift == EnergyRescaler::ERR_UP) {
	energy *= (1+er_up);
      } else if (m_egammaScaleShift == EnergyRescaler::ERR_DOWN) {
	energy *= (1+er_do);
      }

    } 
  } else {
    energy = m_eRescale.applyEnergyCorrection(electron->cluster()->eta(), 
					      electron->cluster()->phi(), 
					      uncorrectedE,
					      uncorrectedEt,
					      m_egammaScaleShift,
					      "ELECTRON");
  }

  double pt = energy/cosh(eta);

  // add this as user data
  if (m_userdatasvc->decorateElement(*electron, std::string("corrPt"), pt)
      != StatusCode::SUCCESS) {
    ATH_MSG_ERROR("Error in electron decoration");
    return false;
  }


  if ( m_isAtlfast ) {
    select = pt > m_electronPt && absClusEta < m_electronEta;
    return select;
  }

  select = electron->isElectron(m_electronIsEM); 

  select = select && pt > m_electronPt && absClusEta <m_electronEta;
  ATH_MSG_DEBUG("select is now " << select);

  // check if electron is in bad eta region

  if(m_doElectronEtaWindCut) {
    const bool isCrack = absClusEta > m_electronEtaWindMin && absClusEta < m_electronEtaWindMax; 
    select = select && !isCrack;
  }

  ATH_MSG_DEBUG("after crack, select is now " << select);

  // check OQ
  const bool badOQ = m_OQ.checkOQClusterElectron(runNum, electron->cluster()->eta(), electron->cluster()->phi())==3;
  select = select && !badOQ;

  if ( m_doElectronIsolation ) {
    const EMShower* egdetail = electron->detail<EMShower>();
    double isol = 1000;
    if(egdetail && pt > 0.0) {
      double etcone = egdetail->etcone20();
      if (m_isMC) {
	etcone += m_mcEtconeShift;
      }
      isol = etcone / uncorrectedEt;
    }
    select = select && isol < m_electronEtcone20ovEt;
  }

  ATH_MSG_DEBUG("after iso, select is now " << select);

  return select;
}

bool gmsbSelectionTool::isSelected( const Analysis::Photon * photon, int runNum ) const 
{
  ATH_MSG_DEBUG("in photon isSelected()");

  bool select = false;
  if ( !photon ) return select;

  const double absClusEta = fabs(photon->cluster()->etaBE(2));

  double energy;
  if (m_isMC) {
    if (m_smearMC) {
      energy = photon->e() * m_eRescale.getSmearingCorrection(photon->cluster()->eta(),
							      photon->e(),
							      EnergyRescaler::NOMINAL,
							      m_MCHasConstantTerm);

      double er_up=-1,er_do=-1;
      m_eRescale.getError(photon->cluster()->eta(), energy/cosh(photon->eta()), er_up, er_do,
			  (photon->conversion()) ? "CONVERTED_PHOTON" : "UNCONVERTED_PHOTON");


      if (m_egammaScaleShift == EnergyRescaler::ERR_UP) {
	energy *= (1+er_up);
      } else if (m_egammaScaleShift == EnergyRescaler::ERR_DOWN) {
	energy *= (1+er_do);
      }
    } else {
      energy = photon->e();
    }
  } else { 
    energy = m_eRescale.applyEnergyCorrection(photon->cluster()->eta(), 
					      photon->cluster()->phi(), 
					      photon->e(),
					      photon->et(),
					      m_egammaScaleShift,
					      (photon->conversion()) ? 
					      "CONVERTED_PHOTON" : "UNCONVERTED_PHOTON");
  }

  double pt = energy/cosh(photon->eta());

  // add this as a barcode
  if (m_userdatasvc->decorateElement(*photon, std::string("corrPt"), pt)
      != StatusCode::SUCCESS) {
    ATH_MSG_ERROR("Error in photon decoration");
    return false;
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
 
  select = pt > m_photonPt && absClusEta < m_photonEta && photon->isPhoton(m_photonIsEM);

  ATH_MSG_DEBUG("after pt, eta, and isEM cut, select = " << select);

  // check if photon is in bad eta region

  if(m_doPhotonEtaWindCut) {
    const bool isCrack = absClusEta > m_photonEtaWindMin && absClusEta < m_photonEtaWindMax; 
    select = select && !isCrack;
  }

  // check OQ
  const bool badOQ = m_OQ.checkOQClusterPhoton(runNum, photon->cluster()->eta(), photon->cluster()->phi())==3;
  select = select && !badOQ;

  if ( m_doPhotonIsolation ) {
    const EMShower* egdetail = photon->detail<EMShower>();
    double isol = 1000;
    if(egdetail && pt > 0.0) {
      double etcone = egdetail->etcone20();
      if (m_isMC) {
	etcone += m_mcEtconeShift;
      }
      isol = etcone / photon->pt();
    }
    select = select && isol < m_photonEtcone20ovEt;
  }

  ATH_MSG_DEBUG("photon select = " << select);

  return select;
}

bool gmsbSelectionTool::isSelected( const Analysis::Muon * muon ) const
{
  bool select = false;
  if ( !muon ) return select;

  if ( m_isAtlfast ) {
    select = muon->pt()>m_muonPt && fabs(muon->eta())<m_muonEta;
    return select;
  }

  select = muon->pt()>m_muonPt && fabs(muon->eta())<m_muonEta;

  // do iso cut
  if (m_do_iso_cut) {
    if (m_do_flat_iso_cut) {
      select = select && muon->parameter(MuonParameters::etcone20)<m_flat_isolation_cut;
    } else {
      select = select && muon->parameter(MuonParameters::etcone20)/muon->et()<m_isolation_cut;
    }
  }

  // do ID cut
  select = select && ((m_sel_combined && muon->isCombinedMuon()) ||
		      (m_sel_seg_tag && muon->isLowPtReconstructedMuon()));
  
  if (m_sel_combined && muon->isCombinedMuon()) {
    select = select && muon->matchChi2() < m_matchChi2Max;
    double ms_pt=0.;
    const Rec::TrackParticle* msTrkParticle = muon->muonSpectrometerTrackParticle();
    if (msTrkParticle) {
      ms_pt = msTrkParticle->pt();
    }
    if (ms_pt < m_ms_pt_limit) {
      double p_ms_ex=0;
      double p_id=0;
      const Rec::TrackParticle* msExTrkParticle = muon->muonExtrapolatedTrackParticle();
      if (msExTrkParticle) {
	p_ms_ex = msExTrkParticle->p();
      }
      const Rec::TrackParticle* idTrkParticle = muon->inDetTrackParticle();
      if (idTrkParticle) {
	p_id = idTrkParticle->p();
      }
      select = select && p_ms_ex-p_id > -0.4*p_id;
    }
  }

  // do track cuts
  select = select && muon->numberOfPixelHits() >= 1 && muon->numberOfSCTHits() >= 6;
  const int nTRTOutliers = muon->numberOfTRTOutliers();
  const int nTRTTotal = muon->numberOfTRTHits() + muon->numberOfTRTOutliers();
  const Rec::TrackParticle* idTrkParticle = muon->inDetTrackParticle();
  const double mu_eta = (idTrkParticle) ? idTrkParticle->eta() : muon->eta();
  if (mu_eta < 1.9) {
    select = select && nTRTTotal > 5;
  }
  if (nTRTTotal > 5) {
    select = select && nTRTOutliers < 0.9*nTRTTotal;
  }

  return select;
}

bool gmsbSelectionTool::isSelected( const Jet* jet ) const
{
  bool select = false;
  if ( !jet ) return select;

  double emJESfactor = jet->getMoment("EMJES");
  if (emJESfactor == 0) {
    emJESfactor = m_jetEMJESfixer.fixAntiKt4H1Topo(jet->pt(), jet->eta());
  }
  Jet::hlv_t jet4MomJES =  emJESfactor * jet->hlv(P4SignalState::JETEMSCALE);

  select = jet4MomJES.perp() > m_jetPt && fabs(jet4MomJES.eta()) < m_jetEta;

  return select;
}

bool gmsbSelectionTool::isSelected( const Rec::TrackParticle * trackParticle ) const
{
  bool select = false;
  if ( !trackParticle ) return select;
  select = trackParticle->pt() > m_trackParticlePt; 

  if ( m_isAtlfast ) return select;

  return select;
}

bool gmsbSelectionTool::isSelected( const CaloCluster* caloCluster ) const
{
  bool select = false;
  if ( !caloCluster ) return select;
  select = caloCluster->e() > m_caloClusterE;

  if ( m_isAtlfast ) return select;

  return select;
}

bool gmsbSelectionTool::isSelected( const Analysis::TauJet * tauJet ) const {

  bool select = false;
  if ( !tauJet ) return select;

  int numTrack = tauJet->numTrack();
  select = tauJet->pt()>m_tauJetPt &&
           fabs(tauJet->eta())<m_tauJetEta &&
           fabs(tauJet->charge())==1.0 &&
           (numTrack==1 || numTrack==3);

  const Analysis::TauPID* tauId = tauJet->tauID();
  if ( tauId ) {
    select = select &&
             tauId->discriminant( TauJetParameters::Likelihood ) > m_tauJetLikelihood &&
             tauId->discriminant( TauJetParameters::TauElTauLikelihood )< m_tauElTauLikelihood;
  }

  return select;

}

bool gmsbSelectionTool::isBJet( const Jet * jet ) const {

  /** first check that it is a selected jet */
  if ( !jet ) return false;

  bool is_bjet = this->isSelected( jet );
  return ( is_bjet && jet->getFlavourTagWeight()>m_bJetLikelihood );
}

