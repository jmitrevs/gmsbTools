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
  : AthAlgTool( type, name, parent ) {
  declareInterface<gmsbSelectionTool>( this );

  declareProperty("IsAtlfastData",          m_isAtlfast=false);

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
  declareProperty("MuonPt",                 m_muonPt=3.0*GeV);
  declareProperty("MuonEta",                m_muonEta=2.7);
  declareProperty("AuthorHighPtOnly",       m_authorHighPtOnly=false);
  declareProperty("DoMuonIsolation",        m_doMuonIsolation=true);
  declareProperty("MuonIsolationConeIndex", m_muonIsolationConeIndex=1);
  declareProperty("MuonIsolationEt",        m_muonIsolationEt=10*GeV);
  declareProperty("UseMatchChi2",           m_useMatchChi2=false);
  declareProperty("MuonMatchChi2",          m_muonMatchChi2=100);
  declareProperty("NormalizedMuonIsolationEt",m_normMuonIsolEt=0.2);

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

bool gmsbSelectionTool::isSelected( const Analysis::Electron * electron ) const
{
  ATH_MSG_DEBUG("in electron isSelected(), with electron = " << electron);

  if ( !electron ) return false;

  bool select = true;

  if ( m_authorEgammaOnly ) select = select && electron->author(egammaParameters::AuthorElectron);

  if (!select) return false;

  const double absClusEta = fabs(electron->cluster()->etaBE(2));
  const double pt = electron->pt();

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


  if ( m_doElectronIsolation ) {
    const EMShower* egdetail = electron->detail<EMShower>();
    double isol = 1000;
    if(egdetail && pt > 0.0) {
      isol = egdetail->etcone20() / pt;
    }
    select = select && isol < m_electronEtcone20ovEt;
  }

  ATH_MSG_DEBUG("after iso, select is now " << select);

  return select;
}

bool gmsbSelectionTool::isSelected( const Analysis::Photon * photon ) const 
{
  ATH_MSG_DEBUG("in photon isSelected()");

  bool select = false;
  if ( !photon ) return select;

  const double absClusEta = fabs(photon->cluster()->etaBE(2));
  const double pt = photon->pt();

  if ( m_isAtlfast ) {
    select = pt >m_photonPt && absClusEta < m_photonEta;
    return select;
  }
 
  select = pt > m_photonPt && absClusEta < m_photonEta && photon->isPhoton(m_photonIsEM);

  // check if photon is in bad eta region

  if(m_doPhotonEtaWindCut) {
    const bool isCrack = absClusEta > m_photonEtaWindMin && absClusEta < m_photonEtaWindMax; 
    select = select && !isCrack;
  }


  if ( m_doPhotonIsolation ) {
    const EMShower* egdetail = photon->detail<EMShower>();
    double isol = 1000;
    if(egdetail && pt > 0.0) {
      isol = egdetail->etcone20() / pt;
    }
    select = select && isol < m_photonEtcone20ovEt;
  }

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
  if ( m_authorHighPtOnly ) select = select && muon->isHighPt();
  if ( m_useMatchChi2 && muon->isCombinedMuon() ) select = select && muon->matchChi2() < m_muonMatchChi2;
  if ( m_doMuonIsolation ) {
    double etIsol = muon->parameter( static_cast<MuonParameters::ParamDef>(m_muonIsolationConeIndex) );
    select = select && ( etIsol < m_muonIsolationEt );
    if ( muon->pt() ) select = select && ( (etIsol/muon->pt()) < m_normMuonIsolEt );
 
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

