/*****************************************************************************
Name    : gmsbPreparationTool.cxx
Package : offline/PhysicsAnalysis/AnalysisCommon/gmsbTools
Author  : Ketevi A. Assamagan
Created : November 2007
Purpose : User Analysis Overlap Removal - see gmsbOverlapRemovalTool.h for details
*****************************************************************************/

#include "GaudiKernel/GaudiException.h"
#include "GaudiKernel/Property.h"

// Accessing data:
#include "CLHEP/Units/PhysicalConstants.h"

// User Tools
#include "gmsbTools/gmsbOverlapRemovalTool.h"

#include "muonEvent/Muon.h"
#include "egammaEvent/Electron.h"
#include "egammaEvent/Photon.h"
#include "CaloEvent/CaloCluster.h"
#include "Particle/TrackParticle.h"
#include "tauEvent/TauJet.h"
#include "JetEvent/Jet.h"

#include <sstream>
#include <iomanip>
#include <iostream>

using namespace Analysis;
using namespace Rec;


//------------------------------------------------------------------------------
gmsbOverlapRemovalTool::gmsbOverlapRemovalTool( const std::string& type,
                                                                const std::string& name, 
                                                                const IInterface* parent )
  : AthAlgTool( type, name, parent ),
    m_userSelectionTool ( "gmsbSelectionTool"),
    m_userOverlapCheckingTool ( "gmsbOverlapCheckingTool") {

  declareInterface<gmsbOverlapRemovalTool>( this );

  declareProperty("UserSelectionTool",       m_userSelectionTool);
  declareProperty("UserOverlapCheckingTool", m_userOverlapCheckingTool);
  declareProperty("InputContainerKeys",      m_inputContainerKeys);
  declareProperty("IsAtlfastData",           m_isAtlfast=false);

  declareProperty("OuputObjectKey",         m_outputObjectKey        = "FinalStateObjects");
  declareProperty("OutputLeptonKey",        m_outputLeptonKey        = "FinalStateLeptons");
  declareProperty("OutputPhotonKey",        m_outputPhotonKey        = "FinalStatePhotons");
  declareProperty("OutputElectronKey",      m_outputElectronKey      = "FinalStateElectrons");
  declareProperty("OutputMuonKey",          m_outputMuonKey          = "FinalStateMuons");
  declareProperty("OutputTauJetKey",        m_outputTauJetKey        = "FinalStateTauJets");
  declareProperty("OutputCalloClusterKey",  m_outputCaloClusterKey   = "FinalStateCaloClusters");
  declareProperty("OutputTrackParticleKey", m_outputTrackParticleKey = "FinalStateTrackParticles");
  declareProperty("OutputJetKey",           m_outputJetKey           = "FinalStateJets");
  declareProperty("OutputBJetKey",          m_outputBJetKey          = "FinalStateBJets");
  declareProperty("OutputLightJetKey",      m_outputLightJetKey      = "FinalStateLightJets");
  declareProperty("RemoveOverlapInSameContainer", m_removeOverlapInSameContainer=true);

  /** initialize counters */
  m_numElectrons      = std::make_pair(0,0);
  m_numPhotons        = std::make_pair(0,0);
  m_numMuons          = std::make_pair(0,0);
  m_numTauJets        = std::make_pair(0,0);
  m_numJets           = std::make_pair(0,0);
  m_numBJets          = std::make_pair(0,0);
  m_numLightJets      = std::make_pair(0,0);
  m_numTrackParticles = std::make_pair(0,0);
  m_numCaloClusters   = std::make_pair(0,0);
  
}

//------------------------------------------------------------------------------
StatusCode gmsbOverlapRemovalTool::initialize() {

  ATH_MSG_DEBUG("in initialize()");

  m_debug = msg().level() <= MSG::DEBUG;

  /// get a handle on the selection tools
  StatusCode sc = m_userSelectionTool.retrieve();
  if ( sc.isFailure() ) {
      ATH_MSG_ERROR("Can't get handle on analysis selection tool");
      return sc;
  }

  sc = m_userOverlapCheckingTool.retrieve();
  if ( sc.isFailure() ) {
      ATH_MSG_ERROR("Can't get handle on analysis overlap checking tool");
      return sc;
  }

  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
StatusCode gmsbOverlapRemovalTool::finalize() {

  ATH_MSG_DEBUG("in finalize()");
 
  this->summarize();

  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
gmsbOverlapRemovalTool::~gmsbOverlapRemovalTool()
{}

//-------------------------------------------------------------------------------
StatusCode gmsbOverlapRemovalTool::execute() {
  ATH_MSG_DEBUG("in execute()");

  /** check if the execute is already called or not 
      in one job, execute should be called once for each event */
  if ( this->isExecuted() ) {
    ATH_MSG_WARNING("overlapRemovalTool->execute() already called for the event in this job");
    return StatusCode::SUCCESS; 
  }

  /** prepare the container for selection and overlap removal */
  StatusCode sc = this->prepareContainers();
  if ( sc.isFailure() ) return sc;
 
  /** now object preparation with overlap removal */
  for ( unsigned int i=0; i<m_inputContainerKeys.size(); ++i ) {

    std::string::size_type loc = m_inputContainerKeys[i].find( "Electron", 0);
    if( loc != std::string::npos ) sc = this->electronPreparation( m_inputContainerKeys[i] );

    loc = m_inputContainerKeys[i].find( "Photon", 0);
    if( loc != std::string::npos ) sc = this->photonPreparation( m_inputContainerKeys[i] );

    loc = m_inputContainerKeys[i].find( "Muon", 0);
    if( loc != std::string::npos ) sc = this->muonPreparation( m_inputContainerKeys[i] );

    std::string tau = "Tau";
    if ( m_isAtlfast ) tau = "TauJet";
    loc = m_inputContainerKeys[i].find( tau, 0);
    if( loc != std::string::npos ) sc = this->tauJetPreparation( m_inputContainerKeys[i] );

    std::string jet = "Jet";
    if ( m_isAtlfast ) jet = "Jet";
    loc = m_inputContainerKeys[i].find( jet, 0);
    if( loc != std::string::npos ) sc = this->jetPreparation( m_inputContainerKeys[i] );

    loc = m_inputContainerKeys[i].find( "Track", 0);
    if( loc != std::string::npos ) sc = this->trackParticlePreparation( m_inputContainerKeys[i] );

    loc = m_inputContainerKeys[i].find( "Cluster", 0);
    if( loc != std::string::npos ) sc = this->caloClusterPreparation( m_inputContainerKeys[i] );

    if ( sc.isFailure() ) return sc;

  }

  /** lock the containers so that they are no longer modified */
  sc = this->lockContainers();
  if ( sc.isFailure() ) return sc;

  if ( m_debug ) this->print();

  return StatusCode::SUCCESS;
}

//-------------------------------------------------------------------------------
const INavigable4MomentumCollection * gmsbOverlapRemovalTool::finalStateObjects() {
  ATH_MSG_DEBUG("in finalStateObjects()");
  const INavigable4MomentumCollection * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputObjectKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State objects not found");
  return container;
}

const PhotonContainer * gmsbOverlapRemovalTool::finalStatePhotons() {
  ATH_MSG_DEBUG("in finalStatePhotons()");
  const PhotonContainer * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputPhotonKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State Photons not found");
  return container;
}

const ElectronContainer * gmsbOverlapRemovalTool::finalStateElectrons() {
  ATH_MSG_DEBUG("in finalStateElectrons()");
  const ElectronContainer * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputElectronKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State Electrons not found");
  return container;
}

const MuonContainer * gmsbOverlapRemovalTool::finalStateMuons() {
  ATH_MSG_DEBUG("in finalStateMuons()");
  const MuonContainer * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputMuonKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State Muons not found");
  return container;
}

const INavigable4MomentumCollection * gmsbOverlapRemovalTool::finalStateLeptons() {
  ATH_MSG_DEBUG("in finalStateLeptons()");
  const INavigable4MomentumCollection * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputLeptonKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State Leptons not found");
  return container;
}

const TauJetContainer * gmsbOverlapRemovalTool::finalStateTauJets() {
  ATH_MSG_DEBUG("in finalStateTauJets()");
  const TauJetContainer * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputTauJetKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State TauJets not found");
  return container;
}

const JetCollection * gmsbOverlapRemovalTool::finalStateJets() {
  ATH_MSG_DEBUG("in finalStateJets()");
  const JetCollection * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputJetKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State Jets not found");
  return container;
}

const JetCollection * gmsbOverlapRemovalTool::finalStateBJets() {
  ATH_MSG_DEBUG("in finalStateBJets()");
  const JetCollection * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputBJetKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State BJets not found");
  return container;
}

const JetCollection * gmsbOverlapRemovalTool::finalStateLightJets() {
  ATH_MSG_DEBUG("in finalStateLightJets()");
  const JetCollection * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputLightJetKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State Light Jets not found");
  return container;
}

const TrackParticleContainer * gmsbOverlapRemovalTool::finalStateTrackParticles() {
  ATH_MSG_DEBUG("in finalStateTrackParticles()");
  const TrackParticleContainer * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputTrackParticleKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State TrackParticles not found");
  return container;
}

const CaloClusterContainer * gmsbOverlapRemovalTool::finalStateCaloClusters() {
  ATH_MSG_DEBUG("in finalStateCaloClusters()");
  const CaloClusterContainer * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputCaloClusterKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State CaloClusters not found");
  return container;
}

  /** container preparation */
StatusCode gmsbOverlapRemovalTool::electronPreparation( std::string key ) {
  ATH_MSG_DEBUG("in electronPreparation() with key = " << key);
  StatusCode sc = StatusCode::SUCCESS;

  INavigable4MomentumCollection * particles = this->allParticles();
  if ( !particles ) return sc;

  INavigable4MomentumCollection * leptons = this->allLeptons();
  if ( !leptons ) return sc;

  ElectronContainer * electrons = this->allElectrons();
  if ( !electrons ) return sc;

  const ElectronContainer * aod_electrons = 0;
  sc = evtStore()->retrieve( aod_electrons, key );
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("No Existing ESD/AOD/DPD electron container found: key = " << key);
    return sc;
  }
  ATH_MSG_DEBUG("Initial ElectronContainer size is " << aod_electrons->size());
  m_numElectrons.first += aod_electrons->size();

  /// iterators over the container 
  ElectronContainer::const_iterator elecItr  = aod_electrons->begin();
  ElectronContainer::const_iterator elecItrE = aod_electrons->end();

  for (; elecItr != elecItrE; ++elecItr) {

    /** if this is the first particle, just put it in */ 
    if ( particles->size() == 0 ) {
       particles->push_back( *elecItr );
       leptons->push_back( *elecItr );
       electrons->push_back( *elecItr );
    }   
    /** check for the overlap and save non overlapping ones */
    else {
      INavigable4MomentumCollection::const_iterator nav4MomItr  = particles->begin();
      INavigable4MomentumCollection::const_iterator nav4MomItrE = particles->end();
      bool overlap = false;
      for (; nav4MomItr != nav4MomItrE; ++nav4MomItr) {
          /** overlap checking */
          const Electron * electron = dynamic_cast<const Electron*>(*nav4MomItr);
          const egamma * eg = dynamic_cast<const egamma*>(*nav4MomItr);
          if ( !electron || ( electron && m_removeOverlapInSameContainer ) ) {
	    if (eg) {
	      overlap = m_userOverlapCheckingTool->overlap((*elecItr)->cluster(), eg->cluster());
	    } else {	      
	      overlap = m_userOverlapCheckingTool->overlap((*elecItr)->trackParticle(), *nav4MomItr);
	    }
	  }
          /** get out of the loop as soon as an overlap is found */
          if ( overlap ) break;
      }

      /** if no overlap then save */  
      if ( !overlap ) { 
         particles->push_back( *elecItr ); 
         leptons->push_back( *elecItr ); 
         electrons->push_back( *elecItr );
      }
    }
  }

  m_numElectrons.second += electrons->size();

  return sc;
}

StatusCode gmsbOverlapRemovalTool::photonPreparation( std::string key ) {
  ATH_MSG_DEBUG("in photonPreparation() with key = " << key);
  StatusCode sc = StatusCode::SUCCESS;

  INavigable4MomentumCollection * particles = this->allParticles();
  if ( !particles ) return sc;

  PhotonContainer * photons = this->allPhotons();
  if ( !photons ) return sc;

  const PhotonContainer * aod_photons = 0;
  sc = evtStore()->retrieve( aod_photons, key );
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("No Existing ESD/AOD/DPD photon container found: key = " << key);
    return sc;
  }
  ATH_MSG_DEBUG("Initial PhotonContainer size is " << aod_photons->size());
  m_numPhotons.first += aod_photons->size();

  /// iterators over the container 
  PhotonContainer::const_iterator photItr  = aod_photons->begin();
  PhotonContainer::const_iterator photItrE = aod_photons->end();

  for (; photItr != photItrE; ++photItr) {

    /** if this is the first particle, just put it in */ 
    if ( particles->size() == 0 ) {
       particles->push_back( *photItr );
       photons->push_back( *photItr );
    }   
    /** check for the overlap and save non overlapping ones */
    else {
      INavigable4MomentumCollection::const_iterator nav4MomItr  = particles->begin();
      INavigable4MomentumCollection::const_iterator nav4MomItrE = particles->end();
      bool overlap = false;
      for (; nav4MomItr != nav4MomItrE; ++nav4MomItr) {
          /** overlap checking */
          const Photon * photon = dynamic_cast<const Photon*>(*nav4MomItr);
          const egamma * eg = dynamic_cast<const egamma*>(*nav4MomItr);
          if ( !photon || ( photon && m_removeOverlapInSameContainer ) ) {
	    if (eg) {
	      overlap = m_userOverlapCheckingTool->overlap((*photItr)->cluster(), eg->cluster());
	    } else {
	      overlap = m_userOverlapCheckingTool->overlap(*photItr, *nav4MomItr);
	    }
	  }
          /** get out of the loop as soon as an overlap is found */
          if ( overlap ) break;
      }

      /** if no overlap then save */  
      if ( !overlap ) { 
         particles->push_back( *photItr );  
         photons->push_back( *photItr );
      }
    }
  }

  m_numPhotons.second += photons->size();

  return sc;
}

StatusCode gmsbOverlapRemovalTool::muonPreparation( std:: string key ) {
  ATH_MSG_DEBUG("in muonPreparation() ");
  StatusCode sc = StatusCode::SUCCESS;

  INavigable4MomentumCollection * particles = this->allParticles();
  if ( !particles ) return sc;

  MuonContainer * muons = this->allMuons();
  if ( !muons ) return sc;

  const MuonContainer * aod_muons = 0;
  sc = evtStore()->retrieve( aod_muons, key );
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("No Existing ESD/AOD/DPD muon container found: key = " << key);
    return sc; 
  }
  ATH_MSG_DEBUG("Initial MuonContainer size is " << aod_muons->size());
  m_numMuons.first += aod_muons->size();

  /// iterators over the container 
  MuonContainer::const_iterator muonItr  = aod_muons->begin();
  MuonContainer::const_iterator muonItrE = aod_muons->end();

  for (; muonItr != muonItrE; ++muonItr) {

    /** if this is the first particle, just put it in */ 
    if ( particles->size() == 0 ) {
      particles->push_back( *muonItr );
       muons->push_back( *muonItr );
    }   
    /** check for the overlap and save non overlapping ones */
    else {
      INavigable4MomentumCollection::const_iterator nav4MomItr  = particles->begin();
      INavigable4MomentumCollection::const_iterator nav4MomItrE = particles->end();
      bool overlap = false;
      for (; nav4MomItr != nav4MomItrE; ++nav4MomItr) {
          /** overlap checking */
          const Analysis::Muon * muon = dynamic_cast<const Analysis::Muon*>(*nav4MomItr);
          if ( !muon || ( muon && m_removeOverlapInSameContainer ) )  
             overlap = m_userOverlapCheckingTool->overlap(*muonItr, *nav4MomItr);

          /** get out of the loop as soon as an overlap is found */
          if ( overlap ) break;
      }

      /** if no overlap then save */  
      if ( !overlap ) { 
         particles->push_back( *muonItr ); 
         muons->push_back( *muonItr );
      }
    }
  }

  m_numMuons.second += muons->size();

  return sc;
}

StatusCode gmsbOverlapRemovalTool::tauJetPreparation( std::string key ) {
  ATH_MSG_DEBUG("in tauJetPreparation() with key = " << key);
  StatusCode sc = StatusCode::SUCCESS;

  INavigable4MomentumCollection * particles = this->allParticles();
  if ( !particles ) return sc;

  TauJetContainer * tauJets = this->allTauJets();
  if ( !tauJets ) return sc;

  const TauJetContainer * aod_tauJets = 0;
  sc = evtStore()->retrieve( aod_tauJets, key );
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("No Existing ESD/AOD/DPD tauJet container found: key = " << key);
    return sc;
  }
  ATH_MSG_DEBUG("initial TauJetContainer size is " << aod_tauJets->size());
  m_numTauJets.first += aod_tauJets->size();

  /// iterators over the container 
  TauJetContainer::const_iterator tauJetItr  = aod_tauJets->begin();
  TauJetContainer::const_iterator tauJetItrE = aod_tauJets->end();

  for (; tauJetItr != tauJetItrE; ++tauJetItr) {

    /** if this is the first particle, just put it in */ 
    if ( particles->size() == 0 ) {
      particles->push_back( *tauJetItr );
       tauJets->push_back( *tauJetItr );
    }   
    /** check for the overlap and save non overlapping ones */
    else {
      INavigable4MomentumCollection::const_iterator nav4MomItr  = particles->begin();
      INavigable4MomentumCollection::const_iterator nav4MomItrE = particles->end();
      bool overlap = false;
      for (; nav4MomItr != nav4MomItrE; ++nav4MomItr) {
          /** overlap checking */
          const TauJet * taujet = dynamic_cast<const TauJet*>(*nav4MomItr);
          if ( !taujet || ( taujet && m_removeOverlapInSameContainer ) )  
             overlap = m_userOverlapCheckingTool->overlap(*tauJetItr, *nav4MomItr);

          /** get out of the loop as soon as an overlap is found */
          if ( overlap ) break;
      }

      /** if no overlap then save */  
      if ( !overlap ) { 
         particles->push_back( *tauJetItr ); 
         tauJets->push_back( *tauJetItr );
      }
    }
  }

  m_numTauJets.second += tauJets->size();

  return sc;
}

StatusCode gmsbOverlapRemovalTool::jetPreparation( std::string key ) {
  ATH_MSG_DEBUG("in jetPreparation() with key = " << key);
  StatusCode sc = StatusCode::SUCCESS;

  INavigable4MomentumCollection * particles = this->allParticles();
  if ( !particles ) return sc;

  JetCollection * jets = this->allJets();
  if ( !jets ) return sc;

  JetCollection * bJets = this->allBJets();
  if ( !bJets ) return sc;

  JetCollection * lightJets = this->allLightJets();
  if ( !lightJets ) return sc;

  const JetCollection * aod_jets = 0;
  sc = evtStore()->retrieve( aod_jets, key );
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("No Existing ESD/AOD/DPD jet container found: key = " << key);
    return sc;
  }
  ATH_MSG_DEBUG("Initial JetCollection size is " << aod_jets->size());
  m_numJets.first      += aod_jets->size();
  m_numBJets.first     += aod_jets->size();
  m_numLightJets.first += aod_jets->size();

  /// iterators over the container 
  JetCollection::const_iterator jetItr  = aod_jets->begin();
  JetCollection::const_iterator jetItrE = aod_jets->end();

  for (; jetItr != jetItrE; ++jetItr) {
    /** check if this jet passes pre-selection */
    /** if this is the first particle, just put it in */ 
    if ( particles->size() == 0 ) {
       particles->push_back( *jetItr );
       jets->push_back( *jetItr );
       if ( m_userSelectionTool->isBJet( *jetItr ) ) bJets->push_back( *jetItr);
       else lightJets->push_back( *jetItr );
    }   
    /** check for the overlap and save non overlapping ones */
    else {
      INavigable4MomentumCollection::const_iterator nav4MomItr  = particles->begin();
      INavigable4MomentumCollection::const_iterator nav4MomItrE = particles->end();
      bool overlap = false;
      for (; nav4MomItr != nav4MomItrE; ++nav4MomItr) {
          /** overlap checking */
          const Jet * jet = dynamic_cast<const Jet*>(*nav4MomItr);
          const Electron * electron = dynamic_cast<const Electron*>(*nav4MomItr);
          if ( !jet || ( jet && m_removeOverlapInSameContainer ) ) {
	    if (electron) {
	      overlap = m_userOverlapCheckingTool->overlap(*jetItr, electron->trackParticle());
	    } else {
	      overlap = m_userOverlapCheckingTool->overlap(*jetItr, *nav4MomItr);
	    }
	  }
          /** get out of the loop as soon as an overlap is found */
          if ( overlap ) break;
      }

      /** if no overlap then save */  
      if ( !overlap ) { 
	particles->push_back( *jetItr );
	jets->push_back( *jetItr );
	if ( m_userSelectionTool->isBJet( *jetItr ) ) bJets->push_back( *jetItr);
	else lightJets->push_back( *jetItr );
      }
    }
  }

  m_numJets.second      += jets->size();
  m_numBJets.second     += bJets->size();
  m_numLightJets.second += lightJets->size();

  return sc;
}

StatusCode gmsbOverlapRemovalTool::trackParticlePreparation( std::string key ) {
  ATH_MSG_DEBUG("in trackParticlePreparation() ");
  StatusCode sc = StatusCode::SUCCESS;

  INavigable4MomentumCollection * particles = this->allParticles();
  if ( !particles ) return sc;

  TrackParticleContainer * trackParticles = this->allTrackParticles();
  if ( !trackParticles ) return sc;

  const TrackParticleContainer * aod_trackParticles = 0;
  sc = evtStore()->retrieve( aod_trackParticles, key );
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("No ESD/AOD/DPD TrackParticle container found: key = " << key);
    return sc;
  }
  ATH_MSG_DEBUG("Initial TrackParticleContainer size is " << aod_trackParticles->size());
  m_numTrackParticles.first += aod_trackParticles->size();

  /// iterators over the container 
  TrackParticleContainer::const_iterator trackParticleItr  = aod_trackParticles->begin();
  TrackParticleContainer::const_iterator trackParticleItrE = aod_trackParticles->end();

  for (; trackParticleItr != trackParticleItrE; ++trackParticleItr) {
    /** if this is the first particle, just put it in */ 
    if ( particles->size() == 0 ) {
      particles->push_back( *trackParticleItr );
       trackParticles->push_back( *trackParticleItr );
    }   
    /** check for the overlap and save non overlapping ones */
    else {
      INavigable4MomentumCollection::const_iterator nav4MomItr  = particles->begin();
      INavigable4MomentumCollection::const_iterator nav4MomItrE = particles->end();
      bool overlap = false; 
      for (; nav4MomItr != nav4MomItrE; ++nav4MomItr) {
          /** overlap checking */
          const TrackParticle * trackparticle = dynamic_cast<const TrackParticle*>(*nav4MomItr);
          if ( !trackparticle || ( trackparticle && m_removeOverlapInSameContainer ) )  
             overlap = m_userOverlapCheckingTool->overlap(*trackParticleItr, *nav4MomItr);
          /** get out of the loop as soon as an overlap is found */
          if ( overlap ) break;
      }

      /** if no overlap then save */  
      if ( !overlap ) { 
         particles->push_back( *trackParticleItr ); 
         trackParticles->push_back( *trackParticleItr );
      }
    }
  }
  
  m_numTrackParticles.second += trackParticles->size();

  return sc;
}

StatusCode gmsbOverlapRemovalTool::caloClusterPreparation( std::string key ) {
  ATH_MSG_DEBUG("in caloClusterPreparation() ");
  StatusCode sc = StatusCode::SUCCESS;

  INavigable4MomentumCollection * particles = this->allParticles();
  if ( !particles ) return sc;

  CaloClusterContainer * caloClusters = this->allCaloClusters();
  if ( !caloClusters ) return sc;

  const CaloClusterContainer * aod_caloClusters = 0;
  sc = evtStore()->retrieve( aod_caloClusters, key );
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("No Existing ESD/AOD/DPD CaloCluster container found: key = " << key);
    return sc;
  }
  ATH_MSG_DEBUG("Initial CaloClusterContainer size is " << aod_caloClusters->size());
  m_numCaloClusters.first += aod_caloClusters->size();

  /// iterators over the container 
  CaloClusterContainer::const_iterator caloClusterItr  = aod_caloClusters->begin();
  CaloClusterContainer::const_iterator caloClusterItrE = aod_caloClusters->end();

  for (; caloClusterItr != caloClusterItrE; ++caloClusterItr) {
    /** check if this caloCluster passes pre-selection */
    if ( !m_userSelectionTool->isSelected( *caloClusterItr ) ) continue;

    /** if this is the first particle, just put it in */ 
    if ( particles->size() == 0 ) {
      particles->push_back( *caloClusterItr );
       caloClusters->push_back( *caloClusterItr );
    }   
    /** check for the overlap and save non overlapping ones */
    else {
      INavigable4MomentumCollection::const_iterator nav4MomItr  = particles->begin();
      INavigable4MomentumCollection::const_iterator nav4MomItrE = particles->end();
      bool overlap = false;
      for (; nav4MomItr != nav4MomItrE; ++nav4MomItr) {
          /** overlap checking */
          const CaloCluster * cluster = dynamic_cast<const CaloCluster*>(*nav4MomItr);
          if ( !cluster || ( cluster && m_removeOverlapInSameContainer ) )  
             overlap = m_userOverlapCheckingTool->overlap(*caloClusterItr, *nav4MomItr);
          /** get out of the loop as soon as an overlap is found */
          if ( overlap ) break;
      }

      /** if no overlap then save */  
      if ( !overlap ) { 
         particles->push_back( *caloClusterItr ); 
         caloClusters->push_back( *caloClusterItr );
      } 
    }
  }

  m_numCaloClusters.second += caloClusters->size();

  return sc;
}

INavigable4MomentumCollection * gmsbOverlapRemovalTool::allParticles() {
  ATH_MSG_DEBUG("in allObjects()");
  INavigable4MomentumCollection * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputObjectKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State objects not found");
  return container;
}

INavigable4MomentumCollection * gmsbOverlapRemovalTool::allLeptons() {
  ATH_MSG_DEBUG("in allLeptons()");
  INavigable4MomentumCollection * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputLeptonKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State leptons not found");
  return container;

}

PhotonContainer * gmsbOverlapRemovalTool::allPhotons() {
  ATH_MSG_DEBUG("in allPhotons()");
  PhotonContainer * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputPhotonKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State Photons not found");
  return container;
}

MuonContainer * gmsbOverlapRemovalTool::allMuons() {
  ATH_MSG_DEBUG("in allMuons()");
  MuonContainer * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputMuonKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State Muons not found");
  return container;
}

ElectronContainer * gmsbOverlapRemovalTool::allElectrons() {
  ATH_MSG_DEBUG("in allElectrons()");
  ElectronContainer * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputElectronKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State Electrons not found");
  return container;
}

TauJetContainer * gmsbOverlapRemovalTool::allTauJets() {
  ATH_MSG_DEBUG("in allTauJets()");
  TauJetContainer * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputTauJetKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State TauJets not found");
  return container;
}

JetCollection * gmsbOverlapRemovalTool::allJets() {
  ATH_MSG_DEBUG("in allJets()");
  JetCollection * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputJetKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State Jets not found");
  return container;
}

JetCollection * gmsbOverlapRemovalTool::allBJets() {
  ATH_MSG_DEBUG("in allBJets()");
  JetCollection * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputBJetKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State BJets not found");
  return container;
}

JetCollection * gmsbOverlapRemovalTool::allLightJets() {
  ATH_MSG_DEBUG("in allLightJets()");
  JetCollection * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputLightJetKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State Light Jets not found");
  return container;
}

TrackParticleContainer * gmsbOverlapRemovalTool::allTrackParticles() {
  ATH_MSG_DEBUG("in allTrackParticles()");
  TrackParticleContainer * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputTrackParticleKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State TrackParticles not found");
  return container;
}

CaloClusterContainer * gmsbOverlapRemovalTool::allCaloClusters() {
  ATH_MSG_DEBUG("in allCaloClusters()");
  CaloClusterContainer * container = 0;
  StatusCode sc = evtStore()->retrieve(container, m_outputCaloClusterKey);
  if ( sc.isFailure() || container ==0 )
     ATH_MSG_ERROR("Final State CaloClusters not found");
  return container;
}

//-------------------------------------------------------------------------------
StatusCode gmsbOverlapRemovalTool::prepareContainers() {
  ATH_MSG_DEBUG("in prepareContainers()");

  /** create an empty container of all particles and record it */
  CaloClusterContainer * caloClusters = new CaloClusterContainer( SG::VIEW_ELEMENTS );
  StatusCode sc = evtStore()->record ( caloClusters, m_outputCaloClusterKey);
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("Not able to create a collection of CaloClusters in StoreGate: key= " << m_outputCaloClusterKey);
     return sc;
  }

  /** create an empty container of TrackParticles and record it */
  TrackParticleContainer * trackParticles = new TrackParticleContainer( SG::VIEW_ELEMENTS );
  sc = evtStore()->record ( trackParticles, m_outputTrackParticleKey);
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("Not able to create a collection of TrackParticles in StoreGate: key= " << m_outputTrackParticleKey);
     return sc;
  }

  /** create an empty container of all particles and record it */
  INavigable4MomentumCollection * particles = new INavigable4MomentumCollection( SG::VIEW_ELEMENTS );
  sc = evtStore()->record ( particles, m_outputObjectKey);
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("Not able to create a collection of particles in StoreGate: key=  " << m_outputObjectKey);
     return sc; 
  }
  
  /** create an empty container of all leptons and record it */
  INavigable4MomentumCollection * leptons = new INavigable4MomentumCollection( SG::VIEW_ELEMENTS );
  sc = evtStore()->record ( leptons, m_outputLeptonKey);
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("Not able to create a collection of leptons in StoreGate: key= " << m_outputLeptonKey);
     return sc;
  }

  /** create an empty container of all electrons and record it */
  ElectronContainer * electrons = new ElectronContainer( SG::VIEW_ELEMENTS );
  sc = evtStore()->record ( electrons, m_outputElectronKey);
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("Not able to create a collection of electrons in StoreGate: key= " << m_outputElectronKey);
     return sc;
  }

  /** create an empty container of all photons and record it */
  PhotonContainer * photons = new PhotonContainer( SG::VIEW_ELEMENTS );
  sc = evtStore()->record ( photons, m_outputPhotonKey);
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("Not able to create a collection of photons in StoreGate: key= " << m_outputPhotonKey);
     return sc;
  }

  /** create an empty container of all muons and record it */
  MuonContainer * muons = new MuonContainer( SG::VIEW_ELEMENTS );
  sc = evtStore()->record ( muons, m_outputMuonKey);
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("Not able to create a collection of muons in StoreGate: key= " << m_outputMuonKey);
     return sc;
  }

  /** create an empty container of all tauJets and record it */
  TauJetContainer * tauJets = new TauJetContainer( SG::VIEW_ELEMENTS );
  sc = evtStore()->record ( tauJets, m_outputTauJetKey);
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("Not able to create a collection of tau jets in StoreGate: key= " << m_outputTauJetKey);
     return sc;
  }

  /** create an empty container of all jets and record it */
  JetCollection * jets = new JetCollection( SG::VIEW_ELEMENTS );
  sc = evtStore()->record ( jets, m_outputJetKey);
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("Not able to create a collection of jets in StoreGate: key= " << m_outputJetKey);
     return sc;
  }

  /** create an empty container of b-jets and record it */
  JetCollection * bjets = new JetCollection( SG::VIEW_ELEMENTS );
  sc = evtStore()->record ( bjets, m_outputBJetKey);
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("Not able to create a collection of b-jets in StoreGate: key= " << m_outputBJetKey);
     return sc;
  }

  /** create an empty container of light (non b-jet) jets and record it */
  JetCollection * lightJets = new JetCollection( SG::VIEW_ELEMENTS );
  sc = evtStore()->record ( lightJets, m_outputLightJetKey);
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("Not able to create a collection of lightJets in StoreGate: key= " << m_outputLightJetKey);
     return sc;
  }

  return StatusCode::SUCCESS;
}

//-------------------------------------------------------------------------------
StatusCode gmsbOverlapRemovalTool::lockContainers() {
  ATH_MSG_DEBUG("in lockContainers()");

  /** lock the contianer so it is not modified downstream by anyone else */
  StatusCode sc = evtStore()->setConst( this->allParticles() );
  if ( sc.isFailure() ) ATH_MSG_WARNING("Not able to lock the container of particles ");

  sc = evtStore()->setConst( this->allLeptons() );
  if ( sc.isFailure() ) ATH_MSG_WARNING("Not able to lock the container of leptons ");

  sc = evtStore()->setConst( this->allElectrons() );
  if ( sc.isFailure() ) ATH_MSG_WARNING("Not able to lock the container of electrons ");

  sc = evtStore()->setConst( this->allPhotons() );
  if ( sc.isFailure() ) ATH_MSG_WARNING("Not able to lock the container of photons ");

  sc = evtStore()->setConst( this->allMuons() );
  if ( sc.isFailure() ) ATH_MSG_WARNING("Not able to lock the container of muons ");

  sc = evtStore()->setConst( this->allTauJets() );
  if ( sc.isFailure() ) ATH_MSG_WARNING("Not able to lock the container of tauJets ");

  sc = evtStore()->setConst( this->allJets() );
  if ( sc.isFailure() ) ATH_MSG_WARNING("Not able to lock the container of jets ");

  sc = evtStore()->setConst( this->allBJets() );
  if ( sc.isFailure() ) ATH_MSG_WARNING("Not able to lock the container of b-jets ");

  sc = evtStore()->setConst( this->allLightJets() );
  if ( sc.isFailure() ) ATH_MSG_WARNING("Not able to lock the container of light Jets ");

  sc = evtStore()->setConst( this->allTrackParticles() );
  if ( sc.isFailure() ) ATH_MSG_WARNING("Not able to lock the container of trackParticles ");

  sc = evtStore()->setConst( this->allCaloClusters() );
  if ( sc.isFailure() ) ATH_MSG_WARNING("Not able to lock the container of calo clusters ");

  return sc;
}

//-----------------------------------------------------------------------------------------------
void gmsbOverlapRemovalTool::print() {
  ATH_MSG_DEBUG("in print() ");

  /** Get the container of pre-selected Electrons */
  const ElectronContainer * electrons = this->finalStateElectrons();
  ATH_MSG_DEBUG("Number of Pre-selected Electrons is " << electrons->size());

  /** Get the container of pre-selected Photons */
  const PhotonContainer * photons = this->finalStatePhotons();
  ATH_MSG_DEBUG("Number of Pre-selected Photons is " << photons->size());

  /** Get the container of pre-selected Muons */
  const MuonContainer * muons = this->finalStateMuons();
  ATH_MSG_DEBUG("Number of Pre-selected Muons is " << muons->size());

  /** Get the container of pre-selected TauJets */
  const TauJetContainer * tauJets = this->finalStateTauJets();
  ATH_MSG_DEBUG("Number of Pre-selected TauJets is " << tauJets->size());

  /** Get the container of pre-selected Jets */
  const JetCollection * jets = this->finalStateJets();
  ATH_MSG_DEBUG("Number of Pre-selected Jets is " << jets->size());

  /** Get the container of pre-selected B-tagged Jets */
  const JetCollection * bjets = this->finalStateBJets();
  ATH_MSG_DEBUG("Number of Pre-selected b-Jets is " << bjets->size());

  /** Get the container of pre-selected non b-jets */
  const JetCollection * lightJets = this->finalStateLightJets();
  ATH_MSG_DEBUG("Number of Pre-selected LightJets is " << lightJets->size());

  /** Get the container of pre-selected TrackParticles */
  const TrackParticleContainer * trackParticles = this->finalStateTrackParticles();
  ATH_MSG_DEBUG("Number of Pre-selected TrackParticles is " << trackParticles->size());

  /** Get the container of pre-selected CaloClusters */
  const CaloClusterContainer * caloClusters = this->finalStateCaloClusters();
  ATH_MSG_DEBUG("Number of Pre-selected CaloClusters is " << caloClusters->size());

  /** Get the container of pre-selected leptons (electrons, muons) */
  const INavigable4MomentumCollection * leptons = this->finalStateLeptons();
  ATH_MSG_DEBUG("Number of Pre-selected Leptons is " << leptons->size());

  /** Get the container of ALL pre-selected objects */
  const INavigable4MomentumCollection * allObjects = this->finalStateObjects();
  ATH_MSG_DEBUG("Number of Pre-selected final State Objects is " << allObjects->size());

}

//---------------------------------------------------------------------------------------------------------
void gmsbOverlapRemovalTool::summarize() {

  ATH_MSG_INFO("Summary Pre-selected/Overlap Removed Events ###################");
  ATH_MSG_INFO("---------------------------------------------------------------");
  ATH_MSG_INFO("Pre-selected Electrons             = " << std::setw(10) << m_numElectrons.first 
                      << "   Overlap-removed Electrons       = " << std::setw(10) << m_numElectrons.second);
  ATH_MSG_INFO("Pre-selected Photons               = " << std::setw(10) << m_numPhotons.first 
                      << "   Overlap-removed Photons         = " << std::setw(10) << m_numPhotons.second);
  ATH_MSG_INFO("Pre-selected pMuons                = " << std::setw(10) << m_numMuons.first 
                      << "    Overlap-removed Muons          = " << std::setw(10) << m_numMuons.second);
  ATH_MSG_INFO("Pre-selected TauJets               = " << std::setw(10) << m_numTauJets.first  
                      << "    Overlap-removed TauJets        = " << std::setw(10) << m_numTauJets.second);
  ATH_MSG_INFO("Pre-selected Jets                  = " << std::setw(10) << m_numJets.first 
                      << "    Overlap-removed Jets           = " << std::setw(10) << m_numJets.second);
  ATH_MSG_INFO("Pre-selected BJets                 = " << std::setw(10) << m_numBJets.first  
                      << "    Overlap-removed BJets          = " << std::setw(10) << m_numBJets.second);
  ATH_MSG_INFO("Pre-selected LightJets             = " << std::setw(10) << m_numLightJets.first
                      << "    Overlpa-removed LightJets      = " << std::setw(10) << m_numLightJets.second);
  ATH_MSG_INFO("Pre-selected TrackParticles        = " << std::setw(10) << m_numTrackParticles.first
                      << "    Overlap-removed TrackParticles = " << std::setw(10) << m_numTrackParticles.second);
  ATH_MSG_INFO("Pre-selected CaloClusters          = " << std::setw(10) << m_numCaloClusters.first
                      << "   Overlap-removed CaloClusters    = " << std::setw(10) << m_numCaloClusters.second);
}

//-----------------------------------------------------------------------------------------------------------
bool gmsbOverlapRemovalTool::isExecuted() {
  ATH_MSG_DEBUG("in isExecuted() ");
  return evtStore()->contains<INavigable4MomentumCollection>( m_outputObjectKey );
}

