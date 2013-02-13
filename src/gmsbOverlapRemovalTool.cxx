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
#include "gmsbD3PDObjects/ElectronD3PDObject.h"
#include "gmsbD3PDObjects/MuonD3PDObject.h"
#include "gmsbD3PDObjects/JetD3PDObject.h"
#include "gmsbD3PDObjects/PhotonD3PDObject.h"

// User Tools
#include "gmsbTools/gmsbOverlapRemovalTool.h"

#include <sstream>
#include <iomanip>
#include <iostream>



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
  declareProperty("OutputJetKey",           m_outputJetKey           = "FinalStateJets");
  // declareProperty("OutputBJetKey",          m_outputBJetKey          = "FinalStateBJets");
  // declareProperty("OutputLightJetKey",      m_outputLightJetKey      = "FinalStateLightJets");
  declareProperty("RemoveOverlapInSameContainer", m_removeOverlapInSameContainer=true);

  /** initialize counters */
  m_numElectrons      = std::make_pair(0,0);
  m_numPhotons        = std::make_pair(0,0);
  m_numMuons          = std::make_pair(0,0);
  m_numJets           = std::make_pair(0,0);
  // m_numBJets          = std::make_pair(0,0);
  // m_numLightJets      = std::make_pair(0,0);
  
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

  /** now object preparation with overlap removal */
  for ( unsigned int i=0; i<m_inputContainerKeys.size(); ++i ) {

    std::string::size_type loc = m_inputContainerKeys[i].find( "el_", 0);
    if( loc != std::string::npos ) sc = this->electronPreparation( m_inputContainerKeys[i] );

    loc = m_inputContainerKeys[i].find( "ph_", 0);
    if( loc != std::string::npos ) sc = this->photonPreparation( m_inputContainerKeys[i] );

    loc = m_inputContainerKeys[i].find( "mu_", 0);
    if( loc != std::string::npos ) sc = this->muonPreparation( m_inputContainerKeys[i] );

    loc = m_inputContainerKeys[i].find( "jet_", 0);
    if( loc != std::string::npos ) sc = this->jetPreparation( m_inputContainerKeys[i] );

    if ( sc.isFailure() ) return sc;

  }

  if ( m_debug ) this->print();

  return StatusCode::SUCCESS;
}

//-------------------------------------------------------------------------------

const PhotonD3PDObject * gmsbOverlapRemovalTool::finalStatePhotons() {
  ATH_MSG_DEBUG("in finalStatePhotons()");
  const PhotonD3PDObject * container = new PhotonD3PDObject(m_outputPhotonKey);
  StatusCode sc = container.retrieve();
  if (!sc.isSuccess()) {
    delete container;
    container = NULL;
  }
  return container;
}

const ElectronD3PDObject * gmsbOverlapRemovalTool::finalStateElectrons() {
  ATH_MSG_DEBUG("in finalStateElectrons()");
  const ElectronD3PDObject * container = new ElectronD3PDObject(m_outputElectronKey);
  StatusCode sc = container.retrieve();
  if (!sc.isSuccess()) {
    delete container;
    container = NULL;
  }
  return container;
}

const MuonD3PDObject * gmsbOverlapRemovalTool::finalStateMuons() {
  ATH_MSG_DEBUG("in finalStateMuons()");
  const MuonD3PDObject * container = new MuonD3PDObject(m_outputMuonKey);
  StatusCode sc = container.retrieve();
  if (!sc.isSuccess()) {
    delete container;
    container = NULL;
  }
}

const JetD3PDObject * gmsbOverlapRemovalTool::finalStateJets() {
  ATH_MSG_DEBUG("in finalStateJets()");
  const JetD3PDObject * container = new JetD3PDObject(m_outputJetKey);
  StatusCode sc = container.retrieve();
  if (!sc.isSuccess()) {
    delete container;
    container = NULL;
  }
}

// const JetD3PDObject * gmsbOverlapRemovalTool::finalStateBJets() {
//   ATH_MSG_DEBUG("in finalStateBJets()");
//   const JetD3PDObject * container = new JetD3PDObject(m_outputBJetKey);
//   StatusCode sc = container.retrieve();
//   if (!sc.isSuccess()) {
//     delete container;
//     container = NULL;
//   }
// }

// const JetD3PDObject * gmsbOverlapRemovalTool::finalStateLightJets() {
//   ATH_MSG_DEBUG("in finalStateLightJets()");
//   const JetD3PDObject * container = new JetD3PDObject(m_outputLightJetKey);
//   StatusCode sc = container.retrieve();
//   if (!sc.isSuccess()) {
//     delete container;
//     container = NULL;
//   }
// }

  /** container preparation */
StatusCode gmsbOverlapRemovalTool::electronPreparation( std::string key ) {
  ATH_MSG_DEBUG("in electronPreparation() with key = " << key);

  ElectronD3PDObject electrons(m_outputElectronKey);
  ATH_CHECK(electrons.record());

  ElectronD3PDObject aod_electrons(key);
  ATH_CHECK(aod_electrons.retrieve());

  ATH_MSG_DEBUG("AOD ElectronD3PDObject size is " << aod_electrons.n());
  m_numElectrons.first += aod_electrons.n();


  for (std::size_t idx = 0; idx < aod_electrons.n(); idx++) {

    ///////////////////////////////// STOPPED HERE /////////////////

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

  PhotonD3PDObject * photons = this->allPhotons();
  if ( !photons ) return sc;

  const PhotonD3PDObject * aod_photons = 0;
  sc = evtStore()->retrieve( aod_photons, key );
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("No Existing ESD/AOD/DPD photon container found: key = " << key);
    return sc;
  }
  ATH_MSG_DEBUG("Initial PhotonD3PDObject size is " << aod_photons->size());
  m_numPhotons.first += aod_photons->size();

  /// iterators over the container 
  PhotonD3PDObject::const_iterator photItr  = aod_photons->begin();
  PhotonD3PDObject::const_iterator photItrE = aod_photons->end();

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

  MuonD3PDObject * muons = this->allMuons();
  if ( !muons ) return sc;

  const MuonD3PDObject * aod_muons = 0;
  sc = evtStore()->retrieve( aod_muons, key );
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("No Existing ESD/AOD/DPD muon container found: key = " << key);
    return sc; 
  }
  ATH_MSG_DEBUG("Initial MuonD3PDObject size is " << aod_muons->size());
  m_numMuons.first += aod_muons->size();

  /// iterators over the container 
  MuonD3PDObject::const_iterator muonItr  = aod_muons->begin();
  MuonD3PDObject::const_iterator muonItrE = aod_muons->end();

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


StatusCode gmsbOverlapRemovalTool::jetPreparation( std::string key ) {
  ATH_MSG_DEBUG("in jetPreparation() with key = " << key);
  StatusCode sc = StatusCode::SUCCESS;

  INavigable4MomentumCollection * particles = this->allParticles();
  if ( !particles ) return sc;

  JetD3PDObject * jets = this->allJets();
  if ( !jets ) return sc;

  JetD3PDObject * bJets = this->allBJets();
  if ( !bJets ) return sc;

  JetD3PDObject * lightJets = this->allLightJets();
  if ( !lightJets ) return sc;

  const JetD3PDObject * aod_jets = 0;
  sc = evtStore()->retrieve( aod_jets, key );
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("No Existing ESD/AOD/DPD jet container found: key = " << key);
    return sc;
  }
  ATH_MSG_DEBUG("Initial JetD3PDObject size is " << aod_jets->size());
  m_numJets.first      += aod_jets->size();
  m_numBJets.first     += aod_jets->size();
  m_numLightJets.first += aod_jets->size();

  /// iterators over the container 
  JetD3PDObject::const_iterator jetItr  = aod_jets->begin();
  JetD3PDObject::const_iterator jetItrE = aod_jets->end();

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


//-----------------------------------------------------------------------------------------------
void gmsbOverlapRemovalTool::print() {
  ATH_MSG_DEBUG("in print() ");

  /** Get the container of pre-selected Electrons */
  const ElectronD3PDObject * electrons = this->finalStateElectrons();
  ATH_MSG_DEBUG("Number of Pre-selected Electrons is " << electrons->n());
  delete electrons;

  /** Get the container of pre-selected Photons */
  const PhotonD3PDObject * photons = this->finalStatePhotons();
  ATH_MSG_DEBUG("Number of Pre-selected Photons is " << photons->n());
  delete photons;

  /** Get the container of pre-selected Muons */
  const MuonD3PDObject * muons = this->finalStateMuons();
  ATH_MSG_DEBUG("Number of Pre-selected Muons is " << muons->n());
  delete muons;

  /** Get the container of pre-selected Jets */
  const JetD3PDObject * jets = this->finalStateJets();
  ATH_MSG_DEBUG("Number of Pre-selected Jets is " << jets->n());
  delete jets;

  // /** Get the container of pre-selected B-tagged Jets */
  // const JetD3PDObject * bjets = this->finalStateBJets();
  // ATH_MSG_DEBUG("Number of Pre-selected b-Jets is " << bjets->n());
  // delete bjets;

  // /** Get the container of pre-selected non b-jets */
  // const JetD3PDObject * lightJets = this->finalStateLightJets();
  // ATH_MSG_DEBUG("Number of Pre-selected LightJets is " << lightJets->n());
  // delete lightJets;

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

