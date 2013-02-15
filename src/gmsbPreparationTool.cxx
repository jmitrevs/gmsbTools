/*****************************************************************************
Name    : gmsbPreparationTool.cxx
Package : offline/PhysicsAnalysis/AnalysisCommon/gmsbTools
Author  : Ketevi A. Assamagan
Created : November 2007
Purpose : User Analysis Preparation - see gmsbPreparationTool.h for details
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
#include "gmsbTools/gmsbPreparationTool.h"
#include "AnalysisUtils/AnalysisMisc.h"

#include <sstream>
#include <iomanip>
#include <iostream>

//------------------------------------------------------------------------------
gmsbPreparationTool::gmsbPreparationTool( const std::string& type,
					  const std::string& name, 
					  const IInterface* parent )
  : AthAlgTool( type, name, parent ),
    m_userSelectionTool ( "gmsbSelectionTool" ) {

  declareInterface<gmsbPreparationTool>( this );

  declareProperty("UserSelectionTool",      m_userSelectionTool);
  declareProperty("InputContainerKeys",     m_inputContainerKeys);
  declareProperty("OutputContainerKeys",    m_outputContainerKeys);
  declareProperty("IsAtlfastData",          m_isAtlfast=false);

  /** initialize counters */
  m_numElectrons      = std::make_pair(0,0);
  m_numPhotons        = std::make_pair(0,0);
  m_numMuons          = std::make_pair(0,0);
  m_numJets           = std::make_pair(0,0);
  m_first             = true;
}

//------------------------------------------------------------------------------
StatusCode gmsbPreparationTool::initialize() {

  ATH_MSG_DEBUG("in initialize()");

  m_debug = msg().level() <= MSG::DEBUG;

  /// get a handle on the selection tools
  StatusCode sc = m_userSelectionTool.retrieve();
  if ( sc.isFailure() ) {
    ATH_MSG_ERROR("Can't get handle on analysis selection tool");
    return sc;
  }

  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
StatusCode gmsbPreparationTool::finalize() {

  ATH_MSG_DEBUG("in finalize()");
  
  this->summarize();
  
  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
gmsbPreparationTool::~gmsbPreparationTool()
{}

//-------------------------------------------------------------------------------
StatusCode gmsbPreparationTool::execute() {
  ATH_MSG_DEBUG("in execute()");

  /** check that the input and the output containers are defined */
  StatusCode sc = StatusCode::SUCCESS;

  if ( m_first ) {
    if ( m_outputContainerKeys.size() != m_inputContainerKeys.size() ) {
      ATH_MSG_FATAL("Input/Output container mis-match: please fix job options");
      return StatusCode::FAILURE;
    }
    if ( m_outputContainerKeys.size() == 0 ) {
      ATH_MSG_ERROR("You should input at least one container : please fix job options");
      return StatusCode::FAILURE;
    }
  }
  

  /** now object preparation with selection */
  for ( unsigned int i=0; i<m_inputContainerKeys.size(); ++i ) {

    std::string::size_type loc = m_inputContainerKeys[i].find( "el_", 0);
    if ( loc != std::string::npos ) { 
      if ( m_first ) m_outputElectronKey = m_outputContainerKeys[i]; 
      sc = this->electronPreparation( m_inputContainerKeys[i]);
    }

    loc = m_inputContainerKeys[i].find( "ph_", 0);
    if ( loc != std::string::npos ) {
      if ( m_first ) m_outputPhotonKey = m_outputContainerKeys[i];       
      sc = this->photonPreparation( m_inputContainerKeys[i]);
    }

    loc = m_inputContainerKeys[i].find( "mu_", 0);
    if ( loc != std::string::npos ) {
      if ( m_first ) m_outputMuonKey = m_outputContainerKeys[i];       
      sc = this->muonPreparation( m_inputContainerKeys[i] );
    }

    loc = m_inputContainerKeys[i].find( "jet_", 0);
    if ( loc != std::string::npos ) { 
      if ( m_first ) m_outputJetKey = m_outputContainerKeys[i];
      sc = this->jetPreparation( m_inputContainerKeys[i] );
    }
    if ( sc.isFailure() ) return sc;
  }

  if ( m_debug ) this->print();
  m_first = false;

  return StatusCode::SUCCESS;
}

//-------------------------------------------------------------------------------
const PhotonD3PDObject * gmsbPreparationTool::selectedPhotons() {
  ATH_MSG_DEBUG("in selectedPhotons()");
  const PhotonD3PDObject * container = new PhotonD3PDObject(m_outputPhotonKey);
  StatusCode sc = container->retrieve();
  if (!sc.isSuccess()) {
    delete container;
    container = NULL;
  }
  return container;
}

const ElectronD3PDObject * gmsbPreparationTool::selectedElectrons() {
  ATH_MSG_DEBUG("in selectedElectrons()");
  const ElectronD3PDObject * container = new ElectronD3PDObject(m_outputElectronKey);
  StatusCode sc = container->retrieve();
  if (!sc.isSuccess()) {
    delete container;
    container = NULL;
  }
  return container;
}

const MuonD3PDObject * gmsbPreparationTool::selectedMuons() {
  ATH_MSG_DEBUG("in selectedMuons()");
  const MuonD3PDObject * container = new MuonD3PDObject(m_outputMuonKey);
  StatusCode sc = container->retrieve();
  if (!sc.isSuccess()) {
    delete container;
    container = NULL;
  }
  return container;
}

const JetD3PDObject * gmsbPreparationTool::selectedJets() {
  ATH_MSG_DEBUG("in selectedJets()");
  const JetD3PDObject * container = new JetD3PDObject(m_outputJetKey);
  StatusCode sc = container->retrieve();
  if (!sc.isSuccess()) {
    delete container;
    container = NULL;
  }
  return container;
}


/** container preparation */
StatusCode gmsbPreparationTool::electronPreparation( std::string key) {
  ATH_MSG_DEBUG("in electronPreparation() ");

  ElectronD3PDObject electrons(m_outputElectronKey);
  ATH_CHECK(electrons.record());

  ElectronD3PDObject aod_electrons(key);
  ATH_CHECK(aod_electrons.retrieve());

  ATH_MSG_DEBUG("AOD ElectronD3PDObject size is " << aod_electrons.n());
  m_numElectrons.first += aod_electrons.n();

  for (std::size_t idx = 0; idx < static_cast<std::size_t>(aod_electrons.n()); idx++) {
    if ( m_userSelectionTool->isSelected(aod_electrons, idx) ) {
      electrons.add_object(aod_electrons, idx);
    }
  }
  m_numElectrons.second += electrons.n();
  
  return StatusCode::SUCCESS;
}

StatusCode gmsbPreparationTool::photonPreparation( std::string key) {
  ATH_MSG_DEBUG("in photonPreparation() ");

  PhotonD3PDObject photons(m_outputPhotonKey);
  ATH_CHECK(photons.record());

  PhotonD3PDObject aod_photons(key);
  ATH_CHECK(aod_photons.retrieve());

  ATH_MSG_DEBUG("AOD PhotonD3PDObject size is " << aod_photons.n());
  m_numPhotons.first += aod_photons.n();

  for (std::size_t idx = 0; idx < static_cast<std::size_t>(aod_photons.n()); idx++) {
    if ( m_userSelectionTool->isSelected(aod_photons, idx) ) {
      photons.add_object(aod_photons, idx);
    }
  }
  m_numPhotons.second += photons.n();
  
  return StatusCode::SUCCESS;
}

StatusCode gmsbPreparationTool::muonPreparation( std::string key) {
  ATH_MSG_DEBUG("in muonPreparation() ");

  MuonD3PDObject muons(m_outputMuonKey);
  ATH_CHECK(muons.record());

  MuonD3PDObject aod_muons(key);
  ATH_CHECK(aod_muons.retrieve());

  ATH_MSG_DEBUG("AOD MuonD3PDObject size is " << aod_muons.n());
  m_numMuons.first += aod_muons.n();

  for (std::size_t idx = 0; idx < static_cast<std::size_t>(aod_muons.n()); idx++) {
    if ( m_userSelectionTool->isSelected(aod_muons, idx) ) {
      muons.add_object(aod_muons, idx);
    }
  }
  m_numMuons.second += muons.n();
  
  return StatusCode::SUCCESS;
}

StatusCode gmsbPreparationTool::jetPreparation( std::string key) {
  ATH_MSG_DEBUG("in jetPreparation() ");

  JetD3PDObject jets(m_outputJetKey);
  ATH_CHECK(jets.record());

  JetD3PDObject aod_jets(key);
  ATH_CHECK(aod_jets.retrieve());

  ATH_MSG_DEBUG("AOD JetD3PDObject size is " << aod_jets.n());
  m_numJets.first += aod_jets.n();

  for (std::size_t idx = 0; idx < static_cast<std::size_t>(aod_jets.n()); idx++) {
    if ( m_userSelectionTool->isSelected(aod_jets, idx) ) {
      jets.add_object(aod_jets, idx);
    }
  }
  m_numJets.second += jets.n();
  
  return StatusCode::SUCCESS;
}


//-----------------------------------------------------------------------------------------------
void gmsbPreparationTool::print() {
  ATH_MSG_DEBUG("in print() ");

  /** Get the container of pre-selected Electrons */
  const ElectronD3PDObject * electrons = this->selectedElectrons();
  if(electrons) ATH_MSG_DEBUG("Number of Pre-selected Electrons is " << electrons->n());
  delete electrons;

  /** Get the container of pre-selected Photons */
  const PhotonD3PDObject * photons = this->selectedPhotons();
  if(photons) ATH_MSG_DEBUG("Number of Pre-selected Photons is " << photons->n());
  delete photons;

  /** Get the container of pre-selected Muons */
  const MuonD3PDObject * muons = this->selectedMuons();
  if(muons)ATH_MSG_DEBUG("Number of Pre-selected Muons is " << muons->n());
  delete muons;

  /** Get the container of pre-selected Jets */
  const JetD3PDObject * jets = this->selectedJets();
  if(jets) ATH_MSG_DEBUG("Number of Pre-selected Jets is " << jets->n());
  delete jets;
}

//---------------------------------------------------------------------------------------------------------
void gmsbPreparationTool::summarize() {
  ATH_MSG_INFO("in summarize() ");

  ATH_MSG_INFO("Summary of Reconstructed Events/pre-selected events ############");
  ATH_MSG_INFO("---------------------------------------------------------------");
  ATH_MSG_INFO("Reconstructed Electrons        = " << std::setw(10) << m_numElectrons.first 
	       << "   Pre-selected Electrons      = " << std::setw(10) << m_numElectrons.second);
  ATH_MSG_INFO("Reconstructed Photons          = " << std::setw(10) << m_numPhotons.first 
	       << "   Pre-selected Photons        = " << std::setw(10) << m_numPhotons.second);
  ATH_MSG_INFO("Reconstructed Muons            = " << std::setw(10) << m_numMuons.first 
	       << "   Pre-selected Muons          = " << std::setw(10) << m_numMuons.second);
  ATH_MSG_INFO("Reconstructed Jets             = " << std::setw(10) << m_numJets.first 
	       << "   Pre-selected Jets           = " << std::setw(10) << m_numJets.second);
}


