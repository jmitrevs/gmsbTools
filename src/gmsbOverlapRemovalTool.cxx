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

#include "gmsbTools/SortHelpers.h"

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

  declareProperty("OutputPhotonKey",        m_outputPhotonKey        = "FinalStatePhotons");
  declareProperty("OutputElectronKey",      m_outputElectronKey      = "FinalStateElectrons");
  declareProperty("OutputMuonKey",          m_outputMuonKey          = "FinalStateMuons");
  declareProperty("OutputJetKey",           m_outputJetKey           = "FinalStateJets");
  // declareProperty("OutputBJetKey",          m_outputBJetKey          = "FinalStateBJets");
  // declareProperty("OutputLightJetKey",      m_outputLightJetKey      = "FinalStateLightJets");
  declareProperty("RemoveOverlapInSameContainer", m_removeOverlapInSameContainer=false);
  declareProperty("ElPhKillMuon", m_elPhKillMuon=false);

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

  StatusCode sc = StatusCode::SUCCESS;

  m_allParticles.clear();

  /** now object preparation with overlap removal */
  for ( unsigned int i=0; i<m_inputContainerKeys.size(); ++i ) {

    ATH_MSG_DEBUG("m_inputContainerKeys[" << i << "] = " << m_inputContainerKeys[i]);

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

  return sc;
}

//-------------------------------------------------------------------------------

PhotonD3PDObject * gmsbOverlapRemovalTool::finalStatePhotons() {
  ATH_MSG_DEBUG("in finalStatePhotons()");
  PhotonD3PDObject * container = new PhotonD3PDObject(m_outputPhotonKey);
  StatusCode sc = container->retrieve();
  if (!sc.isSuccess()) {
    delete container;
    container = NULL;
  }
  return container;
}

ElectronD3PDObject * gmsbOverlapRemovalTool::finalStateElectrons() {
  ATH_MSG_DEBUG("in finalStateElectrons()");
  ElectronD3PDObject * container = new ElectronD3PDObject(m_outputElectronKey);
  StatusCode sc = container->retrieve();
  if (!sc.isSuccess()) {
    delete container;
    container = NULL;
  }
  return container;
}

MuonD3PDObject * gmsbOverlapRemovalTool::finalStateMuons() {
  ATH_MSG_DEBUG("in finalStateMuons()");
  MuonD3PDObject * container = new MuonD3PDObject(m_outputMuonKey);
  StatusCode sc = container->retrieve();
  if (!sc.isSuccess()) {
    delete container;
    container = NULL;
  }
  return container;
}

JetD3PDObject * gmsbOverlapRemovalTool::finalStateJets() {
  ATH_MSG_DEBUG("in finalStateJets()");
  JetD3PDObject * container = new JetD3PDObject(m_outputJetKey);
  StatusCode sc = container->retrieve();
  if (!sc.isSuccess()) {
    delete container;
    container = NULL;
  }
  return container;
}

// JetD3PDObject * gmsbOverlapRemovalTool::finalStateBJets() {
//   ATH_MSG_DEBUG("in finalStateBJets()");
//   JetD3PDObject * container = new JetD3PDObject(m_outputBJetKey);
//   StatusCode sc = container->retrieve();
//   if (!sc.isSuccess()) {
//     delete container;
//     container = NULL;
//   }
// }

// JetD3PDObject * gmsbOverlapRemovalTool::finalStateLightJets() {
//   ATH_MSG_DEBUG("in finalStateLightJets()");
//   JetD3PDObject * container = new JetD3PDObject(m_outputLightJetKey);
//   StatusCode sc = container->retrieve();
//   if (!sc.isSuccess()) {
//     delete container;
//     container = NULL;
//   }
// }

  /** container preparation */
StatusCode gmsbOverlapRemovalTool::electronPreparation( std::string key ) {
  ATH_MSG_DEBUG("in electronPreparation() with key = " << key);

  ElectronD3PDObject electrons = ElectronD3PDObject::create(m_outputElectronKey);
  ATH_CHECK(electrons.record());

  const ElectronD3PDObject aod_electrons(key);
  ATH_CHECK(aod_electrons.retrieve());

  ATH_MSG_DEBUG("AOD ElectronD3PDObject size is " << aod_electrons.n());
  m_numElectrons.first += aod_electrons.n();

  // need to first sort 
  SortHelpers::sl_t sortedList;
  SortHelpers::sort(sortedList, aod_electrons);

  for (SortHelpers::sl_t::const_iterator it = sortedList.begin(); it != sortedList.end(); ++it) {
    const std::size_t idx = it->first;

    bool overlap = false;
    for (std::vector<genericParticle>::iterator particle = m_allParticles.begin();
	 particle != m_allParticles.end(); 
	 ++particle) {
      /** overlap checking */
      switch (particle->type()) {
      case genericParticle::electron:
	if (m_removeOverlapInSameContainer) {
	  overlap = m_userOverlapCheckingTool->overlap(aod_electrons.cl_eta(idx), 
						       aod_electrons.cl_phi(idx),
						       particle->cl_eta(), 
						       particle->cl_phi(), false);
	}
	break;
      case genericParticle::photon:
	overlap = m_userOverlapCheckingTool->overlap(aod_electrons.cl_eta(idx), 
						     aod_electrons.cl_phi(idx),
						     particle->cl_eta(), 
						     particle->cl_phi(), false);
	break;
      case genericParticle::muon:
	overlap = m_userOverlapCheckingTool->overlap(aod_electrons.eta(idx), 
						     aod_electrons.phi(idx),
						     particle->eta(), 
						     particle->phi(), false);
	break;
      default:
	overlap = m_userOverlapCheckingTool->overlap(aod_electrons.cl_eta(idx), 
						     aod_electrons.cl_phi(idx),
						     particle->eta(), 
						     particle->phi(), true);
	break;
      }
      /** get out of the loop as soon as an overlap is found */
      if ( overlap ) break;
    }
    
    /** if no overlap then save */  
    if ( !overlap ) { 
      electrons.add_object(aod_electrons, idx);
      m_allParticles.push_back(genericParticle(aod_electrons.pt(idx), 
					       aod_electrons.eta(idx), 
					       aod_electrons.phi(idx),
					       genericParticle::electron,
					       aod_electrons.cl_eta(idx), 
					       aod_electrons.cl_phi(idx))); 
    }
  }

  m_numElectrons.second += electrons.n();

  return StatusCode::SUCCESS;
}

StatusCode gmsbOverlapRemovalTool::photonPreparation( std::string key ) {
  ATH_MSG_DEBUG("in photonPreparation() with key = " << key);

  PhotonD3PDObject photons = PhotonD3PDObject::create(m_outputPhotonKey);
  ATH_CHECK(photons.record());

  const PhotonD3PDObject aod_photons(key);
  ATH_CHECK(aod_photons.retrieve());

  ATH_MSG_DEBUG("AOD PhotonD3PDObject size is " << aod_photons.n());
  m_numPhotons.first += aod_photons.n();

  // need to first sort 
  SortHelpers::sl_t sortedList;
  SortHelpers::sort(sortedList, aod_photons);

  for (SortHelpers::sl_t::const_iterator it = sortedList.begin(); it != sortedList.end(); ++it) {
    const std::size_t idx = it->first;

    bool overlap = false;
    for (std::vector<genericParticle>::iterator particle = m_allParticles.begin();
	 particle != m_allParticles.end(); 
	 ++particle) {
      /** overlap checking */
      switch (particle->type()) {
      case genericParticle::photon:
	if (m_removeOverlapInSameContainer) {
	  overlap = m_userOverlapCheckingTool->overlap(aod_photons.cl_eta(idx), 
						       aod_photons.cl_phi(idx),
						       particle->cl_eta(), 
						       particle->cl_phi(), false);
	}
	break;
      case genericParticle::electron:
	overlap = m_userOverlapCheckingTool->overlap(aod_photons.cl_eta(idx), 
						     aod_photons.cl_phi(idx),
						     particle->cl_eta(), 
						     particle->cl_phi(), false);
	break;
      case genericParticle::muon:
	overlap = m_userOverlapCheckingTool->overlap(aod_photons.eta(idx), 
						     aod_photons.phi(idx),
						     particle->eta(), 
						     particle->phi(), false);
	break;
      default:
	overlap = m_userOverlapCheckingTool->overlap(aod_photons.cl_eta(idx), 
						     aod_photons.cl_phi(idx),
						     particle->eta(), 
						     particle->phi(), true);
	break;
      }
      /** get out of the loop as soon as an overlap is found */
      if ( overlap ) break;
    }
    
    /** if no overlap then save */  
    if ( !overlap ) { 
      photons.add_object(aod_photons, idx);
      m_allParticles.push_back(genericParticle(aod_photons.pt(idx), 
					       aod_photons.eta(idx), 
					       aod_photons.phi(idx),
					       genericParticle::photon,
					       aod_photons.cl_eta(idx), 
					       aod_photons.cl_phi(idx))); 
    }
  }

  m_numPhotons.second += photons.n();

  return StatusCode::SUCCESS;
}

StatusCode gmsbOverlapRemovalTool::muonPreparation( std:: string key ) {
  ATH_MSG_DEBUG("in muonPreparation() with key = " << key);

  MuonD3PDObject muons = MuonD3PDObject::create(m_outputMuonKey);
  ATH_CHECK(muons.record());

  const MuonD3PDObject aod_muons(key);
  ATH_CHECK(aod_muons.retrieve());

  ATH_MSG_DEBUG("AOD MuonD3PDObject size is " << aod_muons.n());
  m_numMuons.first += aod_muons.n();

  // need to first sort 
  SortHelpers::sl_t sortedList;
  SortHelpers::sort(sortedList, aod_muons);

  for (SortHelpers::sl_t::const_iterator it = sortedList.begin(); it != sortedList.end(); ++it) {
    const std::size_t idx = it->first;

    bool overlap = false;
    for (std::vector<genericParticle>::iterator particle = m_allParticles.begin();
	 particle != m_allParticles.end(); 
	 ++particle) {
      /** overlap checking */
      switch (particle->type()) {
      case genericParticle::muon:
	if (m_removeOverlapInSameContainer) {
	  overlap = m_userOverlapCheckingTool->overlap(aod_muons.eta(idx), 
						       aod_muons.phi(idx),
						       particle->eta(), 
						       particle->phi(), false);
	}
	break;
      case genericParticle::electron:
      case genericParticle::photon:
	if (m_elPhKillMuon) {
	  overlap = m_userOverlapCheckingTool->overlap(aod_muons.eta(idx), 
						       aod_muons.phi(idx),
						       particle->eta(), 
						       particle->phi(), false);
	}
	break;
      default:
	overlap = m_userOverlapCheckingTool->overlap(aod_muons.eta(idx), 
						     aod_muons.phi(idx),
						     particle->eta(), 
						     particle->phi(), true);
	break;
      }
      /** get out of the loop as soon as an overlap is found */
      if ( overlap ) break;
    }
    
    /** if no overlap then save */  
    if ( !overlap ) { 
      muons.add_object(aod_muons, idx);
      m_allParticles.push_back(genericParticle(aod_muons.pt(idx), 
					       aod_muons.eta(idx), 
					       aod_muons.phi(idx),
					       genericParticle::muon));
    }
  }

  m_numMuons.second += muons.n();

  return StatusCode::SUCCESS;
}

StatusCode gmsbOverlapRemovalTool::jetPreparation( std::string key ) {
  ATH_MSG_DEBUG("in jetPreparation() with key = " << key);

  JetD3PDObject jets = JetD3PDObject::create(m_outputJetKey);
  ATH_CHECK(jets.record());

  const JetD3PDObject aod_jets(key);
  ATH_CHECK(aod_jets.retrieve());

  ATH_MSG_DEBUG("AOD JetD3PDObject size is " << aod_jets.n());
  m_numJets.first += aod_jets.n();

  // need to first sort 
  SortHelpers::sl_t sortedList;
  SortHelpers::sort(sortedList, aod_jets);

  for (SortHelpers::sl_t::const_iterator it = sortedList.begin(); it != sortedList.end(); ++it) {
    const std::size_t idx = it->first;

    bool overlap = false;
    for (std::vector<genericParticle>::iterator particle = m_allParticles.begin();
	 particle != m_allParticles.end(); 
	 ++particle) {
      /** overlap checking */
      switch (particle->type()) {
      case genericParticle::jet:
	if (m_removeOverlapInSameContainer) {
	  overlap = m_userOverlapCheckingTool->overlap(aod_jets.eta(idx), 
						       aod_jets.phi(idx),
						       particle->eta(), 
						       particle->phi(), true);
	}
	break;
      case genericParticle::electron:
      case genericParticle::photon:
	overlap = m_userOverlapCheckingTool->overlap(aod_jets.eta(idx), 
						     aod_jets.phi(idx),
						     particle->cl_eta(), 
						     particle->cl_phi(), true);
	break;
      default:
	overlap = m_userOverlapCheckingTool->overlap(aod_jets.eta(idx), 
						     aod_jets.phi(idx),
						     particle->eta(), 
						     particle->phi(), true);
	break;
      }
      /** get out of the loop as soon as an overlap is found */
      if ( overlap ) break;
    }
    
    /** if no overlap then save */  
    if ( !overlap ) { 
      jets.add_object(aod_jets, idx);
      m_allParticles.push_back(genericParticle(aod_jets.pt(idx), 
					       aod_jets.eta(idx), 
					       aod_jets.phi(idx),
					       genericParticle::jet));
    }
  }

  m_numJets.second += jets.n();

  return StatusCode::SUCCESS;
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
  ATH_MSG_INFO("Pre-selected Muons                 = " << std::setw(10) << m_numMuons.first 
                      << "    Overlap-removed Muons          = " << std::setw(10) << m_numMuons.second);
  ATH_MSG_INFO("Pre-selected Jets                  = " << std::setw(10) << m_numJets.first 
                      << "    Overlap-removed Jets           = " << std::setw(10) << m_numJets.second);
}

