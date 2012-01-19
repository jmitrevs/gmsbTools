#include "gmsbTools/TruthStudies.h"

#include "EventInfo/EventInfo.h"
#include "EventInfo/EventID.h"

#include "egammaEvent/ElectronContainer.h"
#include "egammaEvent/Electron.h"
#include "egammaEvent/PhotonContainer.h"
#include "egammaEvent/Photon.h"
#include "egammaEvent/egammaPIDdefs.h"
#include "egammaEvent/EMShower.h"

#include "McParticleEvent/TruthParticleContainer.h"
#include "FourMomUtils/P4Helpers.h"
#include "GeneratorObjects/McEventCollection.h"


/////////////////////////////////////////////////////////////////////////////
TruthStudies::TruthStudies(const std::string& type,
			   const std::string& name, 
			   const IInterface* parent )
  : AthAlgTool( type, name, parent )
{
  declareInterface<TruthStudies>( this );

  declareProperty("McParticleContainer", m_truthParticleContainerName = "SpclMC");
  declareProperty("McEventCollection", m_mcEventCollectionName = "TruthEvent");
  declareProperty("PrintDecayTree", m_printDecayTree = false);
  declareProperty("UseAnnotated", m_useAnnotated = false);
  declareProperty("DumpEntierTree", m_dumpEntireTree = false);

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
StatusCode TruthStudies::initialize(){

  ATH_MSG_DEBUG("initialize()");
 
  return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
StatusCode TruthStudies::execute() 
{
  ATH_MSG_DEBUG("execute");

  StatusCode sc = StatusCode::SUCCESS;
  
  const HepMC::GenEvent *ge = 0;

  /** get the MC truth particle AOD or ESD container from StoreGate */
  const McEventCollection *mcEventCol = 0;
  if (evtStore()->contains<McEventCollection>(m_mcEventCollectionName)) {
    sc=evtStore()->retrieve( mcEventCol, m_mcEventCollectionName);
    if( sc.isFailure()  ||  !mcEventCol) {
      ATH_MSG_ERROR("could not retrieve MC event container");
      return StatusCode::RECOVERABLE;
    }
    ATH_MSG_DEBUG("McEventCollection found with name " <<  m_mcEventCollectionName);
    ge = mcEventCol->at(0);
  }
  
  if (!ge) {
    const TruthParticleContainer*  mcpartTES = 0;
    if (evtStore()->contains<TruthParticleContainer>(m_truthParticleContainerName)) {
      sc=evtStore()->retrieve( mcpartTES, m_truthParticleContainerName);
      if( sc.isFailure()  ||  !mcpartTES ) {
	ATH_MSG_ERROR("could not retrieve MC truth container");
	return StatusCode::RECOVERABLE;
      }
      ATH_MSG_DEBUG("McEventCollection found with name " <<  m_mcEventCollectionName);
      ge=mcpartTES->genEvent();
    }
  }

  // const EventInfo*  evtInfo = 0;
  // sc = evtStore()->retrieve(evtInfo);
  // if(sc.isFailure() || !evtInfo) {
  //   ATH_MSG_ERROR("could not retrieve event info");
  //   return StatusCode::RECOVERABLE;
  // }

  // const unsigned runNum = evtInfo->event_ID()->run_number();
  // //const unsigned lbNum = evtInfo->event_ID()->lumi_block();
  // //const unsigned evNum = evtInfo->event_ID()->event_number();


  //mLog <<MSG::DEBUG << "ge = " << (unsigned int) ge << endreq;

  const HepMC::GenVertex *pvtx = NULL;

  if (ge) {
    pvtx = ge->signal_process_vertex();
    // mLog <<MSG::DEBUG << "pvtx from signal_process_vertex = " << (unsigned int) pvtx << endreq;

    if (!pvtx) {
      pvtx = getMCHardInteraction(ge);
      // mLog <<MSG::DEBUG << "pvtx from getMCHardInteraction = " << (unsigned int) pvtx << endreq;
    }
    
    if (pvtx) {

      if (m_printDecayTree) {
	for (HepMC::GenVertex::particles_in_const_iterator init = pvtx->particles_in_const_begin();
	     init != pvtx->particles_in_const_end();
	     init++) {
	  msg(MSG::INFO) << std::setw(7) << std::left <<  m_pdg.GetParticle((*init)->pdg_id())->GetName();
	}
	msg(MSG::INFO) << " ->\t";
      }
      if (m_useAnnotated) {
	FollowDecayTreeAnnotated(pvtx);
      } else {
	FollowDecayTree(pvtx);
      }
    }
    
    if (m_dumpEntireTree) DumpEntireTree(ge);

  }
  
  return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
StatusCode TruthStudies::finalize() {
    
    ATH_MSG_DEBUG ("finalize()");
    return StatusCode::SUCCESS;
}

// This is how the hard interaction is gotten
// for now it's a heuristic and probably doesn't work in most cases
HepMC::GenVertex* TruthStudies::getMCHardInteraction(const HepMC::GenEvent *const ge) const
{
  for (int i = -1; i > -20; --i) {
    HepMC::GenVertex* vtx = ge->barcode_to_vertex(i);
    if (vtx && vtx->particles_in_size() == 2) {
      return vtx;
    }
  }
  return NULL;
}

// used recursively
// prints the decay products of one vertex, calling itself to to
// further decay of SUSY particles, t, W, Z, higgses, (gamma? not now).
// Doesn't continue into Geant particles 
void TruthStudies::FollowDecayTree(const HepMC::GenVertex *vtx, int extraSpaces)
{
  std::vector<const HepMC::GenVertex *> decayVertices;

  for (HepMC::GenVertex::particles_in_const_iterator outit = vtx->particles_out_const_begin();
       outit != vtx->particles_out_const_end();
       outit++) {
    if (StatusGood((*outit)->status())) {
      msg(MSG::INFO) << std::setw(11) << std::left << m_pdg.GetParticle((*outit)->pdg_id())->GetName();
      const HepMC::GenVertex *nextVertex = FindNextVertex(*outit);
      if (nextVertex) {
	decayVertices.push_back(FindNextVertex(*outit));
      } else if ((*outit)->pdg_id() == 22) {
	//CalcTruthStudies((*outit)->momentum());
      }
    }
  }
  msg(MSG::INFO) << endreq;
  for (int i = decayVertices.size(); i > 0; --i) {
    int index = i-1;
    if (decayVertices.at(index) != NULL) {
      msg(MSG::INFO) << "                 \t";
      for (int j = 0; j < index+extraSpaces; j++) {
	msg(MSG::INFO) << "                ";
      }
      FollowDecayTree(decayVertices.at(index), index+extraSpaces);
    }
  }
}

// used recursively
// prints the decay products of one vertex, calling itself to to
// further decay of SUSY particles, t, W, Z, higgses, (gamma? not now).
// Doesn't continue into Geant particles 
void TruthStudies::FollowDecayTreeAnnotated(const HepMC::GenVertex *vtx, int extraSpaces) const
{
  std::vector<const HepMC::GenVertex *> decayVertices;
  
  //std::cout << "Working on vertex with barcode: " << vtx->barcode() << std::endl;

  for (HepMC::GenVertex::particles_in_const_iterator outit = vtx->particles_out_const_begin();
       outit != vtx->particles_out_const_end();
       outit++) {
    //    if (StatusGood((*outit)->status())) {
    if (1) {

      HepMC::FourVector p = (*outit)->momentum();
      //msg(MSG::INFO) << std::setw(4) << std::right << round(p.perp()/GeV) << " ";
      msg(MSG::INFO) << std::setw(4) << std::right << (*outit)->status() << " ";
      msg(MSG::INFO) << std::setw(11) << std::left << m_pdg.GetParticle((*outit)->pdg_id())->GetName();
      decayVertices.push_back(FindNextVertex(*outit));
    }
  }
  // msg(MSG::INFO) << "\t vertex = " << vtx->barcode();
  msg(MSG::INFO) << endreq;
  for (int i = decayVertices.size(); i > 0; --i) {
    int index = i-1;
    if (decayVertices.at(index) != NULL) {
      msg(MSG::INFO) << "                 \t";
      for (int j = 0; j < index+extraSpaces; j++) {
	msg(MSG::INFO) << "                ";
      }
      FollowDecayTreeAnnotated(decayVertices.at(index), index+extraSpaces);
    }
  }
}

const HepMC::GenVertex *TruthStudies::FindNextVertex(const HepMC::GenParticle *pcl) const
{
  if (pcl->barcode() > 200000) return NULL;

  int pid = abs(pcl->pdg_id());

  if ((pid > 22 && pid < 38) || 
      (pid == 6) ||
      (pid > 1000000 && pid < 1000040) || 
      (pid > 2000000 && pid < 2000016)) { 

    // only show decay products of SUSY and massive guage particles and top

    const HepMC::GenVertex *nextVtx = pcl->end_vertex();
    while (nextVtx != NULL && nextVtx->particles_out_size() == 1) {
      const HepMC::GenParticle *np = *(nextVtx->particles_out_const_begin());
      if (np->barcode() > 200000) return NULL;
      nextVtx = np->end_vertex();
    }
    if (nextVtx && nextVtx->barcode() < -200000) return NULL;
    return nextVtx;
  } else {
    return NULL;
  }
}

void TruthStudies::DumpEntireTree(const HepMC::GenEvent *ge) const
{
  for(HepMC::GenEvent::particle_const_iterator pitr=ge->particles_begin();
      pitr!=ge->particles_end(); ++pitr ){
    if( (*pitr)->status()==1 ) {
      ATH_MSG_INFO("Found particle of type " << m_pdg.GetParticle((*pitr)->pdg_id())->GetName() 
		   << " with pT = " << (*pitr)->momentum().perp());
    }
  }
}


void TruthStudies::FillEventType(decayType d1, decayType d2)
{

  decayType first;
  decayType second;

  if (d1 <= d2) {
    first = d1;
    second = d2;
  } else {
    first = d2;
    second = d1;
  }

  switch (first) {

  case gamma:
    switch (second) {
    case gamma:
      m_type = gammagamma;
      break;
    case Zee:
      m_type = Zeegamma;
      break;
    case Zmumu:
      m_type = Zmumugamma;
      break;
    case Ztautau:
      m_type = Ztautaugamma;
      break;
    case Zjj:
      m_type = Zjjgamma;
      break;
    case Wenu:
      m_type = Wenugamma;
      break;
    case Wmunu:
      m_type = Wmunugamma;
      break;
    case Wtaunu:
      m_type = Wtaunugamma;
      break;
    case Wjj:
      m_type = Wjjgamma;
      break;
    default:
      ATH_MSG_FATAL("Received unexpected decay");
      break;
    }
    break;

  case Zee:
    switch (second) {
    case Zee:
      m_type = ZeeZee;
      break;
    case Zmumu:
      m_type = ZmumuZee;
      break;
    case Ztautau:
      m_type = ZtautauZee;
      break;
    case Zjj:
      m_type = ZjjZee;
      break;
    case Wenu:
      m_type = WenuZee;
      break;
    case Wmunu:
      m_type = WmunuZee;
      break;
    case Wtaunu:
      m_type = WtaunuZee;
      break;
    case Wjj:
      m_type = WjjZee;
      break;
    default:
      ATH_MSG_FATAL("Received unexpected decay");
      break;
    }
    break;

  case Zmumu:
    switch (second) {
    case Zmumu:
      m_type = ZmumuZmumu;
      break;
    case Ztautau:
      m_type = ZtautauZmumu;
      break;
    case Zjj:
      m_type = ZjjZmumu;
      break;
    case Wenu:
      m_type = WenuZmumu;
      break;
    case Wmunu:
      m_type = WmunuZmumu;
      break;
    case Wtaunu:
      m_type = WtaunuZmumu;
      break;
    case Wjj:
      m_type = WjjZmumu;
      break;
    default:
      ATH_MSG_FATAL("Received unexpected decay");
      break;
    }
    break;

  case Ztautau:
    switch (second) {
    case Ztautau:
      m_type = ZtautauZtautau;
      break;
    case Zjj:
      m_type = ZjjZtautau;
      break;
    case Wenu:
      m_type = WenuZtautau;
      break;
    case Wmunu:
      m_type = WmunuZtautau;
      break;
    case Wtaunu:
      m_type = WtaunuZtautau;
      break;
    case Wjj:
      m_type = WjjZtautau;
      break;
    default:
      ATH_MSG_FATAL("Received unexpected decay");
      break;
    }
    break;

  case Zjj:
    switch (second) {
    case Zjj:
      m_type = ZjjZjj;
      break;
    case Wenu:
      m_type = WenuZjj;
      break;
    case Wmunu:
      m_type = WmunuZjj;
      break;
    case Wtaunu:
      m_type = WtaunuZjj;
      break;
    case Wjj:
      m_type = WjjZjj;
      break;
    default:
      ATH_MSG_FATAL("Received unexpected decay");
      break;
    }
    break;

  case Wenu:
    switch (second) {
    case Wenu:
      m_type = WenuWenu;
      break;
    case Wmunu:
      m_type = WmunuWenu;
      break;
    case Wtaunu:
      m_type = WtaunuWenu;
      break;
    case Wjj:
      m_type = WjjWenu;
      break;
    default:
      ATH_MSG_FATAL("Received unexpected decay");
      break;
    }
    break;

  case Wmunu:
    switch (second) {
    case Wmunu:
      m_type = WmunuWmunu;
      break;
    case Wtaunu:
      m_type = WtaunuWmunu;
      break;
    case Wjj:
      m_type = WjjWmunu;
      break;
    default:
      ATH_MSG_FATAL("Received unexpected decay");
      break;
    }
    break;

  case Wtaunu:
    switch (second) {
    case Wtaunu:
      m_type = WtaunuWtaunu;
      break;
    case Wjj:
      m_type = WjjWtaunu;
      break;
    default:
      ATH_MSG_FATAL("Received unexpected decay");
      break;
    }
    break;

  case Wjj:
    switch (second) {
    case Wjj:
      m_type = WjjWtaunu;
      break;
    default:
      ATH_MSG_FATAL("Received unexpected decay");
      break;
    }
    break;

  default:
    ATH_MSG_FATAL("Received unexpected decay");
    break;
  }
}
