#include "gmsbTools/TruthStudies.h"

#include "gmsbD3PDObjects/TruthParticleD3PDObject.h"
//#include "gmsbD3PDObjects/GenEventD3PDObject.h"

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/LorentzVector.h"

#include <iostream>

/////////////////////////////////////////////////////////////////////////////
TruthStudies::TruthStudies(const std::string& type,
			   const std::string& name, 
			   const IInterface* parent )
  : AthAlgTool( type, name, parent )
{
  declareInterface<TruthStudies>( this );

  declareProperty("McParticleContainer", m_truthParticleContainerName = "mc_");
  //declareProperty("McEventCollection", m_mcEventCollectionName = "mcevt_");
  declareProperty("PrintDecayTree", m_printDecayTree = false);
  declareProperty("UseAnnotated", m_useAnnotated = false);
  declareProperty("DumpEntireTree", m_dumpEntireTree = false);

  // when counting photons
  declareProperty("Ptcut",m_Ptmin = 80*GeV);
  declareProperty("Etacut",m_EtaRange = 5.0);
  declareProperty("doDeltaRLepton", m_doDeltaRLepton = false);
  declareProperty("doDeltaRLight", m_doDeltaRLight = false);
  declareProperty("doMInv", m_doMInv = false);
  declareProperty("DeltaRLepton", m_deltaRLepton = 0.1);
  declareProperty("DeltaRLight", m_deltaRLight = 0.1);
  declareProperty("MInv", m_mInv = 5*GeV);
  declareProperty("decayTaus", m_decayTaus = false);
  declareProperty("WptID", m_WptID = 24); // can change it to 23 for Z

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

  m_decays.clear();
  m_type = unknown;
  m_isStrong = false;

  m_leptons.clear();
  m_lightParticles.clear();

  m_Wpt = -999.;

  const TruthParticleD3PDObject truthObj(m_truthParticleContainerName);
  ATH_CHECK(truthObj.retrieve());


  if (truthObj.n()) {
    // int obj = 0;

    // if (m_printDecayTree) {
    //   for (HepMC::GenVertex::particles_in_const_iterator init = pvtx->particles_in_const_begin();
    // 	   init != pvtx->particles_in_const_end();
    // 	   init++) {
    // 	const TParticlePDG* pdgpt = m_pdg.GetParticle((*init)->pdg_id());
    // 	if (pdgpt) {
    // 	  msg(MSG::INFO) << std::setw(7) << std::left <<  pdgpt->GetName();
    // 	} else {
    // 	  msg(MSG::INFO) << std::setw(7) << std::left <<  (*init)->pdg_id();
    // 	}
    //   }
    //   msg(MSG::INFO) << " ->\t";
    // }
    // if (m_useAnnotated) {
    //   FollowDecayTreeAnnotated(pvtx);
    // } else {
    //   FollowDecayTree(pvtx);
    // }
    
    if (m_dumpEntireTree) DumpEntireTree(truthObj);

    // m_nPhotons = findPhotons(ge);

    // findElectrons(ge);

    // FillEventType();
    // ATH_MSG_DEBUG("Truth type = " << m_type << "; strong = " << m_isStrong << ", num truth photons = " << m_nPhotons);

  }
  
  return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
StatusCode TruthStudies::finalize() {
    
    ATH_MSG_DEBUG ("finalize()");
    return StatusCode::SUCCESS;
}

/*
// This is how the hard interaction is gotten
// for now it's a heuristic and probably doesn't work in most cases
HepMC::GenVertex* TruthStudies::getMCHardInteraction(const HepMC::GenEvent *const ge) const
{
  for (int i = -1; i > -20; --i) {
    HepMC::GenVertex* vtx = ge->barcode_to_vertex(i);
    //if (vtx && vtx->particles_in_size() == 2) {
      return vtx;
      //}
  }
  return NULL;
}

// used recursively
// prints the decay products of one vertex, calling itself to to
// further decay of SUSY particles, t, W, Z, higgses, (gamma? not now).
// Doesn't continue into Geant particles 
void TruthStudies::FollowDecayTree(const HepMC::GenVertex *vtx, int extraSpaces, int haveSeen)
{
  std::vector<const HepMC::GenVertex *> decayVertices;
  std::vector<int> pids;

  int newHaveSeen = 0;
  double pt = 0; 		// for precedence, take higher pT one

  for (HepMC::GenVertex::particles_in_const_iterator outit = vtx->particles_out_const_begin();
       outit != vtx->particles_out_const_end();
       outit++) {
    if (StatusGood((*outit)->status())) {
    //if (1) {
      // msg(MSG::INFO) << (*outit)->status() << " ";
      const int pid = (*outit)->pdg_id();
      int abspid = abs(pid);
      if (m_printDecayTree) {
	const TParticlePDG* pdgpt = m_pdg.GetParticle(pid);
	if (pdgpt) {
	  msg(MSG::INFO) << std::setw(11) << std::left <<  pdgpt->GetName();
	} else {
	  msg(MSG::INFO) << std::setw(11) << std::left <<  pid;
	}
      }

      //msg(MSG::INFO) << std::setw(4) << std::right << round(p.perp()/GeV) << " ";

      if (abspid == m_WptID) {
	// maybe add to the Wpt
	const HepMC::FourVector p = (*outit)->momentum();
	const double newPt = p.perp();
	if (newPt > m_Wpt) m_Wpt = newPt;
      }

      if (haveSeen == 1000022) { // neutralino
	const HepMC::FourVector p = (*outit)->momentum();
	const double newPt = p.perp();
	if (abspid == 23) { // Z
	  if (newPt > pt) { // take the highest pT one 
	    newHaveSeen = abspid;
	    pt = newPt;
	  }
	} else if (abspid == 22) {
	  if (newPt > pt) { // take the highest pT one 
	    newHaveSeen = abspid;
	    pt = newPt;
	  }
	}
      } else if (haveSeen == 1000024) {
	if (abspid == 24) {
	  newHaveSeen = abspid;
	}
      } else if (haveSeen == 23 || haveSeen == 24) {
	if (abspid == 11 || abspid == 13 || abspid == 15) {
	  newHaveSeen = abspid;
	  m_leptons.push_back(*outit);
	  m_lightParticles.push_back(*outit);
	} else if ((abspid >= 1 && abspid <= 6) || abspid > 100) {
	  if (abspid < 6) {
	    m_lightParticles.push_back(*outit);
	  }
	  newHaveSeen = 1; // jet
	  // just to make it easier
	  abspid =1;
	}
      } else {
	if (abspid == 1000022 || abspid == 1000024) {
	  newHaveSeen = abspid;
	} else if (abspid == 1000021) {
	  m_isStrong = true;
	}
      }

      const HepMC::GenVertex *nextVertex = FindNextVertex(*outit);
      if (nextVertex) {
	decayVertices.push_back(nextVertex);
	pids.push_back(abspid);
      } else {
	decayVertices.push_back(NULL);
	pids.push_back(0);
      }
    }
  }

  // can now determine what needs to be sent down of the event type
  if (haveSeen == 1000022 && newHaveSeen == 22) {
    m_decays.push_back(gamma);
  } else if (haveSeen == 23) {
    switch (newHaveSeen) {
    case 1:
      m_decays.push_back(Zjj);
      break;
    case 11:
      m_decays.push_back(Zee);
      break;
    case 13:
      m_decays.push_back(Zmumu);
      break;
    case 15:
      m_decays.push_back(Ztautau);
      break;
    case 23:
      break;
    default:
      ATH_MSG_WARNING("Unexpected combination of Z haveSeen and newHaveSeen = " << newHaveSeen);
      break;
    }
  } else if (haveSeen == 24) {
    switch (newHaveSeen) {
    case 1:
      m_decays.push_back(Wjj);
      break;
    case 11:
      m_decays.push_back(Wenu);
      break;
    case 13:
      m_decays.push_back(Wmunu);
      break;
    case 15:
      if (!m_decayTaus) {
	m_decays.push_back(Wtaunu);
      }
      break;
    case 24:
      break;
    default:
      ATH_MSG_WARNING("Unexpected combination of W haveSeen and newHaveSeen = " << newHaveSeen);
      break;
    }
  }

  if (m_printDecayTree) msg(MSG::INFO) << endreq;
  for (int i = decayVertices.size(); i > 0; --i) {
    int index = i-1;
    if (decayVertices.at(index) != NULL) {
      if (m_printDecayTree) {
	msg(MSG::INFO) << "                 \t";
	for (int j = 0; j < index+extraSpaces; j++) {
	  msg(MSG::INFO) << "           ";
	}
      }
      // have to stop the case when you have seen both a gamma and a Z
      if (newHaveSeen == 22 && pids.at(index) == 23) {
	pids[index] = 22;
      }

      if (m_decayTaus && haveSeen == 24 && pids.at(index) == 15) {
	pids[index] = haveSeen; // ignore the tau from W decay
      }

      FollowDecayTree(decayVertices.at(index), index+extraSpaces, pids.at(index));
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
    //if (StatusGood((*outit)->status())) {
    if (1) {
      const int pid = (*outit)->pdg_id();
      if (m_printDecayTree) {
	HepMC::FourVector p = (*outit)->momentum();
	//msg(MSG::INFO) << std::setw(4) << std::right << round(p.perp()/GeV) << " ";
	msg(MSG::INFO) << std::setw(4) << std::right << (*outit)->status() << " ";
	const TParticlePDG* pdgpt = m_pdg.GetParticle(pid);
	if (pdgpt) {
	  msg(MSG::INFO) << std::setw(11) << std::left <<  pdgpt->GetName();
	} else {
	  msg(MSG::INFO) << std::setw(11) << std::left <<  pid;
	}
      }
      const HepMC::GenVertex *nextVertex = FindNextVertex(*outit);
      if (nextVertex) {
	decayVertices.push_back(nextVertex);
      }
    }
  }
  // msg(MSG::INFO) << "\t vertex = " << vtx->barcode();
  if (m_printDecayTree) msg(MSG::INFO) << endreq;
  for (int i = decayVertices.size(); i > 0; --i) {
    int index = i-1;
    if (decayVertices.at(index) != NULL) {
      if (m_printDecayTree) {
	msg(MSG::INFO) << "                 \t";
	for (int j = 0; j < index+extraSpaces; j++) {
	  msg(MSG::INFO) << "                ";
	}
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
      //(pid == 15) ||
      //(pid == 21) ||
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
*/
void TruthStudies::DumpEntireTree(const TruthParticleD3PDObject& truthObj) const
{
  ATH_MSG_INFO("***About to dump particles***");
  for (int i = 0; i < truthObj.n(); i++) {
    if (truthObj.barcode(i) < 200000) {
      msg(MSG::INFO) << i << ": ";
      const TParticlePDG* pdgpt = m_pdg.GetParticle(truthObj.pdgId(i));
      if (pdgpt) {
	msg(MSG::INFO) <<  pdgpt->GetName();
      } else {
	msg(MSG::INFO) <<  truthObj.pdgId(i);
      }
      
      msg(MSG::INFO) << " st: " << truthObj.status(i)
		     << " pt: " << truthObj.pt(i)
	// << " bar: " << truthObj.barcode(i)
		     << " par:";
      
      const std::vector<int>& parents = truthObj.parent_index(i);
      for (std::vector<int>::const_iterator it = parents.begin(); 
	   it != parents.end(); ++it) {
	msg(MSG::INFO) << " " << *it;
      }
      msg(MSG::INFO) << " chld: ";
      const std::vector<int>& children = truthObj.child_index(i);
      for (std::vector<int>::const_iterator it = children.begin(); 
	   it != children.end(); ++it) {
	msg(MSG::INFO) << " " << *it;
      }
      msg(MSG::INFO) << endmsg;
    }
  }
  ATH_MSG_INFO("***Finished dumping particles***");
}


void TruthStudies::FillEventType()
{

  if (m_decays.size() != 2) {
    ATH_MSG_DEBUG("Have the following unexpected number of decay types: " << m_decays.size());
    return;
  }

  decayType first, second;

  if (m_decays.at(0) <= m_decays.at(1)) {
    first = m_decays.at(0);
    second = m_decays.at(1);
  } else {
    first = m_decays.at(1);
    second = m_decays.at(0);
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
/*
int TruthStudies::findPhotons(const HepMC::GenEvent* genEvt)
{
  int NPhotons = 0;
  m_parentPids.clear();

  // first check basic requirements
  if (!passCuts()) {
    return 0;
  }

  for (HepMC::GenEvent::particle_const_iterator pitr = genEvt->particles_begin(); 
       pitr != genEvt->particles_end(); ++pitr) {
    if ((*pitr)->pdg_id() == 22 && (*pitr)->status()==1 && 
	(*pitr)->barcode() < 200000 && // barcode < 200K means it's not from GEANT 
 	(*pitr)->momentum().perp() >= m_Ptmin &&
	fabs((*pitr)->momentum().pseudoRapidity()) <= m_EtaRange) {
      ATH_MSG_DEBUG("Found a photon with pT = " << (*pitr)->momentum().perp() 
		    << ", eta = " << (*pitr)->momentum().pseudoRapidity()
		    << ", phi = " << (*pitr)->momentum().phi()
		    << ", status = " << (*pitr)->status()
		    << ", barcode = " << (*pitr)->barcode()
		    << ", and in particles:");
      const HepMC::GenParticle* parent = findParent(*pitr);
      if (parent) {
	const int pidParent = parent->pdg_id();
	ATH_MSG_DEBUG("  " << m_pdg.GetParticle(pidParent)->GetName() 
		      << ", pT = " << parent->momentum().perp() 
		      << ", eta = " << parent->momentum().pseudoRapidity()
		      << " with status = " << parent->status()
		      << " and barcode = " << parent->barcode());
	m_parentPids.push_back(pidParent);
	
	if (abs(pidParent) < 38) {
	  NPhotons++;
	}
      } else {
	// directly produced photon
	m_parentPids.push_back(0);
	NPhotons++;
      }
    }
  }
  return NPhotons;
}


int TruthStudies::findElectrons(const HepMC::GenEvent* genEvt)
{
  return 0;

  int NPhotons = 0;
  m_parentPids.clear();

  // // first check basic requirements
  // if (!passCuts()) {
  //   return 0;
  // }

  for (HepMC::GenEvent::particle_const_iterator pitr = genEvt->particles_begin(); 
       pitr != genEvt->particles_end(); ++pitr) {
    if (abs((*pitr)->pdg_id()) == 11 && (*pitr)->status()==1 && 
	(*pitr)->barcode() < 200000 ) // && // barcode < 200K means it's not from GEANT 
      //(*pitr)->momentum().perp() >= m_Ptmin &&
      //fabs((*pitr)->momentum().pseudoRapidity()) <= m_EtaRange) 
      {
      ATH_MSG_DEBUG("Found an electron with pT = " << (*pitr)->momentum().perp() 
		    << ", eta = " << (*pitr)->momentum().pseudoRapidity()
		    << ", phi = " << (*pitr)->momentum().phi()
		    << ", status = " << (*pitr)->status()
		    << ", barcode = " << (*pitr)->barcode()
		    << ", and in particles:");
      const HepMC::GenParticle* parent = findParent(*pitr);
      if (parent) {
	const int pidParent = parent->pdg_id();
	ATH_MSG_DEBUG("  " << m_pdg.GetParticle(pidParent)->GetName() 
		      << ", pT = " << parent->momentum().perp() 
		      << ", eta = " << parent->momentum().pseudoRapidity()
		      << " with status = " << parent->status()
		      << " and barcode = " << parent->barcode());
	m_parentPids.push_back(pidParent);
	
	if (abs(pidParent) < 38) {
	  NPhotons++;
	}
      } else {
	// directly produced photon
	m_parentPids.push_back(0);
	NPhotons++;
      }
    }
  }
  return NPhotons;
}

// returns the genparticle
const HepMC::GenParticle* TruthStudies::findParent(const HepMC::GenParticle* pcl) const
{
  ATH_MSG_DEBUG("findParent called with pcl = " << pcl);
  HepMC::GenVertex *prodVx = pcl->production_vertex();
  if (!prodVx) return 0;
  const int pinSize = prodVx->particles_in_size();
  if (pinSize != 1) {
    ATH_MSG_DEBUG("     Size of input particles = " << prodVx->particles_in_size() 
		  << ", if > 1 assuming hard scatter photon");
    return 0;
  }
  HepMC::GenVertex::particles_in_const_iterator pit = prodVx->particles_in_const_begin();
  if (pit != prodVx->particles_in_const_end()) {
    if ((*pit)->pdg_id() != 22) {
      return (*pit);
    } else {
      return findParent(*pit);
    }
  } else {
    return 0;
  }
}

bool TruthStudies::passCuts() const
{
  if (m_doMInv && m_lightParticles.size() >= 2) {
    for (std::vector<const HepMC::GenParticle *>::const_iterator it1 = m_lightParticles.begin();
	 it1 != m_lightParticles.end(); ++it1) {
      const CLHEP::HepLorentzVector lp1_p((*it1)->momentum().x(), 
					  (*it1)->momentum().y(),
					  (*it1)->momentum().z(), 
					  (*it1)->momentum().t());
      for (std::vector<const HepMC::GenParticle *>::const_iterator it2 = it1+1;
	   it2 != m_lightParticles.end(); ++it2) {
	const CLHEP::HepLorentzVector lp2_p((*it2)->momentum().x(), 
					    (*it2)->momentum().y(),
					    (*it2)->momentum().z(), 
					    (*it2)->momentum().t());
	if (lp1_p.invariantMass2(lp2_p) < m_mInv*m_mInv) {
	  return false;
	}
      }
    }
  }

  return true;
}


bool TruthStudies::passCuts(const HepMC::GenParticle* photon) const
{
  const CLHEP::HepLorentzVector ph_p(photon->momentum().x(), 
				     photon->momentum().y(),
				     photon->momentum().z(), 
				     photon->momentum().t());
  if (m_doDeltaRLepton) {
    for (std::vector<const HepMC::GenParticle *>::const_iterator it = m_leptons.begin();
	 it != m_leptons.end(); ++it) {
      const CLHEP::HepLorentzVector lep_p((*it)->momentum().x(), 
					  (*it)->momentum().y(),
					  (*it)->momentum().z(), 
					  (*it)->momentum().t());
      if (ph_p.deltaR(lep_p) < m_deltaRLepton) {
	return false;
      }
    }
  }
  if (m_doMInv) {
    for (std::vector<const HepMC::GenParticle *>::const_iterator it = m_lightParticles.begin();
	 it != m_lightParticles.end(); ++it) {
      const CLHEP::HepLorentzVector lp_p((*it)->momentum().x(), 
					 (*it)->momentum().y(),
					 (*it)->momentum().z(), 
					 (*it)->momentum().t());
      if (ph_p.invariantMass2(lp_p) < m_mInv*m_mInv) {
	return false;
      }
    }
  }

  return true;
}
*/
