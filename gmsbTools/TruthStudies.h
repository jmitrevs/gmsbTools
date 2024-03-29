#ifndef GMSBTOOLS_TRUTHSTUDIES_H
#define GMSBTOOLS_TRUTHSTUDIES_H

#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/ToolHandle.h"

#include "TDatabasePDG.h"
#include <vector>

class TruthParticleD3PDObject;

/** Interface ID */  
static const InterfaceID IID_TruthStudies("TruthStudies", 1, 0);

/////////////////////////////////////////////////////////////////////////////
class TruthStudies:public AthAlgTool {
public:

  enum EventType {
    unknown = 0,
 
    Wenugamma = 1,
    Wmunugamma,
    Wtaunugamma,
    Wjjgamma,
    
    WenuWenu, // = 5
    WmunuWenu,
    WtaunuWenu,
    WjjWenu,
    
    WmunuWmunu, // = 9
    WtaunuWmunu,
    WjjWmunu,
    
    WtaunuWtaunu, // = 12
    WjjWtaunu,
    
    WjjWjj, // = 14
    
    WenuZee, // = 15
    WmunuZee,
    WtaunuZee,
    WjjZee,
    
    WenuZmumu, // = 19
    WmunuZmumu,
    WtaunuZmumu,
    WjjZmumu,
    
    WenuZtautau, // = 23
    WmunuZtautau,
    WtaunuZtautau,
    WjjZtautau,
    
    WenuZjj, // = 27
    WmunuZjj,
    WtaunuZjj,
    WjjZjj,
    
    Zeegamma, // 31
    Zmumugamma,
    Ztautaugamma,
    Zjjgamma,

    ZeeZee, // = 35
    ZmumuZee,
    ZtautauZee,
    ZjjZee,
    
    ZmumuZmumu, // = 39
    ZtautauZmumu,
    ZjjZmumu,
    
    ZtautauZtautau, // = 42
    ZjjZtautau,
    
    ZjjZjj, // = 44

    gammagamma, // 45

    numEventTypes
  };


  TruthStudies (const std::string& type,
		const std::string& name, 
		const IInterface* parent);

  /** AlgTool and IAlgTool interface methods */
  static const InterfaceID& interfaceID() { return IID_TruthStudies; };

  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

  EventType GetEventType() const { return m_type; };

  bool isStrong() const {return m_isStrong; };

  int nPhotons() const {return m_nPhotons; };

  double Wpt() const { return m_Wpt; };

private:

  // this is for one side of the decay
  enum decayType {
    gamma,
    Zee,
    Zmumu,
    Ztautau,
    Zjj,
    Wenu,
    Wmunu,
    Wtaunu,
    Wjj
  };

  void FillEventType();

  /*
  HepMC::GenVertex* getMCHardInteraction(const HepMC::GenEvent *const ge) const;

  // haveSeen is the pdgid of what has been seen 
  // also determines the event type
  void FollowDecayTree(const HepMC::GenVertex *vtx, int extraSpaces=0, int haveSeen = 0);

   // also adds the pT of each particle. (Doesn't deterimine the event type)
  void FollowDecayTreeAnnotated(const HepMC::GenVertex *vtx, int extraSpaces=0) const;
  */
  void DumpEntireTree(const TruthParticleD3PDObject& truthObj) const;

  // some utilities for it
  bool StatusGood(int status) const;
  
  /*
  const HepMC::GenVertex *FindNextVertex(const HepMC::GenParticle *pcl) const;

  int findPhotons(const TruthParticleD3PDObject& truthObj);
  int findElectrons(const TruthParticleD3PDObject& truthObj);

  // returns the PID
  //const HepMC::GenParticle* findParent(const HepMC::GenParticle* pcl) const;

  bool passCuts() const;
  bool passCuts(int photon) const;
  */

  // a PDG database that can be used to get particle properties
  TDatabasePDG m_pdg;

  /** name of the AOD truth particle container to retrieve from StoreGate */
  std::string m_truthParticleContainerName;

  // /** can alternately use the McEventCollection */
  // std::string m_mcEventCollectionName;


  bool m_printDecayTree;
  bool m_useAnnotated;
  bool m_dumpEntireTree;

  double m_Ptmin;
  double m_EtaRange;

  // the overall event type
  EventType m_type;

  bool m_isStrong;
  
  int m_nPhotons;
  std::vector<int> m_parentPids;

  // where we classify the decay types
  std::vector<decayType> m_decays;

  double m_Wpt; // for reweighing


  bool m_doDeltaRLepton;
  double m_deltaRLepton;
  bool m_doDeltaRLight;
  double m_deltaRLight;
  bool m_doMInv;
  double m_mInv;
  bool m_decayTaus; // only for classification, not for results

  int m_WptID; // the W ID (usually 24, but can change to 23 for Z)

  // // to match truth filters
  std::vector<int> m_leptons;
  std::vector<int> m_lightParticles;
  

};

// inline bool TruthStudies::StatusGood(int status) const 
// {
//   return (status == 1 || status == 3);
// } 

inline bool TruthStudies::StatusGood(int status) const 
{
  return (status != 141 && status != 142);
} 


#endif // GMSBTOOLS_TRUTHSTUDIES_H
