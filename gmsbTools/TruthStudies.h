#ifndef GMSBTOOLS_TRUTHSTUDIES_H
#define GMSBTOOLS_TRUTHSTUDIES_H

#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/ToolHandle.h"

#include "TDatabasePDG.h"

class Jet;
namespace Reco  { class ITrackToVertex; }

namespace HepMC {
  class GenVertex;
  class GenParticle;
  class GenEvent;
  class FourVector;
}


/** Interface ID */  
static const InterfaceID IID_TruthStudies("TruthStudies", 1, 0);

/////////////////////////////////////////////////////////////////////////////
class TruthStudies:public AthAlgTool {
public:

  enum EventType {
    unknown = -999,
 
    Wenugamma = 0,
    Wmunugamma,
    Wtaunugamma,
    Wjjgamma,
    
    WenuWenu, // = 4
    WmunuWenu,
    WtaunuWenu,
    WjjWenu,
    
    WmunuWmunu, // = 8
    WtaunuWmunu,
    WjjWmunu,
    
    WtaunuWtaunu, // = 11
    WjjWtaunu,
    
    WjjWjj, // = 13
    
    WenuZee, // = 14
    WmunuZee,
    WtaunuZee,
    WjjZee,
    
    WenuZmumu, // = 18
    WmunuZmumu,
    WtaunuZmumu,
    WjjZmumu,
    
    WenuZtautau, // = 22
    WmunuZtautau,
    WtaunuZtautau,
    WjjZtautau,
    
    WenuZjj, // = 26
    WmunuZjj,
    WtaunuZjj,
    WjjZjj,
    
    Zeegamma, // 30
    Zmumugamma,
    Ztautaugamma,
    Zjjgamma,

    ZeeZee, // = 34
    ZmumuZee,
    ZtautauZee,
    ZjjZee,
    
    ZmumuZmumu, // = 38
    ZtautauZmumu,
    ZjjZmumu,
    
    ZtautauZtautau, // = 41
    ZjjZtautau,
    
    ZjjZjj, // = 43

    gammagamma, // 44

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

  void FillEventType(decayType d1, decayType d2);

  HepMC::GenVertex* getMCHardInteraction(const HepMC::GenEvent *const ge) const;

  void FollowDecayTree(const HepMC::GenVertex *vtx, int extraSpaces=0);

   // also adds the pT of each particle
  void FollowDecayTreeAnnotated(const HepMC::GenVertex *vtx, int extraSpaces=0) const;

  void DumpEntireTree(const HepMC::GenEvent *ge) const;

  // some utilities for it
  bool StatusGood(int status) const;
  const HepMC::GenVertex *FindNextVertex(const HepMC::GenParticle *pcl) const;


  // a PDG database that can be used to get particle properties
  TDatabasePDG m_pdg;

  /** name of the AOD truth particle container to retrieve from StoreGate */
  std::string m_truthParticleContainerName;

  /** can alternately use the McEventCollection */
  std::string m_mcEventCollectionName;


  bool m_printDecayTree;
  bool m_useAnnotated;
  bool m_dumpEntireTree;

  EventType m_type;

};

// inline bool TruthStudies::StatusGood(int status) const 
// {
//   return (status == 1 || status == 3);
// } 

inline bool TruthStudies::StatusGood(int status) const 
{
  return (status != 120 && status != 141 && status != 142);
} 


#endif // GMSBTOOLS_TRUTHSTUDIES_H
