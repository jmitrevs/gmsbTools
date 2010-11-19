#ifndef GMSBTOOLS_JETCLEANINGTOOL_H
#define GMSBTOOLS_JETCLEANINGTOOL_H

#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ITHistSvc.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "StoreGate/StoreGateSvc.h"
#include "gmsbTools/IJetCleaningTool.h"

class JetCollection;
class Jet;

typedef enum  { LooseBad, TightBad  } BadJetCategory;

class JetCleaningTool : public AthAlgTool, public virtual IJetCleaningTool
{
 public:
  
  //default constructor due to Athena interface
  JetCleaningTool(const std::string& t, const std::string& n, const IInterface*  p);
  
  //destructor
  virtual ~JetCleaningTool();
  virtual StatusCode initialize();
  virtual StatusCode finalize();

  bool passCleaningCuts(const JetCollection*);
  bool passCleaningCuts(const JetCollection*,bool);
  
 private:

  bool passJetPtCut(const Jet*) const;
  bool isBadJet(const Jet*,bool) const;  
  bool failsJetCuts(BadJetCategory criteria, double quality, double n90,
		    double emf, double hecf, double time,
		    double fmax, double eta) const;
 
  //      configurable data members
  //------------------------------------------------------------------------
  float m_jet_pt_cut;
  bool m_use_emscale;
  bool m_use_JES;
 
};

#endif
