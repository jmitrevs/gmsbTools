#ifndef GMSBTOOLS_IJETCLEANINGTOOL_H
#define GMSBTOOLS_IJETCLEANINGTOOL_H

#include "GaudiKernel/IAlgTool.h"

class JetCollection;

static const InterfaceID IID_IJetCleaningTool("IJetCleaningTool", 1, 0);

class IJetCleaningTool : public virtual IAlgTool {
 public:

  
  // Retrieve Interface ID
  static const InterfaceID& interfaceID() { return IID_IJetCleaningTool; }

  virtual ~IJetCleaningTool(){};

  virtual StatusCode initialize()=0;
  virtual StatusCode finalize()=0;
  
  virtual bool passCleaningCuts(const JetCollection*)=0;
  virtual bool passCleaningCuts(const JetCollection*,bool)=0;
       
};

#endif
