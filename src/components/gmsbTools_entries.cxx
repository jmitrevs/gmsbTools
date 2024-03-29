#include "gmsbTools/gmsbPreparationTool.h"
#include "gmsbTools/gmsbSelectionTool.h"
#include "gmsbTools/gmsbOverlapCheckingTool.h"
#include "gmsbTools/gmsbOverlapRemovalTool.h"
#include "gmsbTools/TruthStudies.h"

#include "GaudiKernel/DeclareFactoryEntries.h"
 

DECLARE_TOOL_FACTORY( gmsbPreparationTool )
DECLARE_TOOL_FACTORY( gmsbSelectionTool )
DECLARE_TOOL_FACTORY( gmsbOverlapCheckingTool )
DECLARE_TOOL_FACTORY( gmsbOverlapRemovalTool )
DECLARE_TOOL_FACTORY( TruthStudies )

DECLARE_FACTORY_ENTRIES( gmsbTools )
{
    DECLARE_TOOL( UserMuonTool )

    DECLARE_TOOL( gmsbPreparationTool )
    DECLARE_TOOL( gmsbSelectionTool )
    DECLARE_TOOL( gmsbOverlapCheckingTool )
    DECLARE_TOOL( gmsbOverlapRemovalTool )
    DECLARE_ALGORITHM( TruthStudies )

}


