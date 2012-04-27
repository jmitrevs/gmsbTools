
from gmsbTools.gmsbToolsConf import \
    gmsbSelectionTool as ConfiguredUserSelectionTool
gmsbLooseSelectionTool = ConfiguredUserSelectionTool(
    name = "gmsbLooseSelectionTool",
    ElectronID = egammaPID.ElectronIDLoosePP,
    # OutputLevel = DEBUG
)

ToolSvc += gmsbLooseSelectionTool
print      gmsbLooseSelectionTool

gmsbLooseFinalSelectionTool = ConfiguredUserSelectionTool(
    name = "gmsbLooseFinalSelectionTool",
    DoNewElectronIsolation = False,
    DoNewPhotonIsolation = False,
    DoElectronEtaWindowCut = True,
    DoPhotonEtaWindowCut = True,
    DoMuonIsoCut = True,
    MuonPt = 25*GeV,
    Simple = True
    )

ToolSvc += gmsbLooseFinalSelectionTool
print      gmsbLooseFinalSelectionTool


from gmsbTools.gmsbToolsConf import \
     gmsbPreparationTool as ConfiguredgmsbPreparationTool
gmsbLoosePreparationTool = ConfiguredgmsbPreparationTool(
    name = "gmsbLoosePreparationTool",

    # define the pre-selection tools
    UserSelectionTool = gmsbLooseSelectionTool,
    
    # thelist of the input container keys - the order does not matter
    InputContainerKeys=[ "PhotonAODCollection",
                         "ElectronAODCollection",
                         "StacoMuonCollection",
                         "AntiKt4TopoEMJets"
                         ],
    
    
    
    # the list of the output container keys - these containers container the selected objects
    # The order matter::Should follow the same order as the input container keys above
    OutputContainerKeys=[ "SelectedPhotonAODCollection",
                          "SelectedElectronAODCollection",
                          "SelectedMuonCollection",
                          "SelectedAntiKt4TopoEMJets"
                          ],
    #OutputLevel = DEBUG
    
    )

ToolSvc += gmsbLoosePreparationTool
print      gmsbLoosePreparationTool

from gmsbTools.gmsbToolsConf import \
     gmsbOverlapRemovalTool as ConfiguredgmsbOverlapRemovalTool
gmsbLooseOverlapRemovalTool1 = ConfiguredgmsbOverlapRemovalTool(
    name = "gmsbLooseOverlapRemovalTool1",

    # define the pre-selection tools - used here only to check if a jet is a b-jet
    UserSelectionTool = gmsbLooseSelectionTool,
    
    # Whether to check overlap in same container or not. 
    # For example, muon overlapping with muon?
    # Currently when set to False, it applies to all contianers. 
    RemoveOverlapInSameContainer = True,
    
    # define the overlap checking tools
    UserOverlapCheckingTool = gmsbOverlapCheckingTool1,
    
    # thelist of the input container keys - the order is important: the overlap removing will be done in that order
    
    InputContainerKeys=[  "SelectedElectronAODCollection",
                          "SelectedPhotonAODCollection",
                          "SelectedAntiKt4TopoEMJets"
                          ],
    
    
    IsAtlfastData=False, # set this to true if running on Atlfast AOD 
    # Only deltaR overlap removal is done for Atlfast, no cluster/TrackParticle overlap
    # and Cell/Hit overlap not done in the case of Atlfast 
    
    # the list of the output container keys - 
    OuputObjectKey         = "IntermediateObjectCollection",
    OutputLeptonKey        = "IntermediateLeptonCollection",
    OutputPhotonKey        = "IntermediatePhotonCollection",
    OutputElectronKey      = "IntermediateElectronCollection",
    OutputMuonKey          = "IntermediateMuonCollection",
    OutputTauJetKey        = "IntermediateTauCollection",
    OutputCalloClusterKey  = "IntermediateCaloClusterCollection",
    OutputTrackParticleKey = "IntermediateTrackParticleCollection",
    OutputJetKey           = "IntermediateJetCollection",
    OutputBJetKey          = "IntermediateBJetCollection",
    OutputLightJetKey      = "IntermediateLightJetCollection",
    
    # OutputLevel = DEBUG
    )

ToolSvc += gmsbLooseOverlapRemovalTool1
print      gmsbLooseOverlapRemovalTool1

gmsbLooseOverlapRemovalTool2 = ConfiguredgmsbOverlapRemovalTool(
    name = "gmsbLooseOverlapRemovalTool2",

    # define the pre-selection tools - used here only to check if a jet is a b-jet
    UserSelectionTool = gmsbLooseSelectionTool,
    
    # Whether to check overlap in same container or not. 
    # For example, muon overlapping with muon?
    # Currently when set to False, it applies to all contianers. 
    RemoveOverlapInSameContainer = True,
    
    # define the overlap checking tools
    UserOverlapCheckingTool = gmsbOverlapCheckingTool2,
    
    # thelist of the input container keys - the order is important: the overlap removing will be done in that order
    
    InputContainerKeys=[  "IntermediateJetCollection",
                          "IntermediatePhotonCollection",
                          "IntermediateElectronCollection",
                          "SelectedMuonCollection"
                          ],
    
    
    IsAtlfastData=False, # set this to true if running on Atlfast AOD 
    # Only deltaR overlap removal is done for Atlfast, no cluster/TrackParticle overlap
    # and Cell/Hit overlap not done in the case of Atlfast 
    
    # the list of the output container keys - 
    OuputObjectKey         = "FinalStateObjectCollection",
    OutputLeptonKey        = "FinalStateLeptonCollection",
    OutputPhotonKey        = "FinalStatePhotonCollection",
    OutputElectronKey      = "FinalStateElectronCollection",
    OutputMuonKey          = "FinalStateMuonCollection",
    OutputTauJetKey        = "FinalStateTauCollection",
    OutputCalloClusterKey  = "FinalStateCaloClusterCollection",
    OutputTrackParticleKey = "FinalStateTrackParticleCollection",
    OutputJetKey           = "FinalStateJetCollection",
    OutputBJetKey          = "FinalStateBJetCollection",
    OutputLightJetKey      = "FinalStateLightJetCollection",
    
    # OutputLevel = DEBUG
    )

ToolSvc += gmsbLooseOverlapRemovalTool2
print      gmsbLooseOverlapRemovalTool2
