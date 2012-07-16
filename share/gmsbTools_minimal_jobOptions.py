
# from PhotonAnalysisUtils.LowPtJetFinder import LowPtJetFinder
# mygetter2 = LowPtJetFinder() # will insert the alg into the AlgSequence()
# print mygetter2

# from egammaCaloTools.egammaIsoPtCorrectionToolBase import egammaIsoPtCorrectionToolBase
# theegammaisoptcorrection=egammaIsoPtCorrectionToolBase("egammaIsoPtCorrection")
# ToolSvc+=theegammaisoptcorrection

# from PhotonAnalysisUtils.PhotonAnalysisUtilsConf import PAUcaloIsolationTool
# mycaloisolationtool = PAUcaloIsolationTool(DoAreaCorrections = True,
#                                            EMCaloIsoPtCorrectionTool = theegammaisoptcorrection,
#                                            # OutputLevel = DEBUG
#                                            )
# ToolSvc += mycaloisolationtool

from gmsbTools.gmsbToolsConf import \
    gmsbSelectionTool as ConfiguredUserSelectionTool
gmsbSelectionTool = ConfiguredUserSelectionTool(
    name = "gmsbSelectionTool",
    PAUcaloIsolationTool = None,
    DoEDPhotonIsolation = False,
    #OutputLevel = DEBUG
)

ToolSvc += gmsbSelectionTool
print      gmsbSelectionTool

gmsbFinalSelectionTool = ConfiguredUserSelectionTool(
    name = "gmsbFinalSelectionTool",
    DoElectronTrackIsolation = True,
    DoEDPhotonIsolation = False, # already done
    DoMuonIsoCut = True,
    MuonPt = 25*GeV,
    Simple = True
    )

ToolSvc += gmsbFinalSelectionTool
print      gmsbFinalSelectionTool


from gmsbTools.gmsbToolsConf import \
     gmsbPreparationTool as ConfiguredgmsbPreparationTool
gmsbPreparationTool = ConfiguredgmsbPreparationTool(
    name = "gmsbPreparationTool",

    # define the pre-selection tools
    UserSelectionTool = gmsbSelectionTool,
    
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

ToolSvc += gmsbPreparationTool
print      gmsbPreparationTool

from gmsbTools.gmsbToolsConf import \
     gmsbOverlapCheckingTool as ConfiguredgmsbOverlapCheckingTool
gmsbOverlapCheckingTool1 = ConfiguredgmsbOverlapCheckingTool(
    name = "gmsbOverlapCheckingTool1",
    OverlapDeltaRWithJets=0.2
    )

ToolSvc += gmsbOverlapCheckingTool1
print      gmsbOverlapCheckingTool1

gmsbOverlapCheckingTool2 = ConfiguredgmsbOverlapCheckingTool(
    name = "gmsbOverlapCheckingTool2",
    OverlapDeltaRWithJets=0.4
    )

ToolSvc += gmsbOverlapCheckingTool2
print      gmsbOverlapCheckingTool2


from gmsbTools.gmsbToolsConf import \
     gmsbOverlapRemovalTool as ConfiguredgmsbOverlapRemovalTool
gmsbOverlapRemovalTool1 = ConfiguredgmsbOverlapRemovalTool(
    name = "gmsbOverlapRemovalTool1",

    # define the pre-selection tools - used here only to check if a jet is a b-jet
    UserSelectionTool = gmsbSelectionTool,
    
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

ToolSvc += gmsbOverlapRemovalTool1
print      gmsbOverlapRemovalTool1

gmsbOverlapRemovalTool2 = ConfiguredgmsbOverlapRemovalTool(
    name = "gmsbOverlapRemovalTool2",

    # define the pre-selection tools - used here only to check if a jet is a b-jet
    UserSelectionTool = gmsbSelectionTool,
    
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

ToolSvc += gmsbOverlapRemovalTool2
print      gmsbOverlapRemovalTool2
