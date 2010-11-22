
from UserAnalysisUtils.UserAnalysisUtilsConf import \
    UserAnalysisSelectionTool as ConfiguredUserSelectionTool
UserAnalysisSelectionTool = ConfiguredUserSelectionTool(name = "UserAnalysisSelectionTool")

ToolSvc += UserAnalysisSelectionTool
print      UserAnalysisSelectionTool

UserAnalysisCrackSelectionTool = ConfiguredUserSelectionTool(
    name = "UserAnalysisCrackSelectionTool",
    ElectronPt = 10.0*GeV,
    PhotonPt = 10.0*GeV,
    DoElectronIsolation = False,
    DoPhotonIsolation = False,
    PhotonEtaWindowMin = 0,
    ElectronEtaWindowMin = 0,
    PhotonEtaWindowMax = 1.37,
    ElectronEtaWindowMax = 1.37,
    ElectronEta = 1.52,
    PhotonEta = 1.52,
    )

ToolSvc += UserAnalysisCrackSelectionTool
print      UserAnalysisCrackSelectionTool


from UserAnalysisUtils.UserAnalysisUtilsConf import \
     UserAnalysisPreparationTool as ConfiguredUserAnalysisPreparationTool
UserAnalysisPreparationTool = ConfiguredUserAnalysisPreparationTool(
    name = "UserAnalysisPreparationTool",

    # define the pre-selection tools
    UserSelectionTool = UserAnalysisSelectionTool,
    
    # thelist of the input container keys - the order does not matter
    InputContainerKeys=[ "PhotonAODCollection",
                         "ElectronAODCollection",
                         "AntiKt4TopoJets"
                         ],
    
    
    
    # the list of the output container keys - these containers container the selected objects
    # The order matter::Should follow the same order as the input container keys above
    OutputContainerKeys=[ "SelectedPhotonAODCollection",
                          "SelectedElectronAODCollection",
                          "SelectedAntiKt4TopoJets"
                          ]
    
    )

ToolSvc += UserAnalysisPreparationTool
print      UserAnalysisPreparationTool

UserAnalysisCrackPreparationTool = ConfiguredUserAnalysisPreparationTool(
    name = "UserAnalysisCrackPreparationTool",
    
    # define the pre-selection tools
    UserSelectionTool = UserAnalysisCrackSelectionTool,
    
    # thelist of the input container keys - the order does not matter
    InputContainerKeys=[ "PhotonAODCollection",
                         "ElectronAODCollection"
                         ],



    # the list of the output container keys - these containers container the selected objects
    # The order matter::Should follow the same order as the input container keys above
    OutputContainerKeys=[ "CrackPhotonCollection",
                          "CrackElectronCollection"
                          ]
    
  )

ToolSvc += UserAnalysisCrackPreparationTool
print      UserAnalysisCrackPreparationTool

from UserAnalysisUtils.UserAnalysisUtilsConf import \
     UserAnalysisOverlapCheckingTool as ConfiguredUserAnalysisOverlapCheckingTool
UserAnalysisOverlapCheckingTool1 = ConfiguredUserAnalysisOverlapCheckingTool(
    name = "UserAnalysisOverlapCheckingTool1",
    OverlapDeltaRWithJets=0.2 
    )

ToolSvc += UserAnalysisOverlapCheckingTool1
print      UserAnalysisOverlapCheckingTool1

UserAnalysisOverlapCheckingTool2 = ConfiguredUserAnalysisOverlapCheckingTool(
    name = "UserAnalysisOverlapCheckingTool2",
    OverlapDeltaRWithJets=0.4
    )

ToolSvc += UserAnalysisOverlapCheckingTool2
print      UserAnalysisOverlapCheckingTool2


from UserAnalysisUtils.UserAnalysisUtilsConf import \
     UserAnalysisOverlapRemovalTool as ConfiguredUserAnalysisOverlapRemovalTool
UserAnalysisOverlapRemovalTool1 = ConfiguredUserAnalysisOverlapRemovalTool(
    name = "UserAnalysisOverlapRemovalTool1",

  # define the pre-selection tools - used here only to check if a jet is a b-jet
  UserSelectionTool = UserAnalysisSelectionTool,

  # Whether to check overlap in same container or not. 
  # For example, muon overlapping with muon?
  # Currently when set to False, it applies to all contianers. 
  RemoveOverlapInSameContainer = True,

  # define the overlap checking tools
  UserOverlapCheckingTool = UserAnalysisOverlapCheckingTool1,

  # thelist of the input container keys - the order is important: the overlap removing will be done in that order
  
  InputContainerKeys=[  "SelectedPhotonAODCollection",
                        "SelectedElectronAODCollection",
                        "SelectedAntiKt4TopoJets"
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
  OutputTauJetKey        = "IntermediateTauJetCollection",
  OutputCalloClusterKey  = "IntermediateCaloClusterCollection",
  OutputTrackParticleKey = "IntermediateTrackParticleCollection",
  OutputJetKey           = "IntermediateJetCollection",
  OutputBJetKey          = "IntermediateBJetCollection",
  OutputLightJetKey      = "IntermediateLightJetCollection"

    )

ToolSvc += UserAnalysisOverlapRemovalTool1
print      UserAnalysisOverlapRemovalTool1

UserAnalysisOverlapRemovalTool2 = ConfiguredUserAnalysisOverlapRemovalTool(
    name = "UserAnalysisOverlapRemovalTool2",

    # define the pre-selection tools - used here only to check if a jet is a b-jet
    UserSelectionTool = UserAnalysisSelectionTool,
    
    # Whether to check overlap in same container or not. 
    # For example, muon overlapping with muon?
    # Currently when set to False, it applies to all contianers. 
    RemoveOverlapInSameContainer = True,
    
    # define the overlap checking tools
    UserOverlapCheckingTool = UserAnalysisOverlapCheckingTool2,
    
    # thelist of the input container keys - the order is important: the overlap removing will be done in that order
    
    InputContainerKeys=[  "IntermediateJetCollectioin"
                          "IntermediatePhotonCollection",
                          "IntermediateElectronCollection",
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
    OutputTauJetKey        = "FinalStateTauJetCollection",
    OutputCalloClusterKey  = "FinalStateCaloClusterCollection",
    OutputTrackParticleKey = "FinalStateTrackParticleCollection",
    OutputJetKey           = "FinalStateJetCollection",
    OutputBJetKey          = "FinalStateBJetCollection",
    OutputLightJetKey      = "FinalStateLightJetCollection"
    
    )

ToolSvc += UserAnalysisOverlapRemovalTool2
print      UserAnalysisOverlapRemovalTool2
