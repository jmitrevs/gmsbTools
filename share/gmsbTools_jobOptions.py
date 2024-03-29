
NONE, EDIso, TrackIso, FlatTrackIso, LooseIso, MediumIso, TightIso, LooserIso = tuple(range(8))


from gmsbTools.gmsbToolsConf import \
    gmsbSelectionTool as ConfiguredUserSelectionTool
gmsbPreSelectionTool = ConfiguredUserSelectionTool(
    name = "gmsbPreSelectionTool",
    Simple = True,
    ElectronPt = -99*GeV,
    MuonPt = -99*GeV,
    PhotonPt = -99*GeV,
    JetPt = -99*GeV,
    #OutputLevel = DEBUG
)

ToolSvc += gmsbPreSelectionTool
print      gmsbPreSelectionTool

gmsbSelectionTool = ConfiguredUserSelectionTool(
    name = "gmsbSelectionTool",
    #OutputLevel = DEBUG
)

ToolSvc += gmsbSelectionTool
print      gmsbSelectionTool

gmsbFinalSelectionTool = ConfiguredUserSelectionTool(
    name = "gmsbFinalSelectionTool",
    DoElectronIsolation = LooseIso,
    DoEDPhotonIsolation = True,
    DoMuonIsolation = TightIso,
    ElectronPt = 20*GeV,
    MuonPt = 20*GeV,
    PhotonPt = 125*GeV,
    JetPt = 40*GeV,
    JetEta = 2.8,
    Simple = True,
    #OutputLevel = DEBUG
    )

ToolSvc += gmsbFinalSelectionTool
print      gmsbFinalSelectionTool


from gmsbTools.gmsbToolsConf import \
     gmsbPreparationTool as ConfiguredgmsbPreparationTool
gmsbPrePreparationTool = ConfiguredgmsbPreparationTool(
    name = "gmsbPrePreparationTool",

    # define the pre-selection tools
    UserSelectionTool = gmsbPreSelectionTool,
    
    # thelist of the input container keys - the order does not matter
    InputContainerKeys=[ "ph_",
                         "el_",
                         "mu_staco_",
                         "jet_AntiKt4LCTopo_"
                         ],
    
    
    
    # the list of the output container keys - these containers container the selected objects
    # The order matter::Should follow the same order as the input container keys above
    OutputContainerKeys=[ "pre_ph_",
                          "pre_el_",
                          "pre_mu_staco_",
                          "pre_jet_AntiKt4LCTopo_"
                          ],
    #OutputLevel = DEBUG
    
    )

ToolSvc += gmsbPrePreparationTool
print      gmsbPrePreparationTool

gmsbPreparationTool = ConfiguredgmsbPreparationTool(
    name = "gmsbPreparationTool",

    # define the pre-selection tools
    UserSelectionTool = gmsbSelectionTool,
    
    # thelist of the input container keys - the order does not matter
    InputContainerKeys=[ "pre_ph_",
                         "pre_el_",
                         "pre_mu_staco_",
                         "pre_jet_AntiKt4LCTopo_"
                         ],
    
    
    
    # the list of the output container keys - these containers container the selected objects
    # The order matter::Should follow the same order as the input container keys above
    OutputContainerKeys=[ "sl_ph_",
                          "sl_el_",
                          "sl_mu_staco_",
                          "sl_jet_AntiKt4LCTopo_"
                          ],
    #OutputLevel = DEBUG
    
    )

ToolSvc += gmsbPreparationTool
print      gmsbPreparationTool

from gmsbTools.gmsbToolsConf import \
     gmsbOverlapCheckingTool as ConfiguredgmsbOverlapCheckingTool
gmsbOverlapCheckingTool1 = ConfiguredgmsbOverlapCheckingTool(
    name = "gmsbOverlapCheckingTool1",
    OverlapDeltaRWithJets=0.2,
    #OutputLevel = DEBUG   
    )

ToolSvc += gmsbOverlapCheckingTool1
print      gmsbOverlapCheckingTool1

gmsbOverlapCheckingTool2 = ConfiguredgmsbOverlapCheckingTool(
    name = "gmsbOverlapCheckingTool2",
    OverlapDeltaRWithJets=0.4,
    #OutputLevel = DEBUG   
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
    # RemoveOverlapInSameContainer = True,
    
    # define the overlap checking tools
    UserOverlapCheckingTool = gmsbOverlapCheckingTool1,
    
    # thelist of the input container keys - the order is important: the overlap removing will be done in that order
    
    InputContainerKeys=[  "sl_el_",
                          "sl_ph_",
                          "sl_jet_AntiKt4LCTopo_"
                          ],
    
    
    IsAtlfastData=False, # set this to true if running on Atlfast AOD 
    # Only deltaR overlap removal is done for Atlfast, no cluster/TrackParticle overlap
    # and Cell/Hit overlap not done in the case of Atlfast 
    
    # the list of the output container keys - 
    OutputPhotonKey        = "int_ph_",
    OutputElectronKey      = "int_el_",
    OutputMuonKey          = "int_mu_staco_",
    OutputJetKey           = "int_jet_AntiKt4LCTopo_",
    # OutputBJetKey          = "IntermediateBJetCollection",
    # OutputLightJetKey      = "IntermediateLightJetCollection",
    
    #OutputLevel = DEBUG
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
    # RemoveOverlapInSameContainer = True,
    
    # define the overlap checking tools
    UserOverlapCheckingTool = gmsbOverlapCheckingTool2,
    
    # thelist of the input container keys - the order is important: the overlap removing will be done in that order
    
    InputContainerKeys=[  "int_jet_AntiKt4LCTopo_",
                          "int_el_",
                          "int_ph_",
                          "sl_mu_staco_"
                          ],
    
    
    IsAtlfastData=False, # set this to true if running on Atlfast AOD 
    # Only deltaR overlap removal is done for Atlfast, no cluster/TrackParticle overlap
    # and Cell/Hit overlap not done in the case of Atlfast 
    
    # the list of the output container keys - 
    OutputPhotonKey        = "fin_ph_",
    OutputElectronKey      = "fin_el_",
    OutputMuonKey          = "fin_mu_staco_",
    OutputJetKey           = "fin_jet_AntiKt4LCTopo_",
    # OutputBJetKey          = "FinalStateBJetCollection",
    # OutputLightJetKey      = "FinalStateLightJetCollection",
    
    #OutputLevel = DEBUG
    )

ToolSvc += gmsbOverlapRemovalTool2
print      gmsbOverlapRemovalTool2
