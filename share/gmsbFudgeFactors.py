
from gmsbTools.gmsbToolsConf import gmsbFudgeFactors
theGmsbFudgeFactors = gmsbFudgeFactors("gmsbFudgeFactors",
                                       # Master Collections
                                       PhotonInputName   = "PhotonAODCollection",
                                       # rerun the EMPIDBUILDER if requested
                                       PID_Builder = None,
                                       OutputLevel = DEBUG
                                       )
# EMPIDBuilder
from egammaTools.EMPIDBuilderBase import EMPIDBuilderBase
theaodempidbuilder=EMPIDBuilderBase("empidaod")
theaodempidbuilder.EMShowerContainerName="egDetailAOD"
theaodempidbuilder.EMTrackMatchContainerName="egDetailAOD"
theaodempidbuilder.EMConversionContainerName="egDetailAOD"
        
ToolSvc+=theaodempidbuilder
theGmsbFudgeFactors.PID_Builder = theaodempidbuilder

