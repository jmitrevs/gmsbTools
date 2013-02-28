#ifndef GMSBTOOLS_GMSBSYSTERROR_H  
#define GMSBTOOLS_GMSBSYSTERROR_H 

namespace SystErr
{
  typedef enum  {
    NONE,
    JESDOWN,JESUP, // Return the total JES uncertainty as the quadratic sum of the different uncertainties listed below 
    EffectiveNP_1_Down,EffectiveNP_1_Up,
    EffectiveNP_2_Down,EffectiveNP_2_Up,
    EffectiveNP_3_Down,EffectiveNP_3_Up,
    EffectiveNP_4_Down,EffectiveNP_4_Up,
    EffectiveNP_5_Down,EffectiveNP_5_Up,
    EffectiveNP_6_Down,EffectiveNP_6_Up,
    EtaIntercalibration_Modelling_Down,EtaIntercalibration_Modelling_Up, 
    EtaIntercalibration_StatAndMethod_Down,EtaIntercalibration_StatAndMethod_Up, 
    SingleParticle_HighPt_Down,SingleParticle_HighPt_Up,
    RelativeNonClosure_Pythia8_Down,RelativeNonClosure_Pythia8_Up,
    PileupOffsetTermMuDown,PileupOffsetTermMuUp, 
    PileupOffsetTermNPVDown,PileupOffsetTermNPVUp, 
    PileupPtTermDown,PileupPtTermUp,
    PileupRhoTopologyDown,PileupRhoTopologyUp, 
    CloseByDown,CloseByUp,
    FlavorCompUncertDown,FlavorCompUncertUp, 
    FlavorResponseUncertDown,FlavorResponseUncertUp, 
    BJesDown,BJesUp, //Break down of the JES uncertainties  
    JER, 
    EGZEEUP, EGZEEDOWN,EGMATUP,EGMATDOWN,EGPSUP,EGPSDOWN,EGLOWUP,EGLOWDOWN, EGRESDOWN, EGRESUP, EEFFDOWN, EEFFUP, ETRIGDOWN, ETRIGUP,
    MMSLOW, MMSUP, MIDLOW, MIDUP, MSCALELOW, MSCALEUP, MEFFDOWN, MEFFUP,MTRIGDOWN,MTRIGUP,
    BJETDOWN, BJETUP, CJETDOWN, CJETUP, BMISTAGDOWN, BMISTAGUP,
    SCALESTUP,SCALESTDOWN,RESOST,
    TESUP, TESDOWN /// Tau energy scale and resolution
    //    RESOSTUP,RESOSTDOWN, // Just have one RESOST, because the UP and DOWN variations are both symmetric smearings.
    //    SCALEPHUP,SCALEPHDOWN,RESOPHUP,RESOPHDOWN,RESOPHUPDOWN,RESOPHDOWNUP,
  } Syste;
}

#endif
