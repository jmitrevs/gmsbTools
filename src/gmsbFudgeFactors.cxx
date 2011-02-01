/********************************************************************

NAME:     gmsbFudgeFactors.cxx, based on egammaAODRender

AUTHORS:  J. Mitrevski
CREATED:  

PURPOSE:  Apply fudge factors and recalculate isEM

********************************************************************/

// INCLUDE HEADER FILES:
#include "gmsbTools/gmsbFudgeFactors.h"
#include "egammaInterfaces/IegammaBaseTool.h"
#include "egammaInterfaces/IEMPIDBuilder.h"
#include "egammaInterfaces/IEMClusterTool.h"

#include "egammaEvent/Photon.h"
#include "egammaEvent/PhotonContainer.h"
#include "egammaEvent/EMShower.h"

#include <algorithm> 
#include <math.h>

//  END OF HEADER FILES INCLUDE

/////////////////////////////////////////////////////////////////
//  CONSTRUCTOR:
gmsbFudgeFactors::gmsbFudgeFactors(const std::string& name, 
				 ISvcLocator* pSvcLocator): 
  AthAlgorithm(name, pSvcLocator),
  m_pidBuilder("EMPIDBuilder/empidaod")
{
	
  // The following properties are specified at run-time
  // (declared in jobOptions file)
	
  // Name of the input photon collection
  declareProperty("PhotonInputName",   
		  m_PhotonInputName  = "PhotonAODCollection",
		  "Name of the input photon collection");

  // Boolean to apply PID tools
  declareProperty("doEMPID",         
		  m_doEMPID = true,
		  "Boolean to apply PID tools");

  // Boolean to apply PID tools
  declareProperty("doFudgeFactors",         
		  m_doFudgeFactors = true,
		  "Boolean to apply fudge factors");

  // Name of the PID tool
  declareProperty("PID_Builder", 
		  m_pidBuilder,
		  "Handle of the PID tool");

  // dump the result (not yet implemented)
  declareProperty("Dump",
		  m_dump=false,
		  "Boolean to dump the result");

}

// DESTRUCTOR:
gmsbFudgeFactors::~gmsbFudgeFactors()
{  
}

/////////////////////////////////////////////////////////////////
// INITIALIZE METHOD:
StatusCode gmsbFudgeFactors::initialize() 
{

  ATH_MSG_INFO("Initializing gmsbFudgeFactors");
 
  // retrieve PID builder:
  if (m_doEMPID) {

    // a priori this is not useful
    if(m_pidBuilder.retrieve().isFailure()) {
      ATH_MSG_ERROR("Unable to retrieve "<<m_pidBuilder);
      return StatusCode::FAILURE;
    } 
    else ATH_MSG_DEBUG("Retrieved Tool " << m_pidBuilder); 
  }

  return StatusCode::SUCCESS;
}

// FINALIZE METHOD:
StatusCode gmsbFudgeFactors::finalize()
{
  return StatusCode::SUCCESS;
}

/////////////////////////////////////////////////////////////////
// ATHENA EXECUTE METHOD:
StatusCode gmsbFudgeFactors::execute()
{  
	
  ATH_MSG_DEBUG("Executing gmsbFudgeFactors");
	
  StatusCode sc;
  const PhotonContainer* photons;
  sc = evtStore()->retrieve(photons, m_PhotonInputName);
	
  if(sc.isFailure()  ||  !photons) {
    ATH_MSG_ERROR("no PhotonInputContainer found in TDS");
    return sc;
  }
  

  for(PhotonContainer::const_iterator ph = photons->begin(); 
      ph != photons->end();
      ph++) {

    if (m_doFudgeFactors) {
      const EMShower *constShower = (*ph)->detail<EMShower>();
      EMShower *shower = const_cast<EMShower*>(constShower);
      if (!shower) {
	ATH_MSG_ERROR("Could not find the shower detail");
	return StatusCode::FAILURE;
      }

      const double eta2 = fabs((*ph)->cluster()->etaBE(2));
      const int ieta = getEtaBin(eta2);
      const int ipt = getPtBin((*ph)->pt());

      // the hadronic ratios
      double et = 0.;
      if (fabs(eta2)<999.) 
	et = cosh(eta2)!=0. ? (*ph)->cluster()->energy()/cosh(eta2) : 0.;

      const double raphad1 = fabs(et)>0. ? shower->ethad1()/et : 0.;
      const double raphad  = fabs(et)>0. ? shower->ethad()/et : 0.;

      const double raphad1_fix = raphad1 + rhad1_FF[ieta][ipt];
      const double raphad_fix = raphad + rhad_FF[ieta][ipt];

      shower->set_ethad(raphad_fix * et);
      shower->set_ethad1(raphad1_fix * et);

      // reta, rphi

      const double e233 = shower->e233(); 
      const double e237 = shower->e237(); 
      const double e277 = shower->e277(); 

      const double Reta37 = fabs(e277)>0. ? e237/e277 : 0.;
      const double Rphi33 = fabs(e237)>0. ? e233/e237 : 0.;
     
      const double Reta37_fix = Reta37 + reta_FF[ieta][ipt];
      const double Rphi33_fix = Rphi33 + rphi_FF[ieta][ipt];
      
      shower->set_e233(Rphi33_fix * e237);
      if (Reta37_fix != 0) {
	shower->set_e277(1.0/Rphi33_fix * e237);
      }

      // weta2
      shower->set_weta2(shower->weta2() + weta2_FF[ieta][ipt]);

      // w1
      shower->set_weta1(shower->weta1() + w1_FF[ieta][ipt]);
      
      // wtot
      shower->set_wtots1(shower->wtots1() + wtot_FF[ieta][ipt]);

      // fracm
      shower->set_fracs1(shower->fracs1() + fracm_FF[ieta][ipt]);

      // deltae and demaxs1
      // for detae
      shower->set_emins1(shower->emins1() - deltae_FF[ieta][ipt]); // note -

      //  E of 2nd max between max and min in strips
      const double emax2  = shower->e2tsts1();
      // E of 1st max in strips
      const double emax   = shower->emaxs1();

      const double sum = emax + emax2;
      
      const double dmax = - (eratio_FF[ieta][ipt] * sum * sum)/
	(eratio_FF[ieta][ipt] * sum - 2 * emax2);

      // (Emax1-Emax2)/(Emax1+Emax2)
      shower->set_emaxs1(dmax);      

    }

    if (m_doEMPID) {
      sc = m_pidBuilder->execute(*ph);
      if (sc.isFailure()) {
	ATH_MSG_ERROR("problem executing EMPIDBuilder on photon Container");
	return sc;
      }
    }
    
  }

  return StatusCode::SUCCESS;
}

int gmsbFudgeFactors::getEtaBin(double eta2) const {
  if (eta2 < 0.6) return 0;
  if (eta2 < 1.37) return 1;
  if (eta2 < 1.52) return 2;
  if (eta2 < 1.81) return 3;
  return 4;
}

int gmsbFudgeFactors::getPtBin(double pt) const {
  if (pt < 20.0*GeV) return 0;
  if (pt < 25.0*GeV) return 1;
  if (pt < 30.0*GeV) return 2;
  if (pt < 35.0*GeV) return 3;
  if (pt < 40.0*GeV) return 4;
  if (pt < 50.0*GeV) return 5;
  if (pt < 60.0*GeV) return 6;
  return 7;
}

// ////////////////////////////////////////////////////////////
//           F U G D E   F A C T O R S   A R R A Y S         //
// ////////////////////////////////////////////////////////////

// ===========================
double gmsbFudgeFactors::rhad_FF[5][8] = {
  
  {-0.00102181,-0.00019482,-0.000864613,-0.000307945,0.000168496,-0.00023272,-0.0015275,-0.000708431},
  
  {-0.000497362,-0.000345779,-0.000125786,-0.000594985,-0.000380814,9.71204e-06,-0.000218846,-6.30434e-05},
  {0,0,0,0,0,0,0,0},
  
  {-0.000559827,-0.00136083,-0.000714184,-0.00185926,-0.00125061,2.4436e-05,-0.000305616,0.000977649},
  
  {-0.00145807,-0.00115328,-0.00312003,-0.004415,-0.00154483,-0.00186636,-0.00503688,-0.00413539}
};

// ===========================
double gmsbFudgeFactors::rhad1_FF[5][8] = {

  {-0.000234666,2.71568e-06,-0.000122919,7.97061e-06,0.000179197,8.74221e-07,-3.40582e-05,-5.48568e-05},
  
  {-0.000237292,-0.000130856,-0.00011744,-0.000178787,-0.000151164,-3.10185e-05,-6.73619e-05,-7.41849e-05},
  {0,0,0,0,0,0,0,0},
  
  {-0.000177096,-1.43329e-05,-0.000147711,-0.000336169,-0.000174246,-0.000392115,0.00050934,8.06176e-06},
  
  {-0.000447519,-5.70643e-05,-0.000577496,-0.000765514,-0.00025275,-0.000284143,-0.000642911,-0.000623379}
};

// ===========================
double gmsbFudgeFactors::reta_FF[5][8] = {
  
  {-0.00179404,-0.00273508,-0.00260949,-0.00218415,-0.00298005,-0.00309342,-0.00331885,-0.00326723},
  
  {-0.00466764,-0.00523251,-0.00478846,-0.0038119,-0.00576085,-0.00629777,-0.005413,-0.00589406},
  {0,0,0,0,0,0,0,0},
  
  {-0.00716525,-0.00791937,-0.00795412,-0.00783455,-0.00841039,-0.00916111,-0.00774229,-0.0105758},
  
  {-0.0077455,-0.00823462,-0.00865459,-0.00788707,-0.00947398,-0.0102837,-0.00900179,-0.0104206}
};

// ===========================
double gmsbFudgeFactors::rphi_FF[5][8] = {

  {0.0021565,-0.00248921,0.00495529,0.00142074,0.00118774,0.000369906,-0.000450075,-0.00172585},
  
  {0.00204033,0.000802815,0.00464308,0.00961536,0.0010739,-0.00257665,-0.00126821,-0.00272733},
  {0,0,0,0,0,0,0,0},
  
  {0.0096001,0.00643647,0.010591,0.0122002,0.00472486,0.00125474,0.00573468,-0.00465041},
  
  {0.00320339,-0.00358915,-0.00235271,-0.00216746,-0.00630444,-0.0070765,-0.00630009,-0.00838786}
};

// ===========================
double gmsbFudgeFactors::weta2_FF[5][8] = {
  
  {0.000193914,0.000215059,0.000214871,0.000199641,0.000229836,0.000260861,0.000260197,0.000228934},
  
  {0.000241264,0.00026471,0.00024403,0.00022511,0.000280658,0.000288717,0.0002569,0.000289068},
  {0,0,0,0,0,0,0,0},
  
  {0.000433117,0.00044233,0.000434412,0.000432006,0.0004546,0.000468174,0.000352059,0.000522199},
  
  {0.000425259,0.00041444,0.000461924,0.000402658,0.000488657,0.000445713,0.000364995,0.000475395}
};

// ===========================
double gmsbFudgeFactors::w1_FF[5][8] = {
  
  {-0.00377476,-0.00405866,-0.00556296,-0.00641072,-0.00636834,-0.00333548,-0.00738984,-0.00804341},
  
  {0.0016641,0.00231802,0.00136536,-0.00257462,0.00140679,0.0014897,-0.00131679,0.00172132},
  {0,0,0,0,0,0,0,0},
  
  {0.00304049,0.00320143,0.00250214,-4.76241e-05,0.00224686,0.00494897,-0.00144821,0.00498015},
  
  {0.00224859,0.00573522,0.00427294,0.00505984,0.00633717,0.00778621,0.00656009,0.00983226}
};

// ===========================
double gmsbFudgeFactors::wtot_FF[5][8] = {
  
  {0.0496336,0.0598663,0.0541141,0.0451692,0.0515084,0.0526413,0.0472215,0.0531849},
  
  {0.0441499,0.0524347,0.0447123,0.0113165,0.0487127,0.0460203,0.0427406,0.0413196},
  {0,0,0,0,0,0,0,0},
  
  {0.0635982,0.0696092,0.0572159,0.0387144,0.0486252,0.088006,0.0289972,0.0750468},
  
  {0.0292448,0.0485572,0.0307806,0.0461965,0.044089,0.0589918,0.0734041,0.0736362}
};

// ===========================
double gmsbFudgeFactors::fracm_FF[5][8] = {

  {0.00040549,0.000841856,-0.000488058,-0.00271371,-0.000277132,0.000767648,-0.00188413,0.001747},
  
  {0.00158006,0.00151023,-0.000322729,-0.00918156,-0.00161466,-0.00348493,-0.00514555,-0.00249958},
  {0,0,0,0,0,0,0,0},
  
  {0.00930536,0.00856665,0.00785327,-0.000518382,-0.000112146,0.00965384,-0.00836787,0.00363204},
  
  {0.00565352,0.0092155,0.00671946,0.00818056,0.00965606,0.0118695,0.0119997,0.0168972}
};

// ===========================
double gmsbFudgeFactors::deltae_FF[5][8] = {

  {-1.31023,-1.01895,-1.52932,-1.02465,-1.28839,-2.02682,-2.1058,-0.305387},
  
  {-1.78933,-1.58073,-2.14103,-1.86159,-2.36409,-2.38288,-4.24965,-3.24585},
  {0,0,0,0,0,0,0,0},
  
  {-2.03735,-2.37012,-2.76465,-2.07479,-1.80408,-2.48758,-1.73586,-0.813705},
  
  {-4.25115,-3.67255,-5.45253,-6.84063,-5.51274,-3.37129,-6.32979,-7.4421}
};

// ===========================
double gmsbFudgeFactors::eratio_FF[5][8] = {

  {-0.00370008,-0.00219446,-0.0015853,-0.00197202,-0.000516117,-0.00219464,-0.00193435,-0.000662088},
  
  {-0.00261867,-0.00127923,-0.00200546,0.000994205,-0.00181395,-0.00111651,0.00070107,-0.000898659},
  {0,0,0,0,0,0,0,0},
  
  {-0.000429332,-0.000913501,0.00137872,0.000607491,0.00138175,-0.00220394,0.000901401,0.00137615},
  
  {-0.00263643,-0.00136226,-0.000427485,0.000298023,-0.00117576,-0.000718951,0.000125647,-0.000399649}
};

// ===========================
