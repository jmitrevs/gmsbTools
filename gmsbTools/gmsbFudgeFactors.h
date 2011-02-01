/********************************************************************

NAME:     gmsgFudgeFactors, based on egammaAODRender

CREATED:  Feb 2011
********************************************************************/


#ifndef GMSBTOOLS_GMSBFUDGEFACTORS_H
#define GMSBTOOLS_GMSBFUDGEFACTORS_H


// INCLUDE HEADER FILES:

#include "AthenaBaseComps/AthAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"

// forward declarations
class PhotonContainer;
class IEMPIDBuilder;
namespace Analysis 
{
    class Electron;
    class Photon;
}

class gmsbFudgeFactors : public AthAlgorithm 
{
public:
	
	/** @brief Default constructor*/
	gmsbFudgeFactors(const std::string& name, ISvcLocator* pSvcLocator);
	
	/** @brief Destructor*/
	~gmsbFudgeFactors();

	/** @brief initialize method*/
	StatusCode initialize();
	/** @brief finalize method*/
	StatusCode finalize();
	/** @brief execute method*/
	StatusCode execute();
	
private:

	int getEtaBin(double eta2) const;
	int getPtBin(double pt) const;
	
	/** @brief photon collection input name*/
	std::string m_PhotonInputName;

	/** @brief whether to do the fudge factors or not */
	bool m_doFudgeFactors;

	/** @brief whether to do the EMPID */
	bool m_doEMPID;

	// others:

	/** @brief */
	bool          m_dump ;

	/** @brief */
	ToolHandle<IEMPIDBuilder> m_pidBuilder;

	// fudge factors
	static double rhad_FF[5][8];
	static double rhad1_FF[5][8];
	static double reta_FF[5][8];
	static double rphi_FF[5][8];
	static double weta2_FF[5][8];
	static double w1_FF[5][8];
	static double wtot_FF[5][8];
	static double fracm_FF[5][8];
	static double deltae_FF[5][8];
	static double eratio_FF[5][8];
};

#endif

