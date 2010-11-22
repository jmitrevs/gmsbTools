/*****************************************************************************
Name    : UserAnalysisOverlapCheckingTool.cxx
Package : offline/PhysicsAnalysis/AnalysisCommon/UserAnalysisUtils
Author  : Ketevi A. Assamagan
Created : November 2007
Purpose : User Analysis Overlap Checking - see UserAnalysisOverlapCheckingTool.h for details
*****************************************************************************/

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/GaudiException.h"
#include "GaudiKernel/Property.h"

// Accessing data:
#include "CLHEP/Units/PhysicalConstants.h"

// User Tools
#include "UserAnalysisUtils/UserAnalysisOverlapCheckingTool.h"

#include <sstream>
#include <iomanip>
#include <iostream>

//using namespace Analysis;
//using namespace Rec;
//using namespace std;

//------------------------------------------------------------------------------
UserAnalysisOverlapCheckingTool::UserAnalysisOverlapCheckingTool( const std::string& type,
                                                                  const std::string& name, 
                                                                  const IInterface* parent )
  : AthAlgTool( type, name, parent ) {
  
  declareInterface<UserAnalysisOverlapCheckingTool>( this );

  declareProperty("OverlapDeltaR",          m_deltaR=0.1);
  declareProperty("OverlapDeltaRWithJets",  m_deltaRWithJets=0.2);

}

//------------------------------------------------------------------------------
StatusCode UserAnalysisOverlapCheckingTool::initialize() {

  ATH_MSG_DEBUG("in initialize()");

  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
StatusCode UserAnalysisOverlapCheckingTool::finalize() {

  ATH_MSG_DEBUG("in finalize()");
  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
UserAnalysisOverlapCheckingTool::~UserAnalysisOverlapCheckingTool()
{}



