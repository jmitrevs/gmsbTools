/*****************************************************************************
Name    : gmsbOverlapCheckingTool.cxx
Package : offline/PhysicsAnalysis/AnalysisCommon/gmsbTools
Author  : Ketevi A. Assamagan
Created : November 2007
Purpose : User Analysis Overlap Checking - see gmsbOverlapCheckingTool.h for details
*****************************************************************************/

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/GaudiException.h"
#include "GaudiKernel/Property.h"

// Accessing data:
#include "CLHEP/Units/PhysicalConstants.h"

// User Tools
#include "gmsbTools/gmsbOverlapCheckingTool.h"

#include <sstream>
#include <iomanip>
#include <iostream>


//------------------------------------------------------------------------------
gmsbOverlapCheckingTool::gmsbOverlapCheckingTool( const std::string& type,
						  const std::string& name, 
						  const IInterface* parent )
  : AthAlgTool( type, name, parent ) {
  
  declareInterface<gmsbOverlapCheckingTool>( this );

  declareProperty("OverlapDeltaR",          m_deltaR=0.01);
  declareProperty("OverlapDeltaRWithJets",  m_deltaRWithJets=0.2);

}

//------------------------------------------------------------------------------
StatusCode gmsbOverlapCheckingTool::initialize() {

  ATH_MSG_DEBUG("in initialize()");

  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
StatusCode gmsbOverlapCheckingTool::finalize() {

  ATH_MSG_DEBUG("in finalize()");
  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
gmsbOverlapCheckingTool::~gmsbOverlapCheckingTool()
{}



