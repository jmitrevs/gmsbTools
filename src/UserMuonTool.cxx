/*****************************************************************************
Name    : UserMuonTool.cxx
Package : offline/PhysicsAnalysis/AnalysisCommon/UserAnalysisUtils
Author  : Ketevi A. Assamagan
Created : June 2007
Purpose : Muon stage 0 analysis tools - see UserMuonTool.h for details
*****************************************************************************/

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/GaudiException.h"
#include "GaudiKernel/Property.h"

// Accessing data:
#include "CLHEP/Units/PhysicalConstants.h"

// User Tools
#include "UserAnalysisUtils/UserMuonTool.h"

#include <sstream>

//------------------------------------------------------------------------------
UserMuonTool::UserMuonTool( const std::string& type,
                  const std::string& name, 
                  const IInterface* parent )
  : AlgTool( type, name, parent ) {
  declareInterface<UserMuonTool>( this );

  declareProperty("ObjectKeys", m_keys);
}

//------------------------------------------------------------------------------
StatusCode UserMuonTool::initialize() {
  MsgStream mLog(msgSvc(), name());
  mLog << MSG::DEBUG << "intialize() has been called" << endreq;

  /** get a handle of StoreGate for access to the Event Store */
  StatusCode sc = service("StoreGateSvc", m_storeGate);
  if (sc.isFailure()) {
     mLog << MSG::ERROR
          << "Unable to retrieve pointer to StoreGateSvc"
          << endreq;
     return sc;
  }

  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
StatusCode UserMuonTool::finalize() {
  MsgStream mLog(msgSvc(), name());
  mLog << MSG::DEBUG << "finalize() has been called" << endreq;
  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
UserMuonTool::~UserMuonTool()
{}

/** refit one Muon
      - pass the the muon to be refitted
      - get back a new refitted muon with all its associated new tracks, trackParticles, segments and CaloClusters
      - the isolation is also redone automatically
      - the energies in the calorimeter layers are recalculated */
StatusCode UserMuonTool::refit (const Analysis::Muon* oldMuon, Analysis::Muon* newMuon, 
                                std::vector<const Rec::TrackParticle*>& newTrackParticles,
                                std::vector<const Trk::Track*>& newTracks,
                                std::vector<const Trk::Segment>& newSegments,
                                CaloCluster * newCluster)
{

  return StatusCode::SUCCESS;
}

/** refit a collection of muons
      - it will call the refit for each muon in the collection
      - you will pass the container of the muons
      - you will also pass a vector<std::string> it should contains the StoreGates of the new TrackParticleContainer,
        TrackCollection, SegmentCollection and CaloClusterContainer and MuonContainer respectively and in that order. 
        Thus, the size of the vector must be 5
      - The new collections are written to StoreGate: new MuonContainer of refitted Muons, and new Containers of the
        associated TrackParticles, Tracks, Segments and CaloClusters */
StatusCode UserMuonTool::refit (const Analysis::MuonContainer * oldMuonContainer, std::vector<std::string>& newKeys)
{

  return StatusCode::SUCCESS;
}

/** redo calo isolation either on a refitted muon track or with different parameters 
      if you pass a list of deltaR, you get a list of results back */
StatusCode UserMuonTool::redoCaloIsolation(const Analysis::Muon*, const float deltaR, const bool onlyEM, float& value)
{

   return StatusCode::SUCCESS;
}

StatusCode UserMuonTool::redoCaloIsolation(const Analysis::Muon*, const std::vector<float>& deltaR, 
                                           const bool onlyEM, 
                                           std::vector<double>& values)
{

   return StatusCode::SUCCESS;
}

/** redo track isolation ether on a refitted muon track or with different parameters 
      if you pass a list of deltaR you will get back a list of values */
StatusCode UserMuonTool::redoTrackIsolation(const Analysis::Muon*, const float deltaR, float& value)
{

   return StatusCode::SUCCESS;
}

StatusCode UserMuonTool::redoTrackIsolation(const Analysis::Muon*, const std::vector<float>& deltaR, 
                                            std::vector<double>& values)
{

   return StatusCode::SUCCESS;
}

/** redo combinations -> MuonSpectrometer/Inner Detector after refitting for example
      - you pass an std::vector<std::string>: it should contain the StoreGate keys of:
         # the Inner Detector TrackParticles,
         # the key for the Muon TrackParticles extrapolated to the vertex, 
      - on output, you get a new muon collection - the list of StorGate keys of the associated containers:
         # Inner Detector - same as on input
         # Muon Spectrometer to vertex - same as on input
         # CaloClusterContainer keys
         # and the new MuonCollection keys 
      - the new Muon are completely dressed with energy/isolation/calo information */
StatusCode UserMuonTool::recombine(const std::vector<std::string>& oldKeys, std::vector<std::string>& newKeys)
{

   return StatusCode::SUCCESS;
}

/** refit the common vertex to the list of Muons 
       need to think about this in more details */
StatusCode UserMuonTool::refitVertex( const std::vector<const Analysis::Muon*> )
{

   return StatusCode::SUCCESS;
}





