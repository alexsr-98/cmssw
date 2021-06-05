#include "HLTrigger/Configuration/plugins/MuonMatcher.h"


#include "iostream"

using namespace std;

MuonMatcher::MuonMatcher(const edm::ParameterSet& iCfg) : 
    Cands_(consumes<std::vector<reco::GenParticle>>(iCfg.getParameter<edm::InputTag>("inputCollection"))),
    BeamSpot_(consumes<reco::BeamSpot>(iCfg.getParameter<edm::InputTag>("beamSpot")))
{
    //produces<std::vector<pat::PackedGenParticle>>();
    produces<std::vector<reco::GenParticle>>();
}

MuonMatcher::~MuonMatcher(){
}

void MuonMatcher::beginJob(){
}

void MuonMatcher::endJob(){
}

void MuonMatcher::produce(edm::Event& iEv, const edm::EventSetup& eventSetup){
  eventSetup.get<GlobalTrackingGeometryRecord>().get(globalGeometry);
  eventSetup.get<IdealMagneticFieldRecord>().get(magField);
  eventSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorOpposite", propagator); //Opposite propagator to extrapolate to the beam spot
  
  
  edm::Handle<std::vector<reco::GenParticle>> cands;
  iEv.getByToken(Cands_, cands); //Candidate collection
  
  
  //--BeamSpot--
  reco::BeamSpot BeamSpot;
  edm::Handle<reco::BeamSpot> beamSpot;
  iEv.getByToken(BeamSpot_, beamSpot);

  if ( beamSpot.isValid() )
  {
      BeamSpot = *beamSpot;
  
  } else
  {
      edm::LogInfo("MyAnalyzer")
        << "No beam spot available from EventSetup \n";
  }

  auto outPtrP = std::make_unique<std::vector<reco::GenParticle>>(); //New collection to be filled, in this case with muons 
  for (unsigned int ic = 0, nc = cands->size(); ic < nc; ++ic) {
    reco::GenParticle cand = (*cands)[ic];
    if (abs(cand.pdgId()) == 13 && (cand.status() == 1 || cand.status() == 23)) { //Selection of muons with status 23 or 1
    outPtrP->push_back(reco::GenParticle(cand.charge(), propagateGenPartMomentum(cand, BeamSpot), propagateGenPartVertex(cand, BeamSpot), cand.pdgId(), cand.status(), true));
   }
   }
   iEv.put(std::move(outPtrP)); //Put the new collection in the event
}

math::XYZTLorentzVectorD MuonMatcher::propagateGenPartMomentum(reco::GenParticle gP, reco::BeamSpot BeamSpot){
  int charge = gP.charge();
  GlobalPoint r3GV(gP.vx(), gP.vy(), gP.vz());
  GlobalVector p3GV(gP.px(), gP.py(), gP.pz());
  GlobalTrajectoryParameters tPars(r3GV, p3GV, charge, &*magField);
  FreeTrajectoryState fts = FreeTrajectoryState(tPars); 
  std::pair<FreeTrajectoryState, double> pairPropag = propagator->propagateWithPath(fts, BeamSpot);
  FreeTrajectoryState trackAtBeamSpot = pairPropag.first;
//  if (trackAtBeamSpot.position().eta() == 0.0) return -10; //Something broke when propagating
//  else return trackAtBeamSpot.position().eta();
  double px = trackAtBeamSpot.momentum().x();
  double py = trackAtBeamSpot.momentum().y();
  double pz = trackAtBeamSpot.momentum().z();
  double E = sqrt(px*px+py*py+pz*pz + gP.mass()*gP.mass());
  math::XYZTLorentzVectorD FourMomentum(px,py,pz,E);
  return FourMomentum;
}

math::XYZPoint MuonMatcher::propagateGenPartVertex(reco::GenParticle gP, reco::BeamSpot BeamSpot){
  int charge = gP.charge();
  GlobalPoint r3GV(gP.vx(), gP.vy(), gP.vz());
  GlobalVector p3GV(gP.px(), gP.py(), gP.pz());
  GlobalTrajectoryParameters tPars(r3GV, p3GV, charge, &*magField);
  FreeTrajectoryState fts = FreeTrajectoryState(tPars); 
  std::pair<FreeTrajectoryState, double> pairPropag = propagator->propagateWithPath(fts, BeamSpot);
  FreeTrajectoryState trackAtBeamSpot = pairPropag.first;
//  if (trackAtBeamSpot.position().eta() == 0.0) return -10; //Something broke when propagating
//  else return trackAtBeamSpot.position().eta();
  double x = trackAtBeamSpot.position().x();
  double y = trackAtBeamSpot.position().y();
  double z = trackAtBeamSpot.position().z();
  math::XYZPoint vertex(x, y, z);
  return vertex;
}





DEFINE_FWK_MODULE(MuonMatcher);
