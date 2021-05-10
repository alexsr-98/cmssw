#include "L1Trigger/L1TMuonOverlapPhase1/plugins/MuonMatcher.h"

MuonMatcher::MuonMatcher(const edm::ParameterSet& iCfg) {
    src_ = iCfg.getParameter<edm::InputTag>("src");
    genParts_token = consumes<std::vector<reco::GenParticle>>(src_);
    produces<edm::ValueMap<float>>("genParticledispEta");
    produces<edm::ValueMap<float>>("genParticledispPhi");
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
  eventSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", propagator);
  
  std::vector<float> valuesEta;
  std::vector<float> valuesPhi;
  edm::Handle<std::vector<reco::GenParticle>> genParticles;  

  iEv.getByToken(genParts_token, genParticles);
  valuesEta.reserve(genParticles->size());
  valuesPhi.reserve(genParticles->size());

  for(std::vector<reco::GenParticle>::const_iterator genPart = genParticles->begin() ; genPart != genParticles->end() ; ++genPart){
    float theEta = 0;
    if (abs(genPart->pdgId()) == 13) theEta = propagateGenPart(genPart);
    float thePhi = 0;
    if (abs(genPart->pdgId()) == 13) thePhi = propagateGenPartPhi(genPart);
    valuesEta.push_back(theEta);
    valuesPhi.push_back(thePhi);
  }
  std::unique_ptr<edm::ValueMap<float>>  outEta(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerEta(*outEta);
  fillerEta.insert(genParticles, valuesEta.begin(), valuesEta.end());
  fillerEta.fill();
  iEv.put(std::move(outEta),"genParticledispEta");

  std::unique_ptr<edm::ValueMap<float>>  outPhi(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerPhi(*outPhi);
  fillerPhi.insert(genParticles, valuesPhi.begin(), valuesPhi.end());
  fillerPhi.fill();
  iEv.put(std::move(outPhi),"genParticledispPhi");
}

float MuonMatcher::propagateGenPart(std::vector<reco::GenParticle>::const_iterator gP){
  int charge = gP->charge();
  GlobalPoint r3GV(gP->vx(), gP->vy(), gP->vz());
  GlobalVector p3GV(gP->px(), gP->py(), gP->pz());
  GlobalTrajectoryParameters tPars(r3GV, p3GV, charge, &*magField);
  FreeTrajectoryState fts = FreeTrajectoryState(tPars); 
  ReferenceCountingPointer<Surface> rpc;
  float preeta = gP->eta();
  if (preeta < -1.24)     rpc = ReferenceCountingPointer<Surface>(new  BoundDisk( GlobalPoint(0.,0.,-790.),  TkRotation<float>(), SimpleDiskBounds( 300., 810., -10., 10. ) ) ); 
  else if (preeta < 1.24) rpc = ReferenceCountingPointer<Surface>(new  BoundCylinder( GlobalPoint(0.,0.,0.), TkRotation<float>(), SimpleCylinderBounds( 500, 500, -900, 900 ) ) );
  else                    rpc = ReferenceCountingPointer<Surface>(new  BoundDisk( GlobalPoint(0.,0.,790.),   TkRotation<float>(), SimpleDiskBounds( 300., 810., -10., 10. ) ) );
  TrajectoryStateOnSurface trackAtRPC = propagator->propagate(fts, *rpc);
  if (!trackAtRPC.isValid()) return -999; //Something broke when propagating
  else return trackAtRPC.globalPosition().eta();
}


float MuonMatcher::propagateGenPartPhi(std::vector<reco::GenParticle>::const_iterator gP){
  int charge = gP->charge();
  GlobalPoint r3GV(gP->vx(), gP->vy(), gP->vz());
  GlobalVector p3GV(gP->px(), gP->py(), gP->pz());
  GlobalTrajectoryParameters tPars(r3GV, p3GV, charge, &*magField);
  FreeTrajectoryState fts = FreeTrajectoryState(tPars);
  ReferenceCountingPointer<Surface> rpc;
  float preeta = gP->eta();
  if (preeta < -1.24)     rpc = ReferenceCountingPointer<Surface>(new  BoundDisk( GlobalPoint(0.,0.,-790.),  TkRotation<float>(), SimpleDiskBounds( 300., 810., -10., 10. ) ) );
  else if (preeta < 1.24) rpc = ReferenceCountingPointer<Surface>(new  BoundCylinder( GlobalPoint(0.,0.,0.), TkRotation<float>(), SimpleCylinderBounds( 500, 500, -900, 900 ) ) );
  else                    rpc = ReferenceCountingPointer<Surface>(new  BoundDisk( GlobalPoint(0.,0.,790.),   TkRotation<float>(), SimpleDiskBounds( 300., 810., -10., 10. ) ) );
  TrajectoryStateOnSurface trackAtRPC = propagator->propagate(fts, *rpc);
  if (!trackAtRPC.isValid()) return -999; //Something broke when propagating
  else return trackAtRPC.globalPosition().phi();
}

DEFINE_FWK_MODULE(MuonMatcher);
