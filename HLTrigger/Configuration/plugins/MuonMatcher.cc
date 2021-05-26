#include "HLTrigger/Configuration/plugins/MuonMatcher.h"


#include "iostream"

using namespace std;

MuonMatcher::MuonMatcher(const edm::ParameterSet& iCfg) : 
    genParts_token(consumes<std::vector<reco::GenParticle>>(iCfg.getParameter<edm::InputTag>("inputOriginal"))),
    Cands_(consumes<reco::GenParticleCollection>(iCfg.getParameter<edm::InputTag>("inputCollection"))),
    BeamSpot_(consumes<reco::BeamSpot>(iCfg.getParameter<edm::InputTag>("beamSpot")))
{
    src_ = iCfg.getParameter<edm::InputTag>("src");
    produces<std::vector<pat::PackedGenParticle>>();
    produces<edm::Association<std::vector<pat::PackedGenParticle>>>();
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
  
  edm::Handle<std::vector<reco::GenParticle>> genParticles; 
  iEv.getByToken(genParts_token, genParticles); //Get the genParticles from the event info 
  
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
  std::vector<int> mapping(genParticles->size(), -1);
  auto outPtrP = std::make_unique<std::vector<pat::PackedGenParticle>>(); //New collection to be filled, in this case with muons 
  for (unsigned int ic = 0, nc = cands->size(); ic < nc; ++ic) {
    reco::GenParticle cand = (*cands)[ic];
    if (abs(cand.pdgId()) == 13 && (cand.status() == 1 || cand.status() == 23)) { //Selection of muons with status 23 or 1
    math::PtEtaPhiMLorentzVector FourMomentum(cand.pt(),propagateGenPart(cand, BeamSpot),propagateGenPartPhi(cand, BeamSpot),cand.mass()); //TLorentzVector with the eta and phi propagated
    cand.setP4(FourMomentum); //The new 4-momenta of the muon
    outPtrP->push_back(pat::PackedGenParticle(cand, edm::Ref<std::vector<reco::GenParticle>>()));
   }
   }
   edm::OrphanHandle<std::vector<pat::PackedGenParticle>> oh = iEv.put(std::move(outPtrP)); //Put the new collection in the event
   auto gp2pgp = std::make_unique<edm::Association<std::vector<pat::PackedGenParticle>>>(oh);
   edm::Association<std::vector<pat::PackedGenParticle>>::Filler gp2pgpFiller(*gp2pgp);
   gp2pgpFiller.insert(genParticles, mapping.begin(), mapping.end());
   gp2pgpFiller.fill();
   iEv.put(std::move(gp2pgp));
}

float MuonMatcher::propagateGenPart(reco::GenParticle gP, reco::BeamSpot BeamSpot){
  int charge = gP.charge();
  GlobalPoint r3GV(gP.vx(), gP.vy(), gP.vz());
  GlobalVector p3GV(gP.px(), gP.py(), gP.pz());
  GlobalTrajectoryParameters tPars(r3GV, p3GV, charge, &*magField);
  FreeTrajectoryState fts = FreeTrajectoryState(tPars); 
  FreeTrajectoryState trackAtBeamSpot = propagator->propagate(fts, BeamSpot);
  if (trackAtBeamSpot.position().eta() == 0) return -10; //Something broke when propagating
  else return trackAtBeamSpot.position().eta();

}


float MuonMatcher::propagateGenPartPhi(reco::GenParticle gP, reco::BeamSpot BeamSpot){
  int charge = gP.charge();
  GlobalPoint r3GV(gP.vx(), gP.vy(), gP.vz());
  GlobalVector p3GV(gP.px(), gP.py(), gP.pz());
  GlobalTrajectoryParameters tPars(r3GV, p3GV, charge, &*magField);
  FreeTrajectoryState fts = FreeTrajectoryState(tPars);
  FreeTrajectoryState trackAtBeamSpot = propagator->propagate(fts, BeamSpot);
  if (trackAtBeamSpot.position().eta() == 0) return -10; //Something broke when propagating
  else return trackAtBeamSpot.position().phi();  

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonMatcher);
