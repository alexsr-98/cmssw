  #include "HLTrigger/Configuration/plugins/MuonMatcher.h"


#include "iostream"

using namespace std;

MuonMatcher::MuonMatcher(const edm::ParameterSet& iCfg) : 
    genParts_token(consumes<reco::GenParticleCollection>(iCfg.getParameter<edm::InputTag>("inputOriginal"))),
    Cands_(consumes<reco::GenParticleCollection>(iCfg.getParameter<edm::InputTag>("inputCollection"))),
    Asso_(consumes<edm::Association<reco::GenParticleCollection>>(iCfg.getParameter<edm::InputTag>("Map"))),
    AssoOriginal_(consumes<edm::Association<reco::GenParticleCollection>>(iCfg.getParameter<edm::InputTag>("inputCollection")))
{
    src_ = iCfg.getParameter<edm::InputTag>("src");
    produces<edm::ValueMap<float>>("genParticledispEta");
    produces<edm::ValueMap<float>>("genParticledispPhi");
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
  eventSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorOpposite", propagator);
  
  std::vector<float> valuesEta;
  std::vector<float> valuesPhi;
  edm::Handle<reco::GenParticleCollection> genParticles; 
  reco::BeamSpot BeamSpot; //check this

  iEv.getByToken(genParts_token, genParticles);
  valuesEta.reserve(genParticles->size());
  valuesPhi.reserve(genParticles->size());
  
  edm::Handle<reco::GenParticleCollection> cands;
  iEv.getByToken(Cands_, cands);
  
  edm::Handle<edm::Association<reco::GenParticleCollection>> asso;
  iEv.getByToken(Asso_, asso);
  
  edm::Handle<edm::Association<reco::GenParticleCollection>> assoOriginal;
  iEv.getByToken(AssoOriginal_, assoOriginal);
  
  std::vector<int> mapping(genParticles->size(), -1);
  //invert the value map from Orig2New to New2Orig
  std::map<edm::Ref<reco::GenParticleCollection>, edm::Ref<reco::GenParticleCollection>> reverseMap;
  
for (unsigned int ic = 0, nc = genParticles->size(); ic < nc; ++ic) {
cout << nc  << endl;
  edm::Ref<reco::GenParticleCollection> originalRef = edm::Ref<reco::GenParticleCollection>(genParticles, ic);
  edm::Ref<reco::GenParticleCollection> newRef = (*assoOriginal)[originalRef];
  reverseMap.insert(
      std::pair<edm::Ref<reco::GenParticleCollection>, edm::Ref<reco::GenParticleCollection>>(newRef, originalRef));
   }
  auto outPtrP = std::make_unique<std::vector<pat::PackedGenParticle>>();
  unsigned int packed = 0;
  for (unsigned int ic = 0, nc = cands->size(); ic < nc; ++ic) {
    const reco::GenParticle& cand = (*cands)[ic];
    if (true) {
    cout << "QUE PASOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO" << "\n";
    // Obtain original gen particle collection reference from input reference and map
    edm::Ref<reco::GenParticleCollection> inputRef = edm::Ref<reco::GenParticleCollection>(cands, ic);
    edm::Ref<reco::GenParticleCollection> originalRef = reverseMap[inputRef];
    edm::Ref<reco::GenParticleCollection> finalPrunedRef = (*asso)[inputRef];
    mapping[originalRef.key()] = packed;
    packed++;
    if (finalPrunedRef.isNonnull()) {  //this particle exists also in the final pruned
    outPtrP->push_back(pat::PackedGenParticle(cand, finalPrunedRef));
    } else {
    if (cand.numberOfMothers() > 0) {
    edm::Ref<reco::GenParticleCollection> newRef = (*asso)[cand.motherRef(0)];
    outPtrP->push_back(pat::PackedGenParticle(cand, newRef));
    } else {
    outPtrP->push_back(pat::PackedGenParticle(cand, edm::Ref<reco::GenParticleCollection>()));
    }
   }
   }
   }
   edm::OrphanHandle<std::vector<pat::PackedGenParticle>> oh = iEv.put(std::move(outPtrP));
   auto gp2pgp = std::make_unique<edm::Association<std::vector<pat::PackedGenParticle>>>(oh);
   edm::Association<std::vector<pat::PackedGenParticle>>::Filler gp2pgpFiller(*gp2pgp);
   gp2pgpFiller.insert(genParticles, mapping.begin(), mapping.end());
   gp2pgpFiller.fill();
   iEv.put(std::move(gp2pgp));
   
  
   
   
   
   
   
  for(reco::GenParticleCollection::const_iterator genPart = genParticles->begin() ; genPart != genParticles->end() ; ++genPart){
    float theEta = -5;
    if (abs(genPart->pdgId()) == 13) theEta = propagateGenPart(genPart, BeamSpot);
    float thePhi = -5;
    if (abs(genPart->pdgId()) == 13) thePhi = propagateGenPartPhi(genPart, BeamSpot);
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

float MuonMatcher::propagateGenPart(reco::GenParticleCollection::const_iterator gP, reco::BeamSpot BeamSpot){
  int charge = gP->charge();
  GlobalPoint r3GV(gP->vx(), gP->vy(), gP->vz());
  GlobalVector p3GV(gP->px(), gP->py(), gP->pz());
  GlobalTrajectoryParameters tPars(r3GV, p3GV, charge, &*magField);
  FreeTrajectoryState fts = FreeTrajectoryState(tPars); 
  FreeTrajectoryState trackAtBeamSpot = propagator->propagate(fts, BeamSpot);
  //if (!trackAtBeamSpot.isValid()) return -999; //Something broke when propagating
  //else return trackAtBeamSpot.globalPosition().eta();
  return trackAtBeamSpot.position().eta();
}


float MuonMatcher::propagateGenPartPhi(reco::GenParticleCollection::const_iterator gP, reco::BeamSpot BeamSpot){
  int charge = gP->charge();
  GlobalPoint r3GV(gP->vx(), gP->vy(), gP->vz());
  GlobalVector p3GV(gP->px(), gP->py(), gP->pz());
  GlobalTrajectoryParameters tPars(r3GV, p3GV, charge, &*magField);
  FreeTrajectoryState fts = FreeTrajectoryState(tPars);
  FreeTrajectoryState trackAtBeamSpot = propagator->propagate(fts, BeamSpot);
//  if (!trackAtRPC.isValid()) return -999; //Something broke when propagating
//  else return trackAtRPC.globalPosition().phi();
  return trackAtBeamSpot.position().phi();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonMatcher);
