#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/GeometrySurface/interface/ReferenceCounted.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Matrix/Vector.h"
#include <string>
#include "DataFormats/Common/interface/Association.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "TrackingTools/IPTools/interface/IPTools.h"


class MuonMatcher : public edm::EDProducer{
  public:
    MuonMatcher(const edm::ParameterSet&);
    ~MuonMatcher();

  private:
    virtual void beginJob();
    float propagateGenPart(std::vector<reco::GenParticle>::const_iterator, reco::BeamSpot BeamSpot);
    float propagateGenPartPhi(std::vector<reco::GenParticle>::const_iterator, reco::BeamSpot BeamSpot);
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob();
    edm::InputTag src_;
    const edm::EDGetTokenT<reco::GenParticleCollection> genParts_token;
    edm::ESHandle<GlobalTrackingGeometry> globalGeometry;
    edm::ESHandle<MagneticField> magField;
    edm::ESHandle<Propagator> propagator;
    const edm::EDGetTokenT<reco::GenParticleCollection> Cands_;
    const edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> Asso_;
    const edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> AssoOriginal_;
};
