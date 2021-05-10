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

class MuonMatcher : public edm::EDProducer{
  public:
    MuonMatcher(const edm::ParameterSet&);
    ~MuonMatcher();

  private:
    virtual void beginJob();
    float propagateGenPart(std::vector<reco::GenParticle>::const_iterator);
    float propagateGenPartPhi(std::vector<reco::GenParticle>::const_iterator);
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob();
    edm::InputTag src_;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genParts_token;
    edm::ESHandle<GlobalTrackingGeometry> globalGeometry;
    edm::ESHandle<MagneticField> magField;
    edm::ESHandle<Propagator> propagator;

};
