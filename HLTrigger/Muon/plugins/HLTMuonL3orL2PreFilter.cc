/** \class HLTMuonL3orL2PreFilter
 *
 * See header file for documentation
 *
 *
 */

#include "HLTMuonL3orL2PreFilter.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/Candidate/interface/OverlapChecker.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include <iostream>

//
// constructors and destructor
//
using namespace std;
using namespace edm;
using namespace trigger;
using namespace reco;

HLTMuonL3orL2PreFilter::HLTMuonL3orL2PreFilter(const edm::ParameterSet& iConfig)
    : HLTFilter(iConfig),
      candTag_(iConfig.getParameter<edm::InputTag>("CandTag")),
      previousCandTag_(iConfig.getParameter<edm::InputTag>("PreviousCandTag")),
      L2CandTag_(iConfig.getParameter<edm::InputTag>("L2CandTag")),
      beamspotTag_(iConfig.getParameter<edm::InputTag>("BeamSpotTag")),
      min_N_(iConfig.getParameter<int>("MinN")),
      max_Eta_(iConfig.getParameter<double>("MaxEta")),
      min_Nhits_(iConfig.getParameter<int>("MinNhits")),
      max_Dz_(iConfig.getParameter<double>("MaxDz")),
      min_DxySig_(iConfig.getParameter<double>("MinDxySig")),
      min_Pt_(iConfig.getParameter<double>("MinPt")),
      nsigma_Pt_(iConfig.getParameter<double>("NSigmaPt")),
      max_NormalizedChi2_(iConfig.getParameter<double>("MaxNormalizedChi2")),
      max_DXYBeamSpot_(iConfig.getParameter<double>("MaxDXYBeamSpot")),
      min_DXYBeamSpot_(iConfig.getParameter<double>("MinDXYBeamSpot")),
      min_NmuonHits_(iConfig.getParameter<int>("MinNmuonHits")),
      max_PtDifference_(iConfig.getParameter<double>("MaxPtDifference")),
      min_TrackPt_(iConfig.getParameter<double>("MinTrackPt")),
      matchPreviousCand_(iConfig.getParameter<bool>("MatchToPreviousCand")) {
  candToken_ = consumes<reco::RecoChargedCandidateCollection>(candTag_);
  L2CandToken_ = consumes<reco::RecoChargedCandidateCollection>(L2CandTag_);
  previousCandToken_ = consumes<trigger::TriggerFilterObjectWithRefs>(previousCandTag_);
  beamspotToken_ = consumes<reco::BeamSpot>(beamspotTag_);
}

HLTMuonL3orL2PreFilter::~HLTMuonL3orL2PreFilter() = default;

//
// member functions
//
void HLTMuonL3orL2PreFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  makeHLTFilterDescription(desc);
  desc.add<edm::InputTag>("BeamSpotTag", edm::InputTag("hltOfflineBeamSpot"));
  desc.add<edm::InputTag>("CandTag", edm::InputTag("hltL3MuonCandidates"));
  desc.add<edm::InputTag>("L2CandTag", edm::InputTag("hltL2MuonCandidates"));
  desc.add<edm::InputTag>("PreviousCandTag", edm::InputTag(""));
  desc.add<int>("MinN", 1);
  desc.add<double>("MaxEta", 2.5);
  desc.add<int>("MinNhits", 0);
  desc.add<double>("MaxDz", 9999.0);
  desc.add<double>("MinDxySig", -1.0);
  desc.add<double>("MinPt", 3.0);
  desc.add<double>("NSigmaPt", 0.0);
  desc.add<double>("MaxNormalizedChi2", 9999.0);
  desc.add<double>("MaxDXYBeamSpot", 9999.0);
  desc.add<double>("MinDXYBeamSpot", -1.0);
  desc.add<int>("MinNmuonHits", 0);
  desc.add<double>("MaxPtDifference", 9999.0);
  desc.add<double>("MinTrackPt", 0.0);
  desc.add<bool>("MatchToPreviousCand", true);
  descriptions.add("hltMuonL3orL2PreFilter", desc);
}


// ------------ method called to produce the data  ------------
bool HLTMuonL3orL2PreFilter::hltFilter(edm::Event& iEvent,
                                         const edm::EventSetup& iSetup,
                                         trigger::TriggerFilterObjectWithRefs& filterproduct) const {
  // All HLT filters must create and fill an HLT filter object,
  // recording any reconstructed physics objects satisfying (or not)
  // this HLT filter, and place it in the Event.

  if (saveTags())
    filterproduct.addCollectionTag(candTag_);

  // get hold of trks
  Handle<RecoChargedCandidateCollection> mucands;
  iEvent.getByToken(candToken_, mucands);
  Handle<TriggerFilterObjectWithRefs> previousLevelCands;
  iEvent.getByToken(previousCandToken_, previousLevelCands);
  Handle<RecoChargedCandidateCollection> L2mucands;
  iEvent.getByToken(L2CandToken_, L2mucands);
  vector<RecoChargedCandidateRef> vcands;
  if (previousLevelCands.isValid()) {
    previousLevelCands->getObjects(TriggerMuon, vcands);
  }

  Handle<BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(beamspotToken_, recoBeamSpotHandle);
  OverlapChecker overlaps;
  // Number of objects passing the L3 Trigger or L2:
  int n = 0;
  bool firstIteration = true; 
  for (unsigned int iMuL2 = 0; iMuL2 < L2mucands->size(); iMuL2++) {
    RecoChargedCandidateRef L2cand(L2mucands, iMuL2);
    LogDebug("HLTMuonL3orL2PreFilter") << "cand isNonnull " << L2cand.isNonnull();

//    //did this candidate triggered at previous stage.
//    if (matchPreviousCand_ && !triggerdByPreviousLevel(cand, vcands))
//      continue;
    
    if (std::abs(L2cand->eta()) > max_Eta_)
      continue;

    TrackRef L2tk = L2cand->track();
    LogDebug("HLTMuonL3orL2PreFilter") << " Muon in loop, q*pt= " << L2tk->charge() * L2tk->pt() << " ("
                                         << L2cand->charge() * L2cand->pt() << ") "
                                         << ", eta= " << L2tk->eta() << " (" << L2cand->eta() << ") "
                                         << ", hits= " << L2tk->numberOfValidHits() << ", d0= " << L2tk->d0()
                                         << ", dz= " << L2tk->dz();

    // cut on number of hits
    if (L2tk->numberOfValidHits() < min_Nhits_)
      continue;

    //normalizedChi2 cut
    if (L2tk->normalizedChi2() > max_NormalizedChi2_)
      continue;

    if (recoBeamSpotHandle.isValid()) {
      const BeamSpot& beamSpot = *recoBeamSpotHandle;

      //dz cut
      if (std::abs((L2cand->vz() - beamSpot.z0()) -
                   ((L2cand->vx() - beamSpot.x0()) * L2cand->px() + (L2cand->vy() - beamSpot.y0()) * L2cand->py()) /
                       L2cand->pt() * L2cand->pz() / L2cand->pt()) > max_Dz_)
        continue;

      // dxy significance cut (safeguard against bizarre values)
      if (min_DxySig_ > 0 &&
          (L2tk->dxyError() <= 0 || std::abs(L2tk->dxy(beamSpot.position()) / L2tk->dxyError()) < min_DxySig_))
        continue;

      //dxy beamspot cut
      float absDxy = std::abs(L2tk->dxy(beamSpot.position()));
      if (absDxy > max_DXYBeamSpot_ || absDxy < min_DXYBeamSpot_)
        continue;
    }

    //min muon hits cut
    const reco::HitPattern& L2trackHits = L2tk->hitPattern();
    if (L2trackHits.numberOfValidMuonHits() < min_NmuonHits_)
      continue;

    //pt difference cut
    double L2candPt = L2cand->pt();
    double L2trackPt = L2tk->pt();

    if (std::abs(L2candPt - L2trackPt) > max_PtDifference_)
      continue;

    //track pt cut
    if (L2trackPt < min_TrackPt_)
      continue;

    // Pt threshold cut
    double L2pt = L2cand->pt();
    double L2err0 = L2tk->error(0);
    double L2abspar0 = std::abs(L2tk->parameter(0));
    double L2ptLx = L2pt;
    // convert 50% efficiency threshold to 90% efficiency threshold
    if (L2abspar0 > 0)
      L2ptLx += nsigma_Pt_ * L2err0 / L2abspar0 * L2pt;
    LogTrace("HLTMuonL3orL2PreFilter") << " ...Muon in loop, trackkRef pt= " << L2tk->pt() << ", ptLx= " << L2ptLx
                                         << " cand pT " << L2cand->pt();
    if (L2ptLx < min_Pt_)
      continue;

    bool matchedL3 = false;
    
    
    
////////////////////////////////////////L3///////////////////////////////////////////////////////////
    for (unsigned int iMu = 0; iMu < mucands->size(); iMu++) {
    RecoChargedCandidateRef cand(mucands, iMu);
    LogDebug("HLTMuonL3SimplePreFilter") << "cand isNonnull " << cand.isNonnull();

    //did this candidate triggered at previous stage.
    if (matchPreviousCand_ && !triggerdByPreviousLevel(cand, vcands))
      continue;

    if (std::abs(cand->eta()) > max_Eta_)
      continue;

    TrackRef tk = cand->track();
    LogDebug("HLTMuonL3SimplePreFilter") << " Muon in loop, q*pt= " << tk->charge() * tk->pt() << " ("
                                         << cand->charge() * cand->pt() << ") "
                                         << ", eta= " << tk->eta() << " (" << cand->eta() << ") "
                                         << ", hits= " << tk->numberOfValidHits() << ", d0= " << tk->d0()
                                         << ", dz= " << tk->dz();

    // cut on number of hits
    if (tk->numberOfValidHits() < min_Nhits_)
      continue;

    //normalizedChi2 cut
    if (tk->normalizedChi2() > max_NormalizedChi2_)
      continue;

    if (recoBeamSpotHandle.isValid()) {
      const BeamSpot& beamSpot = *recoBeamSpotHandle;

      //dz cut
      if (std::abs((cand->vz() - beamSpot.z0()) -
                   ((cand->vx() - beamSpot.x0()) * cand->px() + (cand->vy() - beamSpot.y0()) * cand->py()) /
                       cand->pt() * cand->pz() / cand->pt()) > max_Dz_)
        continue;

      // dxy significance cut (safeguard against bizarre values)
      if (min_DxySig_ > 0 &&
          (tk->dxyError() <= 0 || std::abs(tk->dxy(beamSpot.position()) / tk->dxyError()) < min_DxySig_))
        continue;

      //dxy beamspot cut
      float absDxy = std::abs(tk->dxy(beamSpot.position()));
      if (absDxy > max_DXYBeamSpot_ || absDxy < min_DXYBeamSpot_)
        continue;
    }

    //min muon hits cut
    const reco::HitPattern& trackHits = tk->hitPattern();
    if (trackHits.numberOfValidMuonHits() < min_NmuonHits_)
      continue;

    //pt difference cut
    double candPt = cand->pt();
    double trackPt = tk->pt();

    if (std::abs(candPt - trackPt) > max_PtDifference_)
      continue;

    //track pt cut
    if (trackPt < min_TrackPt_)
      continue;

    // Pt threshold cut
    double pt = cand->pt();
    double err0 = tk->error(0);
    double abspar0 = std::abs(tk->parameter(0));
    double ptLx = pt;
    // convert 50% efficiency threshold to 90% efficiency threshold
    if (abspar0 > 0)
      ptLx += nsigma_Pt_ * err0 / abspar0 * pt;
    LogTrace("HLTMuonL3SimplePreFilter") << " ...Muon in loop, trackkRef pt= " << tk->pt() << ", ptLx= " << ptLx
                                         << " cand pT " << cand->pt();
    if (ptLx < min_Pt_)
      continue;

    if (firstIteration){
      n++;
      filterproduct.addObject(TriggerMuon, cand);}
    ////CHECK HERE IF L2 AND L3 MUONS ARE THE SAME
    
    
    if (overlaps(*L2cand, *cand)){
    matchedL3 = true;}  
  }  //for iMu

    firstIteration = false;
    if (!matchedL3){
        n++;
        filterproduct.addObject(TriggerMuon, L2cand);
    } 
  }  //for iMuL2
  


  // filter decision
  const bool accept(n >= min_N_);

  LogDebug("HLTMuonL3orL2PreFilter") << " >>>>> Result of HLTMuonL3orL2PreFilter is " << accept
                                       << ", number of muons passing thresholds= " << n;

  return accept;
}

bool HLTMuonL3orL2PreFilter::triggerdByPreviousLevel(const reco::RecoChargedCandidateRef& candref,
                                                       const std::vector<reco::RecoChargedCandidateRef>& vcands) {
  unsigned int i = 0;
  unsigned int i_max = vcands.size();
  for (; i != i_max; ++i) {
    if (candref == vcands[i])
      return true;
  }

  return false;
}

// declare this class as a framework plugin
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HLTMuonL3orL2PreFilter);


