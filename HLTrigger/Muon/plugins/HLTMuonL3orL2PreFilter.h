#ifndef HLTMuonL3orL2PreFilter_h
#define HLTMuonL3orL2PreFilter_h

/** \class HLTMuonL3orL2PreFilter
 *
 *
 *  This class is an HLTFilter (-> EDFilter) implementing
 *  a simple filter to select L3 muons and L2 if they dont match with any L3
 * 
 *  Original author:  A. Soto <alejandro.soto.rodriguez@cern.ch>
 *
 */

#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/Track.h"

class HLTMuonL3orL2PreFilter : public HLTFilter {
public:
  explicit HLTMuonL3orL2PreFilter(const edm::ParameterSet &);
  ~HLTMuonL3orL2PreFilter() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  bool hltFilter(edm::Event &,
                 const edm::EventSetup &,
                 trigger::TriggerFilterObjectWithRefs &filterproduct) const override;
  
  bool checkOverlap(const reco::RecoChargedCandidate r1, const reco::RecoChargedCandidate r2) const {
    return (deltaR(r1,r2)<0.05);}
private:
  static bool triggerdByPreviousLevel(const reco::RecoChargedCandidateRef &,
                                      const std::vector<reco::RecoChargedCandidateRef> &);

  edm::InputTag candTag_;                                             // input tag identifying muon container
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> candToken_;  // token identifying muon container
  edm::InputTag previousCandTag_;  // input tag identifying product contains muons passing the previous level
  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs>
      previousCandToken_;  // token identifying product contains muons passing the previous level
  edm::InputTag L2CandTag_;                                             // input tag identifying muon container
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> L2CandToken_;  // token identifying muon container
  edm::InputTag beamspotTag_;
  edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;

  const int min_N_;                  // minimum number of muons to fire the trigger
  const double max_Eta_;             // Eta cut
  const int min_Nhits_;              // threshold on number of hits on muon
  const double max_Dz_;              // dz cut
  const double min_DxySig_;          // dxy significance cut
  const double min_Pt_;              // pt threshold in GeV
  const double nsigma_Pt_;           // pt uncertainty margin (in number of sigmas)
  const double max_NormalizedChi2_;  // cutoff in normalized chi2
  const double max_DXYBeamSpot_;     // cutoff in dxy from the beamspot
  const double min_DXYBeamSpot_;     // minimum cut on dxy from the beamspot
  const int min_NmuonHits_;          // cutoff in minumum number of chi2 hits
  const double max_PtDifference_;    // cutoff in maximum different between global track and tracker track
  const double min_TrackPt_;         // cutoff in tracker track pt
  bool matchPreviousCand_;
};
#endif  //HLTMuonL3orL2PreFilter_h