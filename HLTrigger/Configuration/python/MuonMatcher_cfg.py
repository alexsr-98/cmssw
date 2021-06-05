import FWCore.ParameterSet.Config as cms
dispGenEta  = cms.EDProducer("MuonMatcher",
                             inputCollection = cms.InputTag("genParticles","","HLT"),
                             beamSpot = cms.InputTag("hltOnlineBeamSpot"),
                            )
