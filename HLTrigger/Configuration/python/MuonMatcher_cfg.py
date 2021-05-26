import FWCore.ParameterSet.Config as cms
dispGenEta  = cms.EDProducer("MuonMatcher",
                             src=cms.InputTag("genParticles","","HLT"),
                             inputCollection = cms.InputTag("genParticles","","HLT"),
                             inputOriginal = cms.InputTag("genParticles","","HLT"),
                             beamSpot = cms.InputTag("hltOnlineBeamSpot"),
                            )
