import FWCore.ParameterSet.Config as cms
dispGenEta  = cms.EDProducer("MuonMatcher",
                             src=cms.InputTag("genParticles","","HLT"),
                            )
