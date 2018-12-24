import FWCore.ParameterSet.Config as cms

tHqLeptonicFilter  = cms.EDFilter("tHqLeptonicFilter",
                           #genParticleTag = cms.InputTag( "prunedGenParticles" )
                           genParticleTag = cms.InputTag( "flashggPrunedGenParticles" )
)
