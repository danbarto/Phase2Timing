import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
process.load("Geometry.CaloEventSetup.CaloTopology_cfi");
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("Configuration.Geometry.GeometryExtended2026D76Reco_cff")
process.load("Validation.MtdValidation.btlDigiHits_cfi")
process.load("Validation.MtdValidation.btlLocalReco_cfi")
process.load("Validation.MtdValidation.etlDigiHits_cfi")
process.load("Validation.MtdValidation.etlLocalReco_cfi")
process.load('Geometry.MTDGeometryBuilder.mtdGeometry_cfi')

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
                                'file:/ceph/cms//store/user/mcitron/ProjectMetis/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-10000mm_privateMC_11X_RECOMINI_v2_generationForPhase2/output_10.root'
    )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple_phase2timing.root")
)

process.demo = cms.EDAnalyzer('Phase2TimingAnalyzer',
                              genParticles    = cms.InputTag("genParticles", "", "HLT"),
                              recoak4PFJets    = cms.InputTag("ak4PFJets", "", "RECO"),
                              ticlTrackstersEM = cms.InputTag("ticlTrackstersEM"),
                              ticlTrackstersMerge = cms.InputTag("ticlTrackstersMerge"),
                              ticlTrackstersHAD = cms.InputTag("ticlTrackstersHAD"),
                              ticlTrackstersTrk = cms.InputTag("ticlTrackstersTrk"),
                              ticlTrackstersTrkEM = cms.InputTag("ticlTrackstersTrkEM"),
                              ebRecHitsColl = cms.InputTag( 'ecalRecHit','EcalRecHitsEB',"RECO" ),
                              mtdBTLRecHitsColl = cms.InputTag( 'mtdRecHits','FTLBarrel',"RECO" ),
                              mtdETLRecHitsColl = cms.InputTag( 'mtdRecHits','FTLEndcap',"RECO" ),
)


process.p = cms.Path(process.demo)
