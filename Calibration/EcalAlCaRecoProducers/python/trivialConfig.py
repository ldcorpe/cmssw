import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
process = cms.Process("test")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'POSTLS170_V5::All'
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )
process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/data/Run2012A/DoubleElectron/AOD/15Apr2014-v2/00000/007A5890-D4DD-E311-88B1-001E673976ED.root"))
process.TFileService = cms.Service("TFileService",fileName = cms.string("flashggTreeWithTags.root"))
process.test = cms.EDAnalyzer('louieTest')
process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('myOutputFile.root'),
 outputCommands = cms.untracked.vstring("drop *" ) )

process.p = cms.Path(process.test)
process.e = cms.EndPath(process.out)

