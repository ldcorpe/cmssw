import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackerGeometryWriter")
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.load("Configuration.Geometry.GeometryExtended_cff")

index = process.XMLIdealGeometryESSource.geomXMLFiles.index('Geometry/TrackerCommonData/data/pixfwdMaterials.xml')
process.XMLIdealGeometryESSource.geomXMLFiles.insert(index, 'Geometry/TrackerCommonData/data/trackerParameters.xml')
    
process.source = cms.Source("EmptyIOVSource",
                            lastValue = cms.uint64(1),
                            timetype = cms.string('runnumber'),
                            firstValue = cms.uint64(1),
                            interval = cms.uint64(1)
                            )

process.TrackerGeometricDetExtraESModule = cms.ESProducer( "TrackerGeometricDetExtraESModule",
                                                           fromDDD = cms.bool( True )
                                                           )

process.TrackerGeometryWriter = cms.EDAnalyzer("PGeometricDetBuilder")
process.TrackerGeometryExtraWriter = cms.EDAnalyzer("PGeometricDetExtraBuilder")
process.TrackerParametersWriter = cms.EDAnalyzer("PTrackerParametersDBBuilder")

process.CondDBCommon.BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService')
process.CondDBCommon.timetype = cms.untracked.string('runnumber')
process.CondDBCommon.connect = cms.string('sqlite_file:myfile.db')
process.PoolDBOutputService = cms.Service("PoolDBOutputService",
                                          process.CondDBCommon,
                                          toPut = cms.VPSet(cms.PSet(record = cms.string('IdealGeometryRecord'),tag = cms.string('TKRECO_Geometry_Test02')),
                                                            cms.PSet(record = cms.string('PGeometricDetExtraRcd'),tag = cms.string('TKExtra_Geometry_Test02')),
                                                            cms.PSet(record = cms.string('PTrackerParametersRcd'),tag = cms.string('TK_Parameters_Test02'))
                                                            )
                                          )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
    )

process.p1 = cms.Path(process.TrackerGeometryWriter*process.TrackerGeometryExtraWriter*process.TrackerParametersWriter)

