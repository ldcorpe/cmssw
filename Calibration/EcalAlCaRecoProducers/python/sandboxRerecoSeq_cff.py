import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Reconstruction_Data_cff import ecalLocalRecoSequence, pfClusteringPS, pfClusteringECAL, ecalClusters

from RecoLocalCalo.Configuration.RecoLocalCalo_cff import *
recoECAL = cms.Sequence(ecalLocalRecoSequence)


#globalreco = cms.Sequence(particleFlowCluster*ecalClusters*egammaGlobalReco*

#ecalRecHit.EBuncalibRecHitCollection = cms.InputTag("ecalGlobalUncalibRecHit","EcalUncalibRecHitsEB","ALCASKIM")
#ecalRecHit.EEuncalibRecHitCollection = cms.InputTag("ecalGlobalUncalibRecHit","EcalUncalibRecHitsEE","ALCASKIM")
electronRecoSeq = cms.Sequence( ecalLocalRecoSequence)

from RecoParticleFlow.PFClusterProducer.particleFlowCluster_cff import *
#pfClusteringECAL = cms.Sequence(particleFlowRecHitECAL*particleFlowClusterECAL)
rerecoPFClusteringSeq = cms.Sequence(pfClusteringPS + pfClusteringECAL)

from  RecoEcal.Configuration.RecoEcal_cff import *
#correctedHybridSuperClusters.corectedSuperClusterCollection = 'recalibSC'
#correctedMulti5x5SuperClustersWithPreshower.corectedSuperClusterCollection = 'endcapRecalibSC'
#    if(re.match("CMSSW_5_.*",CMSSW_VERSION)):
#        multi5x5PreshowerClusterShape.endcapSClusterProducer = "correctedMulti5x5SuperClustersWithPreshower:endcapRecalibSC"

#    process.load("Calibration.EcalCalibAlgos.electronRecalibSCAssociator_cfi")
#from Calibration.EcalAlCaRecoProducers.electronRecalibSCAssociatorSH_cfi import *
from Calibration.EcalCalibAlgos.electronRecalibSCAssociator_cfi import *

#electronRecalibSCAssociator.scIslandCollection = cms.string('endcapRecalibSC')
#electronRecalibSCAssociator.scIslandProducer = cms.string('correctedMulti5x5SuperClustersWithPreshower')
#electronRecalibSCAssociator.scProducer = cms.string('correctedHybridSuperClusters')
#electronRecalibSCAssociator.scCollection = cms.string('recalibSC')
#electronRecalibSCAssociator.electronProducer = 'gedGsfElectrons'
electronClusteringSeq = cms.Sequence(ecalClusters * electronRecalibSCAssociator)


sandboxRerecoSeq = cms.Sequence(electronRecoSeq * electronClusteringSeq)
sandboxPFRerecoSeq = cms.Sequence(electronRecoSeq * rerecoPFClusteringSeq * electronClusteringSeq)
