# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: step3 -s RAW2DIGI,RECO,RECOSIM,PAT,VALIDATION:@phase2Validation+@miniAODValidation,DQM:@phase2+@miniAODDQM --conditions auto:phase2_realistic_T21 --datatier GEN-SIM-RECO,MINIAODSIM,DQMIO -n 10 --eventcontent FEVTDEBUGHLT,MINIAODSIM,DQM --geometry Extended2026D95 --era Phase2C17I13M9 --no_exec --filein file:step2.root --fileout file:step3.root
import numpy as np
from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
from Configuration.AlCa.GlobalTag import GlobalTag
import FWCore.ParameterSet.Config as cms
# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
# Reconstruction
from RecoHGCal.TICL.iterativeTICL_cff import *
from RecoLocalCalo.HGCalRecProducers.hgcalLayerClusters_cff import hgcalLayerClustersEE, hgcalLayerClustersHSi, hgcalLayerClustersHSci
from RecoLocalCalo.HGCalRecProducers.hgcalMergeLayerClusters_cfi import hgcalMergeLayerClusters
from RecoHGCal.TICL.ticlDumperSimple_cfi import ticlDumperSimple
# Validation
#from Validation.HGCalValidation.HGCalValidator_cfi import *
from RecoLocalCalo.HGCalRecProducers.recHitMapProducer_cfi import recHitMapProducer 

# Load DNN ESSource
from RecoTracker.IterativeTracking.iterativeTk_cff import trackdnn_source

# Automatic addition of the customisation function from RecoHGCal.Configuration.RecoHGCal_EventContent_cff
from RecoHGCal.Configuration.RecoHGCal_EventContent_cff import customiseHGCalOnlyEventContent
from SimCalorimetry.HGCalAssociatorProducers.simTracksterAssociatorByEnergyScore_cfi import simTracksterAssociatorByEnergyScore as simTsAssocByEnergyScoreProducer
# from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import
# from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociationByLCs_cfi import tracksterSimTracksterFromCPsAssociationLinking, tracksterSimTracksterAssociationLinking, tracksterSimTracksterFromCPsAssociationPR, tracksterSimTracksterAssociationPR
from RecoHGCal.TICL.TICLSeedingRegions_cff import ticlSeedingGlobal
#from RecoHGCal.TICL.simpleValidation_cfi import simpleValidation

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5
from FWCore.ParameterSet.VarParsing import VarParsing
from utils import read_csv
import ROOT
ROOT.EnableImplicitMT(0)  # Ensures only CMSSW manages threads


# VarParsing instance
options = VarParsing('analysis')

# Custom options
options.register('parametersFile',
                 'default/default_params.csv',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 'Name of parameters file')

options.register('nEvents',
                 100,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 'Number of events')

# options.register('outputFile',
#              'temp/simple_validation.root',
#              VarParsing.multiplicity.singleton,
#              VarParsing.varType.string,
#              'output file validation')

# options.register('inputFile',
#               'file:input/step2.root',
#               VarParsing.multiplicity.singleton,
#               VarParsing.varType.string,
#               'Name of input file')

options.parseArguments()


process = cms.Process('RECO3', Phase2C17I13M9, ticl_v5)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMServices.Core.DQMStoreNonLegacy_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
print(f"Running over {options.nEvents}")
process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(options.nEvents),
    output=cms.optional.untracked.allowed(cms.int32, cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(
                                [
                               'file:/data/wredjeb/TestCP/CMSSW_15_0_0_pre1/src/29896.203_CloseByPGun_CE_E_Front_120um+Run4D110PU_ticl_v5/step3.root'
                            #'file:/data/wredjeb/Optimizer/CMSSW_15_0_0_pre2/src/29696.203_CloseByPGun_CE_E_Front_120um+Run4D110_ticl_v5/step3.root'
                                ]
                            ),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
                            bypassVersionCheck = cms.untracked.bool(True)
                            )

process.options = cms.untracked.PSet(
    IgnoreCompletely=cms.untracked.vstring(),
    Rethrow=cms.untracked.vstring(),
    accelerators=cms.untracked.vstring('*'),
    allowUnscheduled=cms.obsolete.untracked.bool,
    canDeleteEarly=cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules=cms.untracked.bool(True),
    dumpOptions=cms.untracked.bool(False),
    emptyRunLumiMode=cms.obsolete.untracked.string,
    eventSetup=cms.untracked.PSet(
        forceNumberOfConcurrentIOVs=cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs=cms.untracked.uint32(0)
    ),
    fileMode=cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun=cms.untracked.bool(False),
    holdsReferencesToDeleteEarly=cms.untracked.VPSet(),
    makeTriggerResults=cms.obsolete.untracked.bool,
    modulesToIgnoreForDeleteEarly=cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(0),
    numberOfConcurrentRuns=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(0),
    numberOfThreads=cms.untracked.uint32(1),
    printDependencies=cms.untracked.bool(False),
    sizeOfStackForThreadsInKB=cms.optional.untracked.uint32,
    throwIfIllegalParameter=cms.untracked.bool(True),
    wantSummary=cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation=cms.untracked.string('step3 nevts:10'),
    name=cms.untracked.string('Applications'),
    version=cms.untracked.string('$Revision: 1.19 $')
)


# Additional output definition

# Other statements
process.mix.playback = True
process.mix.digitizers = cms.PSet()
for a in process.aliases:
    delattr(process, a)
process.RandomNumberGeneratorService.restoreStateLabel = cms.untracked.string(
    "randomEngineStateProducer")
process.GlobalTag = GlobalTag(
    process.GlobalTag, 'auto:phase2_realistic_T33', '')

# Path and EndPath definitions
# process.filteredLayerClustersCLUE3DEM = cms.EDProducer("FilteredLayerClustersProducer",
#  LayerClusters = cms.InputTag("hgcalMergeLayerClusters"),
#  LayerClustersInputMask = cms.InputTag("hgcalMergeLayerClusters","InitialLayerClustersMask"),
#  algo_number = cms.vint32(6, 7),
#  clusterFilter = cms.string('ClusterFilterByAlgoAndSizeAndLayerRange'),
#  iteration_label = cms.string('CLUE3DEM'),
#  max_cluster_size = cms.int32(9999),
#  max_layerId = cms.int32(99999),
#  mightGet = cms.optional.untracked.vstring,
#  min_cluster_size = cms.int32(2),
#  min_layerId = cms.int32(0)
# )
process.ticlSeedingRegionProducer = cms.Task(ticlSeedingGlobal)
# Path and EndPath definitions
process.filteredLayerClustersCLUE3DEM = cms.EDProducer("FilteredLayerClustersProducer",
                                                       LayerClusters=cms.InputTag(
                                                           "hgcalMergeLayerClusters"),
                                                       LayerClustersInputMask=cms.InputTag(
                                                           "hgcalMergeLayerClusters", "InitialLayerClustersMask"),
                                                       algo_number=cms.vint32(
                                                           6, 7, 8),
                                                       clusterFilter=cms.string(
                                                           'ClusterFilterByAlgoAndSize'),
                                                       iteration_label=cms.string(
                                                           'CLUE3DEM'),
                                                       max_cluster_size=cms.int32(
                                                           9999),
                                                       max_layerId=cms.int32(
                                                           9999),
                                                       mightGet=cms.optional.untracked.vstring,
                                                       min_cluster_size=cms.int32(
                                                           2),
                                                       min_layerId=cms.int32(0)
                                                       )


v = [9.00000000e-01, 9.00000000e-01, 6.13107466e-01, 1.00000000e-01,
     3.55379138e-01, 1.00000000e-01, 1.00000000e+00, 2.50000000e+00,
     2.50000000e+00, 7.00000000e+00, 7.00000000e+00, 4.66002373e+00,
     5.00000000e+00, 2.00000000e+00, 5.00000000e+00, 4.16576578e+00,
     5.00000000e+00, 5.00000000e+00, 1.00000000e-01, 3.07354051e-01,
     4.00000000e-01, 4.00000000e+00, 4.00000000e+00, 3.88817765e+00,
     1.00000000e+00, 1.00000000e+00, 1.37603455e+00, 5.27500000e-01,
     9.02173913e-01, 5.96401028e-01, 5.42743539e-01, 1.69369816e-01,
     2.44698206e-03, 7.26500000e-01]


params = read_csv(options.parametersFile)
totalTask = len(params)
for i, p in enumerate(params):
    setattr(process, 'ticlTrackstersCLUE3D', cms.EDProducer('TrackstersProducer',
                                                                     detector=cms.string(
                                                                         'HGCAL'),
                                                                     filtered_mask=cms.InputTag(
                                                                         "filteredLayerClustersCLUE3DHigh", "CLUE3DHigh"),
                                                                     inferenceAlgo=cms.string(
                                                                         'TracksterInferenceByCNNv4'),
                                                                     itername=cms.string(
                                                                         'CLUE3DHigh'),
                                                                     layer_clusters=cms.InputTag(
                                                                         "hgcalMergeLayerClusters"),
                                                                     layer_clusters_hfnose_tiles=cms.InputTag(
                                                                         "ticlLayerTileHFNose"),
                                                                     layer_clusters_tiles=cms.InputTag(
                                                                         "ticlLayerTileProducer"),
                                                                     mightGet=cms.optional.untracked.vstring,
                                                                     original_mask=cms.InputTag(
                                                                         "hgcalMergeLayerClusters", "InitialLayerClustersMask"),
                                                                     patternRecognitionBy=cms.string(
                                                                         'CLUE3D'),
                                                                     pluginInferenceAlgoTracksterInferenceByANN=cms.PSet(
                                                                         algo_verbosity=cms.int32(
                                                                             0),
                                                                         type=cms.string(
                                                                             'TracksterInferenceByANN')
                                                                     ),
                                                                     pluginInferenceAlgoTracksterInferenceByCNNv4=cms.PSet(
                                                                         algo_verbosity=cms.int32(
                                                                             0),
                                                                         doPID=cms.int32(
                                                                             1),
                                                                         doRegression=cms.int32(
                                                                             0),
                                                                         eid_min_cluster_energy=cms.double(
                                                                             1),
                                                                         eid_n_clusters=cms.int32(
                                                                             10),
                                                                         eid_n_layers=cms.int32(
                                                                             50),
                                                                         inputNames=cms.vstring(
                                                                             'input:0'),
                                                                         onnxModelPath=cms.FileInPath(
                                                                             'RecoHGCal/TICL/data/ticlv4/onnx_models/energy_id_v0.onnx'),
                                                                         outputNames=cms.vstring(
                                                                             'output/regressed_energy:0',
                                                                             'output/id_probabilities:0'
                                                                         ),
                                                                         type=cms.string(
                                                                             'TracksterInferenceByCNNv4')
                                                                     ),
                                                                     pluginInferenceAlgoTracksterInferenceByDNN=cms.PSet(
                                                                         algo_verbosity=cms.int32(
                                                                             0),
                                                                         doPID=cms.int32(
                                                                             1),
                                                                         doRegression=cms.int32(
                                                                             0),
                                                                         eid_min_cluster_energy=cms.double(
                                                                             1),
                                                                         eid_n_clusters=cms.int32(
                                                                             10),
                                                                         eid_n_layers=cms.int32(
                                                                             50),
                                                                         inputNames=cms.vstring(
                                                                             'input'),
                                                                         onnxEnergyModelPath=cms.FileInPath(
                                                                             'RecoHGCal/TICL/data/ticlv5/onnx_models/patternrecognition/energy_v0.onnx'),
                                                                         onnxPIDModelPath=cms.FileInPath(
                                                                             'RecoHGCal/TICL/data/ticlv5/onnx_models/patternrecognition/id_v0.onnx'),
                                                                         output_en=cms.vstring(
                                                                             'enreg_output'),
                                                                         output_id=cms.vstring(
                                                                             'pid_output'),
                                                                         type=cms.string(
                                                                             'TracksterInferenceByDNN')
                                                                     ),
                                                                     pluginPatternRecognitionByCA=cms.PSet(
                                                                         algo_verbosity=cms.int32(
                                                                             0),
                                                                         computeLocalTime=cms.bool(
                                                                             False),
                                                                         energy_em_over_total_threshold=cms.double(
                                                                             -1),
                                                                         etaLimitIncreaseWindow=cms.double(
                                                                             2.1),
                                                                         filter_on_categories=cms.vint32(
                                                                             0),
                                                                         max_delta_time=cms.double(
                                                                             3),
                                                                         max_longitudinal_sigmaPCA=cms.double(
                                                                             9999),
                                                                         max_missing_layers_in_trackster=cms.int32(
                                                                             9999),
                                                                         max_out_in_hops=cms.int32(
                                                                             10),
                                                                         min_cos_pointing=cms.double(
                                                                             -1),
                                                                         min_cos_theta=cms.double(
                                                                             0.915),
                                                                         min_layers_per_trackster=cms.int32(
                                                                             10),
                                                                         oneTracksterPerTrackSeed=cms.bool(
                                                                             False),
                                                                         out_in_dfs=cms.bool(
                                                                             True),
                                                                         pid_threshold=cms.double(
                                                                             0),
                                                                         promoteEmptyRegionToTrackster=cms.bool(
                                                                             False),
                                                                         root_doublet_max_distance_from_seed_squared=cms.double(
                                                                             9999),
                                                                         shower_start_max_layer=cms.int32(
                                                                             9999),
                                                                         siblings_maxRSquared=cms.vdouble(
                                                                             0.0006, 0.0006, 0.0006),
                                                                         skip_layers=cms.int32(
                                                                             0),
                                                                         type=cms.string(
                                                                             'CA')
                                                                     ),
                                                                     pluginPatternRecognitionByCLUE3D=cms.PSet(
                                                                         algo_verbosity=cms.int32(
                                                                             0),
                                                                         computeLocalTime=cms.bool(
                                                                             True),
                                                                         criticalDensity=cms.vdouble(
                                                                             0.6, 0.6, 0.6),
                                                                         criticalEtaPhiDistance=cms.vdouble(
                                                                             0.025, 0.025, 0.025),
                                                                         criticalSelfDensity=cms.vdouble(
                                                                             0.15, 0.15, 0.15),
                                                                         criticalXYDistance=cms.vdouble(
                                                                             1.8, 1.8, 1.8),
                                                                         criticalZDistanceLyr=cms.vint32(
                                                                             5, 5, 5),
                                                                         cutHadProb=cms.double(
                                                                             999),
                                                                         densityEtaPhiDistanceSqr=cms.vdouble(
                                                                             0.0008, 0.0008, 0.0008),
                                                                         densityOnSameLayer=cms.bool(
                                                                             False),
                                                                         densitySiblingLayers=cms.vint32(
                                                                             3, 3, 3),
                                                                         densityXYDistanceSqr=cms.vdouble(
                                                                             3.24, 3.24, 3.24),
                                                                         doPidCut=cms.bool(
                                                                             True),
                                                                         kernelDensityFactor=cms.vdouble(
                                                                             0.2, 0.2, 0.2),
                                                                         minNumLayerCluster=cms.vint32(
                                                                             2, 2, 2),
                                                                         nearestHigherOnSameLayer=cms.bool(
                                                                             False),
                                                                         outlierMultiplier=cms.vdouble(
                                                                             2, 2, 2),
                                                                         rescaleDensityByZ=cms.bool(
                                                                             False),
                                                                         type=cms.string(
                                                                             'CLUE3D'),
                                                                         useAbsoluteProjectiveScale=cms.bool(
                                                                             True),
                                                                         useClusterDimensionXY=cms.bool(
                                                                             False),
                                                                         usePCACleaning=cms.bool(
                                                                             True)
                                                                     ),
                                                                     pluginPatternRecognitionByFastJet=cms.PSet(
                                                                         algo_verbosity=cms.int32(
                                                                             0),
                                                                         antikt_radius=cms.double(
                                                                             0.09),
                                                                         computeLocalTime=cms.bool(
                                                                             False),
                                                                         minNumLayerCluster=cms.int32(
                                                                             5),
                                                                         type=cms.string(
                                                                             'FastJet')
                                                                     ),
                                                                     pluginPatternRecognitionByRecovery=cms.PSet(
                                                                         algo_verbosity=cms.int32(
                                                                             0),
                                                                         type=cms.string(
                                                                             'Recovery')
                                                                     ),
                                                                     seeding_regions=cms.InputTag(
                                                                         "ticlSeedingGlobal"),
                                                                     time_layerclusters=cms.InputTag(
                                                                         "hgcalMergeLayerClusters", "timeLayerCluster")
                                                                     )
            )
    setattr(process, 'filteredLayerClustersRecovery', cms.EDProducer('FilteredLayerClustersProducer',
                                                                              LayerClusters=cms.InputTag(
                                                                                  "hgcalMergeLayerClusters"),
                                                                              LayerClustersInputMask=cms.InputTag(
                                                                                  "ticlTrackstersCLUE3D"),
                                                                              algo_number=cms.vint32(
                                                                                  6, 7, 8),
                                                                              clusterFilter=cms.string(
                                                                                  'ClusterFilterBySize'),
                                                                              iteration_label=cms.string(
                                                                                  'CLUE3DHigh'),
                                                                              max_cluster_size=cms.int32(
                                                                                  9999),
                                                                              max_layerId=cms.int32(
                                                                                  9999),
                                                                              mightGet=cms.optional.untracked.vstring,
                                                                              min_cluster_size=cms.int32(
                                                                                  2),
                                                                              min_layerId=cms.int32(
                                                                                  0)
                                                                              )
            )
    setattr(process, 'ticlTrackstersRecovery', cms.EDProducer('TrackstersProducer',
                                                                     detector=cms.string(
                                                                         'HGCAL'),
                                                                     filtered_mask=cms.InputTag(
                                                                         "filteredLayerClustersRecovery", "CLUE3DHigh"),
                                                                     inferenceAlgo=cms.string(
                                                                         'TracksterInferenceByCNNv4'),
                                                                     itername=cms.string(
                                                                         'Recovery'),
                                                                     layer_clusters=cms.InputTag(
                                                                         "hgcalMergeLayerClusters"),
                                                                     layer_clusters_hfnose_tiles=cms.InputTag(
                                                                         "ticlLayerTileHFNose"),
                                                                     layer_clusters_tiles=cms.InputTag(
                                                                         "ticlLayerTileProducer"),
                                                                     mightGet=cms.optional.untracked.vstring,
                                                                     original_mask=cms.InputTag(
                                                                         "hgcalMergeLayerClusters", "InitialLayerClustersMask"),
                                                                     patternRecognitionBy=cms.string(
                                                                         'Recovery'),
                                                                     pluginInferenceAlgoTracksterInferenceByANN=cms.PSet(
                                                                         algo_verbosity=cms.int32(
                                                                             0),
                                                                         type=cms.string(
                                                                             'TracksterInferenceByANN')
                                                                     ),
                                                                     pluginInferenceAlgoTracksterInferenceByCNNv4=cms.PSet(
                                                                         algo_verbosity=cms.int32(
                                                                             0),
                                                                         doPID=cms.int32(
                                                                             1),
                                                                         doRegression=cms.int32(
                                                                             0),
                                                                         eid_min_cluster_energy=cms.double(
                                                                             1),
                                                                         eid_n_clusters=cms.int32(
                                                                             10),
                                                                         eid_n_layers=cms.int32(
                                                                             50),
                                                                         inputNames=cms.vstring(
                                                                             'input:0'),
                                                                         onnxModelPath=cms.FileInPath(
                                                                             'RecoHGCal/TICL/data/ticlv4/onnx_models/energy_id_v0.onnx'),
                                                                         outputNames=cms.vstring(
                                                                             'output/regressed_energy:0',
                                                                             'output/id_probabilities:0'
                                                                         ),
                                                                         type=cms.string(
                                                                             'TracksterInferenceByCNNv4')
                                                                     ),
                                                                     pluginInferenceAlgoTracksterInferenceByDNN=cms.PSet(
                                                                         algo_verbosity=cms.int32(
                                                                             0),
                                                                         doPID=cms.int32(
                                                                             1),
                                                                         doRegression=cms.int32(
                                                                             0),
                                                                         eid_min_cluster_energy=cms.double(
                                                                             1),
                                                                         eid_n_clusters=cms.int32(
                                                                             10),
                                                                         eid_n_layers=cms.int32(
                                                                             50),
                                                                         inputNames=cms.vstring(
                                                                             'input'),
                                                                         onnxEnergyModelPath=cms.FileInPath(
                                                                             'RecoHGCal/TICL/data/ticlv5/onnx_models/patternrecognition/energy_v0.onnx'),
                                                                         onnxPIDModelPath=cms.FileInPath(
                                                                             'RecoHGCal/TICL/data/ticlv5/onnx_models/patternrecognition/id_v0.onnx'),
                                                                         output_en=cms.vstring(
                                                                             'enreg_output'),
                                                                         output_id=cms.vstring(
                                                                             'pid_output'),
                                                                         type=cms.string(
                                                                             'TracksterInferenceByDNN')
                                                                     ),
                                                                     pluginPatternRecognitionByCA=cms.PSet(
                                                                         algo_verbosity=cms.int32(
                                                                             0),
                                                                         computeLocalTime=cms.bool(
                                                                             False),
                                                                         energy_em_over_total_threshold=cms.double(
                                                                             -1),
                                                                         etaLimitIncreaseWindow=cms.double(
                                                                             2.1),
                                                                         filter_on_categories=cms.vint32(
                                                                             0),
                                                                         max_delta_time=cms.double(
                                                                             3),
                                                                         max_longitudinal_sigmaPCA=cms.double(
                                                                             9999),
                                                                         max_missing_layers_in_trackster=cms.int32(
                                                                             9999),
                                                                         max_out_in_hops=cms.int32(
                                                                             10),
                                                                         min_cos_pointing=cms.double(
                                                                             -1),
                                                                         min_cos_theta=cms.double(
                                                                             0.915),
                                                                         min_layers_per_trackster=cms.int32(
                                                                             10),
                                                                         oneTracksterPerTrackSeed=cms.bool(
                                                                             False),
                                                                         out_in_dfs=cms.bool(
                                                                             True),
                                                                         pid_threshold=cms.double(
                                                                             0),
                                                                         promoteEmptyRegionToTrackster=cms.bool(
                                                                             False),
                                                                         root_doublet_max_distance_from_seed_squared=cms.double(
                                                                             9999),
                                                                         shower_start_max_layer=cms.int32(
                                                                             9999),
                                                                         siblings_maxRSquared=cms.vdouble(
                                                                             0.0006, 0.0006, 0.0006),
                                                                         skip_layers=cms.int32(
                                                                             0),
                                                                         type=cms.string(
                                                                             'CA')
                                                                     ),
                                                                     pluginPatternRecognitionByCLUE3D=cms.PSet(
                                                                         algo_verbosity=cms.int32(
                                                                             0),
                                                                         computeLocalTime=cms.bool(
                                                                             True),
                                                                         criticalDensity=cms.vdouble(
                                                                             0.6, 0.6, 0.6),
                                                                         criticalEtaPhiDistance=cms.vdouble(
                                                                             0.025, 0.025, 0.025),
                                                                         criticalSelfDensity=cms.vdouble(
                                                                             0.15, 0.15, 0.15),
                                                                         criticalXYDistance=cms.vdouble(
                                                                             1.8, 1.8, 1.8),
                                                                         criticalZDistanceLyr=cms.vint32(
                                                                             5, 5, 5),
                                                                         cutHadProb=cms.double(
                                                                             999),
                                                                         densityEtaPhiDistanceSqr=cms.vdouble(
                                                                             0.0008, 0.0008, 0.0008),
                                                                         densityOnSameLayer=cms.bool(
                                                                             False),
                                                                         densitySiblingLayers=cms.vint32(
                                                                             3, 3, 3),
                                                                         densityXYDistanceSqr=cms.vdouble(
                                                                             3.24, 3.24, 3.24),
                                                                         doPidCut=cms.bool(
                                                                             True),
                                                                         kernelDensityFactor=cms.vdouble(
                                                                             0.2, 0.2, 0.2),
                                                                         minNumLayerCluster=cms.vint32(
                                                                             2, 2, 2),
                                                                         nearestHigherOnSameLayer=cms.bool(
                                                                             False),
                                                                         outlierMultiplier=cms.vdouble(
                                                                             2, 2, 2),
                                                                         rescaleDensityByZ=cms.bool(
                                                                             False),
                                                                         type=cms.string(
                                                                             'CLUE3D'),
                                                                         useAbsoluteProjectiveScale=cms.bool(
                                                                             True),
                                                                         useClusterDimensionXY=cms.bool(
                                                                             False),
                                                                         usePCACleaning=cms.bool(
                                                                             True)
                                                                     ),
                                                                     pluginPatternRecognitionByFastJet=cms.PSet(
                                                                         algo_verbosity=cms.int32(
                                                                             0),
                                                                         antikt_radius=cms.double(
                                                                             0.09),
                                                                         computeLocalTime=cms.bool(
                                                                             False),
                                                                         minNumLayerCluster=cms.int32(
                                                                             5),
                                                                         type=cms.string(
                                                                             'FastJet')
                                                                     ),
                                                                     pluginPatternRecognitionByRecovery=cms.PSet(
                                                                         algo_verbosity=cms.int32(
                                                                             0),
                                                                         type=cms.string(
                                                                             'Recovery')
                                                                     ),
                                                                     seeding_regions=cms.InputTag(
                                                                         "ticlSeedingGlobal"),
                                                                     time_layerclusters=cms.InputTag(
                                                                         "hgcalMergeLayerClusters", "timeLayerCluster")
                                                                     )
            )
    setattr(process, 'tracksterLinks' + str(i), cms.EDProducer("TracksterLinksProducer",
                                                                         detector=cms.string(
                                                                             'HGCAL'),
                                                                         inferenceAlgo=cms.string(
                                                                             'TracksterInferenceByDNN'),
                                                                         layer_clusters=cms.InputTag(
                                                                             "hgcalMergeLayerClusters"),
                                                                         layer_clustersTime=cms.InputTag(
                                                                             "hgcalMergeLayerClusters", "timeLayerCluster"),
                                                                         linkingPSet=cms.PSet(
                                                                             algo_verbosity=cms.int32(
                                                                                 0),
                                                                             cylinder_radius_sqr=cms.vdouble(
                                                                                 p[0], p[1]),
                                                                             cylinder_radius_sqr_split=cms.double(
                                                                                 p[2]),
                                                                             deltaRxy=cms.double(
                                                                                 p[3]),
                                                                             dot_prod_th=cms.double(
                                                                                 p[4]),
                                                                             lower_boundary=cms.vdouble(
                                                                                 p[5], p[6]),
                                                                             lower_distance_projective_sqr=cms.vdouble(
                                                                                 p[7], p[8]),
                                                                             lower_distance_projective_sqr_closest_points=cms.vdouble(
                                                                                 p[9], p[10]),
                                                                             max_z_distance_closest_points=cms.vdouble(
                                                                                 35, 35),
                                                                             min_distance_z=cms.vdouble(
                                                                                 35, 35),
                                                                             min_num_lcs=cms.uint32(
                                                                                 int(p[11])),
                                                                             min_trackster_energy=cms.double(
                                                                                 p[12]),
                                                                             pca_quality_th=cms.double(
                                                                                 0.85),
                                                                             proj_distance_split=cms.double(
                                                                                 5),
                                                                             track_time_quality_threshold=cms.double(
                                                                                 0.5),
                                                                             type=cms.string(
                                                                                 'Skeletons'),
                                                                             upper_boundary=cms.vdouble(
                                                                                 p[13], p[14]),
                                                                             upper_distance_projective_sqr=cms.vdouble(
                                                                                 p[15], p[16]),
                                                                             upper_distance_projective_sqr_closest_points=cms.vdouble(
                                                                                 5, 30)
                                                                         ),
                                                                         mightGet=cms.optional.untracked.vstring,
                                                                         original_masks=cms.VInputTag(
                                                                             "hgcalMergeLayerClusters:InitialLayerClustersMask"),
                                                                         pluginInferenceAlgoTracksterInferenceByCNNv4=cms.PSet(
                                                                             algo_verbosity=cms.int32(
                                                                                 0),
                                                                             doPID=cms.int32(
                                                                                 1),
                                                                             doRegression=cms.int32(
                                                                                 0),
                                                                             eid_min_cluster_energy=cms.double(
                                                                                 1),
                                                                             eid_n_clusters=cms.int32(
                                                                                 10),
                                                                             eid_n_layers=cms.int32(
                                                                                 50),
                                                                             inputNames=cms.vstring(
                                                                                 'input:0'),
                                                                             onnxModelPath=cms.FileInPath(
                                                                                 'RecoHGCal/TICL/data/ticlv4/onnx_models/energy_id_v0.onnx'),
                                                                             outputNames=cms.vstring(
                                                                                 'output/regressed_energy:0',
                                                                                 'output/id_probabilities:0'
                                                                             ),
                                                                             type=cms.string(
                                                                                 'TracksterInferenceByCNNv4')
                                                                         ),
                                                                         pluginInferenceAlgoTracksterInferenceByDNN=cms.PSet(
                                                                             algo_verbosity=cms.int32(
                                                                                 0),
                                                                             doPID=cms.int32(
                                                                                 1),
                                                                             doRegression=cms.int32(
                                                                                 1),
                                                                             eid_min_cluster_energy=cms.double(
                                                                                 1),
                                                                             eid_n_clusters=cms.int32(
                                                                                 10),
                                                                             eid_n_layers=cms.int32(
                                                                                 50),
                                                                             inputNames=cms.vstring(
                                                                                 'input'),
                                                                             onnxEnergyModelPath=cms.FileInPath(
                                                                                 'RecoHGCal/TICL/data/ticlv5/onnx_models/linking/energy_v0.onnx'),
                                                                             onnxPIDModelPath=cms.FileInPath(
                                                                                 'RecoHGCal/TICL/data/ticlv5/onnx_models/linking/id_v0.onnx'),
                                                                             output_en=cms.vstring(
                                                                                 'enreg_output'),
                                                                             output_id=cms.vstring(
                                                                                 'pid_output'),
                                                                             type=cms.string(
                                                                                 'TracksterInferenceByDNN')
                                                                         ),
                                                                         propagator=cms.string(
                                                                             'PropagatorWithMaterial'),
                                                                         regressionAndPid=cms.bool(
                                                                             True),
                                                                         tracksters_collections=cms.VInputTag(
                                                                             f"ticlTrackstersCLUE3D", f"ticlTrackstersRecovery")
                                                                         )
            )

    setattr(process, f"allHitToTracksterAssociations{str(i)}", cms.EDProducer("AllHitToTracksterAssociatorsProducer",
                                                                                 hitMapTag=cms.InputTag(
                                                                                     "recHitMapProducer", "hgcalRecHitMap"),
                                                                                 hits=cms.VInputTag("HGCalRecHit:HGCEERecHits",
                                                                                                    "HGCalRecHit:HGCHEFRecHits", "HGCalRecHit:HGCHEBRecHits"),
                                                                                 layerClusters=cms.InputTag(
                                                                                     "hgcalMergeLayerClusters"),
                                                                                 mightGet=cms.optional.untracked.vstring,
                                                                                 tracksterCollections=cms.VInputTag(
                                                                                     cms.InputTag(f"tracksterLinks{str(i)}"), 
                                                                                    cms.InputTag(
                                                                                    "ticlSimTracksters"), cms.InputTag("ticlSimTracksters", "fromCPs")
                                                                                 )
                                                                                 )
            )
#    setattr(process, f"hitToSimClusterCaloParticleAssociator", cms.EDProducer("HitToSimClusterCaloParticleAssociatorProducer",
#                                                                              caloParticles=cms.InputTag(
#                                                                                  "mix", "MergedCaloTruth"),
#                                                                              hitMap=cms.InputTag(
#                                                                                  "recHitMapProducer", "hgcalRecHitMap"),
#                                                                              hits=cms.VInputTag("HGCalRecHit:HGCEERecHits",
#                                                                                                 "HGCalRecHit:HGCHEFRecHits", "HGCalRecHit:HGCHEBRecHits"),
#                                                                              mightGet=cms.optional.untracked.vstring,
#                                                                              simClusters=cms.InputTag(
#                                                                                  "mix", "MergedCaloTruth")
#                                                                              )
#            )
    setattr(process, f"allLCToTrackster{str(i)}", cms.EDProducer("AllLayerClusterToTracksterAssociatorsProducer",
                                                                layer_clusters = cms.InputTag('hgcalMergeLayerClusters'),
                                                                tracksterCollections = cms.VInputTag(
                                                                    "tracksterLinks"+str(i),
                                                                    cms.InputTag("ticlSimTracksters"),
                                                                    cms.InputTag("ticlSimTracksters", "fromCPs"),
                                                                )
                                                                              )
            
            )
#    setattr(process, f"allTrackstersToSimTrackstersAssociationsByHits{str(i)}", cms.EDProducer("AllTracksterToSimTracksterAssociatorsByHitsProducer",
#                                                                                                  prefixAssoc = cms.string(f"allHitToTracksterAssociations{str(i)}"),
#                                                                                                  caloParticles=cms.InputTag(
#                                                                                                      "mix", "MergedCaloTruth"),
#                                                                                                  hitToCaloParticleMap=cms.InputTag(
#                                                                                                      "hitToSimClusterCaloParticleAssociator", "hitToCaloParticleMap"),
#                                                                                                  hitToSimClusterMap=cms.InputTag(
#                                                                                                      "hitToSimClusterCaloParticleAssociator", "hitToSimClusterMap"),
#                                                                                                  hits=cms.VInputTag("HGCalRecHit:HGCEERecHits",
#                                                                                                                     "HGCalRecHit:HGCHEFRecHits", "HGCalRecHit:HGCHEBRecHits"),
#                                                                                                  mightGet=cms.optional.untracked.vstring,
#                                                                                                  simTracksterCollections=cms.VInputTag(
#                                                                                                      "ticlSimTracksters", "ticlSimTracksters:fromCPs"),
#                                                                                                  tracksterCollections=cms.VInputTag(
#                                                                                                      cms.InputTag(f"tracksterLinks{str(i)}"))
#                                                                                                  )
#            )
    setattr(process, f"allTrackstersToSimTrackstersAssociationsByLCs{str(i)}", cms.EDProducer("AllTracksterToSimTracksterAssociatorsByLCsProducer",
                                                                                                  prefixAssoc = cms.string(f"allLCToTrackster{str(i)}"),
                                                                                                  layerClusters = cms.InputTag('hgcalMergeLayerClusters'),
                                                                                                  simTracksterCollections=cms.VInputTag(
                                                                                                      "ticlSimTracksters", "ticlSimTracksters:fromCPs"),
                                                                                                  tracksterCollections=cms.VInputTag(
                                                                                                      cms.InputTag(f"tracksterLinks{str(i)}"))
                                                                                                  )
            )
    setattr(process, f"ticlDumperSimple{str(i)}", cms.EDAnalyzer("TICLDumperSimple",
                                                                    associators=cms.VPSet(
                                                                        cms.PSet(
                                                                            associatorRecoToSimInputTag=cms.InputTag(
                                                                                f"allTrackstersToSimTrackstersAssociationsByLCs{str(i)}", f"tracksterLinks{str(i)}ToticlSimTracksters"),
                                                                            associatorSimToRecoInputTag=cms.InputTag(
                                                                                f"allTrackstersToSimTrackstersAssociationsByLCs{str(i)}", f"ticlSimTrackstersTotracksterLinks{str(i)}"),

                                                                            branchName=cms.string(
                                                                                'ticlTracksterLinks'),
                                                                            suffix=cms.string(
                                                                                'SC')
                                                                        ),
                                                                        cms.PSet(
                                                                            associatorRecoToSimInputTag=cms.InputTag(
                                                                                f"allTrackstersToSimTrackstersAssociationsByLCs{str(i)}", f"tracksterLinks{str(i)}ToticlSimTrackstersfromCPs"),
                                                                            associatorSimToRecoInputTag=cms.InputTag(
                                                                                f"allTrackstersToSimTrackstersAssociationsByLCs{str(i)}", f"ticlSimTrackstersfromCPsTotracksterLinks{str(i)}"),
                                                                            branchName=cms.string(
                                                                                'ticlTracksterLinks'),
                                                                            suffix=cms.string(
                                                                                'CP')
                                                                        ),
                                                                    ),
                                                                    caloparticles=cms.InputTag(
                                                                        "mix", "MergedCaloTruth"),
                                                                    mightGet=cms.optional.untracked.vstring,
                                                                    saveLCs=cms.bool(
                                                                        False),
                                                                    saveRecoSuperclusters=cms.bool(
                                                                        False),
                                                                    saveSimTICLCandidate=cms.bool(
                                                                        False),
                                                                    saveSuperclustering=cms.bool(
                                                                        False),
                                                                    saveTICLCandidate=cms.bool(
                                                                        False),
                                                                    saveTracks=cms.bool(
                                                                        False),
                                                                    simclusters=cms.InputTag(
                                                                        "mix", "MergedCaloTruth"),
                                                                    tracksterCollections=cms.VPSet(
                                                                        cms.PSet(
                                                                            inputTag=cms.InputTag(
                                                                                f"tracksterLinks{str(i)}"),
                                                                            treeName=cms.string(
                                                                                'ticlTracksterLinks')
                                                                        ),
                                                                        cms.PSet(
                                                                            inputTag=cms.InputTag(
                                                                                "ticlSimTracksters"),
                                                                            tracksterType=cms.string(
                                                                                'SimTracksterSC'),
                                                                            treeName=cms.string(
                                                                                'simtrackstersSC')
                                                                        ),
                                                                        cms.PSet(
                                                                            inputTag=cms.InputTag(
                                                                                "ticlSimTracksters", "fromCPs"),
                                                                            tracksterType=cms.string(
                                                                                'SimTracksterCP'),
                                                                            treeName=cms.string(
                                                                                'simtrackstersCP')
                                                                        )
                                                                    )
                                                                    )
            )

#    setattr(process, "simpleValidation" + str(i), cms.EDAnalyzer('SimpleValidation',
#      trackstersclue3d=cms.InputTag('ticlTrackstersCLUE3D' + str(i)),
#      trackstersMerged=cms.InputTag('tracksterLinks' + str(i)),
#      simtrackstersCP=cms.InputTag('ticlSimTracksters', 'fromCPs'),
#      caloParticles=cms.InputTag('mix', 'MergedCaloTruth'),
#      layerClusters=cms.InputTag('hgcalMergeLayerClusters'),
#      recoToSimAssociatorCP=cms.InputTag(
#          'tracksterSimTracksterFromCPsAssociationPR'+str(i), 'tracksterToSimTracksterMap'),
#      simToRecoAssociatorCP=cms.InputTag(
#          'tracksterSimTracksterFromCPsAssociationPR'+str(i), 'simTracksterToTracksterMap'),
#      MergerecoToSimAssociatorCP=cms.InputTag(
#          'tracksterSimTracksterFromCPsAssociationLinking'+str(i), 'tracksterToSimTracksterMap'),
#      MergesimToRecoAssociatorCP=cms.InputTag(
#          'tracksterSimTracksterFromCPsAssociationLinking'+str(i), 'simTracksterToTracksterMap'),
#      mightGet=cms.optional.untracked.vstring
#    )
#    )


taskListTrackstersCLUE3D = []
taskListTrackstersCLUE3D.extend(
    [getattr(process, 'ticlTrackstersCLUE3D')])
taskListTrackstersCLUE3D.extend([getattr(
    process, 'filteredLayerClustersRecovery')])
taskListTrackstersCLUE3D.extend(
    [getattr(process, 'ticlTrackstersRecovery')])

tracksterLinksTasks = []
tracksterLinksTasks.extend(
    [getattr(process, 'tracksterLinks' + str(i)) for i in range(len(params))])

#TaskAssociations = [getattr(process, 'recHitMapProducer'), getattr(process, 'hitToSimClusterCaloParticleAssociator')]
TaskAssociations = [] 
TaskAssociations.extend([getattr(process, "allLCToTrackster" + str(i)) for i in range(len(params))])
TaskAssociations.extend([getattr(
    process, 'allTrackstersToSimTrackstersAssociationsByLCs' + str(i)) for i in range(len(params))])
taskSimpleValidation = [getattr(process, 'ticlDumperSimple' + str(i))
                        for i in range(len(params))]

process.TFESSource = cms.Task(process.trackdnn_source)
process.hgcalLayerClustersTask = cms.Task(process.hgcalLayerClustersEE,
                                          process.hgcalLayerClustersHSi,
                                          process.hgcalLayerClustersHSci,
                                          process.hgcalMergeLayerClusters)

process.trackstersProducersTask = cms.Task(
    process.ticlSeedingRegionProducer, process.filteredLayerClustersCLUE3DHigh, *taskListTrackstersCLUE3D, *tracksterLinksTasks)
# process.Tracer = cms.Service('Tracer')
process.TFileService = cms.Service('TFileService', fileName=cms.string(options.outputFile)
                                   if cms.string(options.outputFile) else 'default.root')
# Path and EndPath definitions
process.TICL = cms.Path(process.TFESSource,
                        process.ticlLayerTileTask,
                        process.trackstersProducersTask)

process.TICLAssociators = cms.Task(*TaskAssociations)
process.TICLValidation = cms.Path(process.TICLAssociators)
process.consume_step = cms.EndPath()
for t in taskSimpleValidation:
    process.consume_step += t

process.schedule = cms.Schedule(
    process.TICL, process.TICLValidation,  process.consume_step)

process.options.wantSummary = True
process.options.numberOfThreads =  4 
process.options.numberOfStreams = 0 

process = setCrossingFrameOn(process)


# Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
process = customiseEarlyDelete(process)
# End adding early deletion
