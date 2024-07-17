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
from RecoHGCal.TICL.ticlDumper_cfi import ticlDumper
# Validation
from Validation.HGCalValidation.HGCalValidator_cfi import *
#from RecoLocalCalo.HGCalRecProducers.hgcalRecHitMapProducer_cfi import hgcalRecHitMapProducer

# Load DNN ESSource
from RecoTracker.IterativeTracking.iterativeTk_cff import trackdnn_source

# Automatic addition of the customisation function from RecoHGCal.Configuration.RecoHGCal_EventContent_cff
from RecoHGCal.Configuration.RecoHGCal_EventContent_cff import customiseHGCalOnlyEventContent
from SimCalorimetry.HGCalAssociatorProducers.simTracksterAssociatorByEnergyScore_cfi import simTracksterAssociatorByEnergyScore as simTsAssocByEnergyScoreProducer
#from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import 
#from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociationByHits_cfi import tracksterSimTracksterFromCPsAssociationLinking, tracksterSimTracksterAssociationLinking, tracksterSimTracksterFromCPsAssociationPR, tracksterSimTracksterAssociationPR
from RecoHGCal.TICL.TICLSeedingRegions_cff import ticlSeedingGlobal
from RecoHGCal.TICL.simpleValidation_cfi import simpleValidation

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
from FWCore.ParameterSet.VarParsing import VarParsing
from utils import read_csv

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

process = cms.Process('RECO3', Phase2C17I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D95Reco_cff')
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
                    'file:/data/user/wredjeb/Optimizer/subSet/step3_385_15.root',
                    'file:/data/user/wredjeb/Optimizer/subSet/step3_385_26.root',
                    'file:/data/user/wredjeb/Optimizer/subSet/step3_385_22.root',
                    'file:/data/user/wredjeb/Optimizer/subSet/step3_385_19.root',
                    'file:/data/user/wredjeb/Optimizer/subSet/step3_385_11.root',
                    'file:/data/user/wredjeb/Optimizer/subSet/step3_385_3.root',
                    'file:/data/user/wredjeb/Optimizer/subSet/step3_385_33.root'
                    ]
                            ),
                            secondaryFileNames=cms.untracked.vstring()
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
    process.GlobalTag, 'auto:phase2_realistic_T21', '')

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
    LayerClusters=cms.InputTag("hgcalMergeLayerClusters"),
    LayerClustersInputMask=cms.InputTag(
        "hgcalMergeLayerClusters", "InitialLayerClustersMask"),
    algo_number=cms.vint32(6, 7, 8),
    clusterFilter=cms.string('ClusterFilterByAlgoAndSize'),
    iteration_label=cms.string('CLUE3DEM'),
    max_cluster_size=cms.int32(9999),
    max_layerId=cms.int32(9999),
    mightGet=cms.optional.untracked.vstring,
    min_cluster_size=cms.int32(2),
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
# params = [[i for i in range(100)]]
totalTask = len(params)
print(params)
for i, p in enumerate(params):
    setattr(process, 'ticlTrackstersCLUE3D' + str(i), cms.EDProducer('TrackstersProducer',
    detector=cms.string('HGCAL'),
    layer_clusters=cms.InputTag('hgcalMergeLayerClusters'),
    filtered_mask=cms.InputTag('filteredLayerClustersCLUE3DHigh', 'CLUE3DHigh'),
    original_mask=cms.InputTag(
        'hgcalMergeLayerClusters', 'InitialLayerClustersMask'),
    time_layerclusters=cms.InputTag(
        'hgcalMergeLayerClusters', 'timeLayerCluster'),
    layer_clusters_tiles=cms.InputTag('ticlLayerTileProducer'),
    layer_clusters_hfnose_tiles=cms.InputTag('ticlLayerTileHFNose'),
    seeding_regions=cms.InputTag('ticlSeedingGlobal'),
    patternRecognitionBy=cms.string('CLUE3D'),
    itername=cms.string('unknown'),
    tfDnnLabel=cms.string('tracksterSelectionTf'),
    pluginPatternRecognitionByCA=cms.PSet(
      algo_verbosity=cms.int32(0),
      oneTracksterPerTrackSeed=cms.bool(False),
      promoteEmptyRegionToTrackster=cms.bool(False),
      out_in_dfs=cms.bool(True),
      max_out_in_hops=cms.int32(10),
      min_cos_theta=cms.double(0.915),
      min_cos_pointing=cms.double(-1),
      root_doublet_max_distance_from_seed_squared=cms.double(9999),
      etaLimitIncreaseWindow=cms.double(2.1),
      skip_layers=cms.int32(0),
      max_missing_layers_in_trackster=cms.int32(9999),
      shower_start_max_layer=cms.int32(9999),
      min_layers_per_trackster=cms.int32(10),
      filter_on_categories=cms.vint32(0),
      pid_threshold=cms.double(0),
      energy_em_over_total_threshold=cms.double(-1),
      max_longitudinal_sigmaPCA=cms.double(9999),
      max_delta_time=cms.double(3),
      eid_input_name=cms.string('input'),
      eid_output_name_energy=cms.string('output/regressed_energy'),
      eid_output_name_id=cms.string('output/id_probabilities'),
      eid_min_cluster_energy=cms.double(1),
      eid_n_layers=cms.int32(50),
      eid_n_clusters=cms.int32(10),
      computeLocalTime=cms.bool(False),
      siblings_maxRSquared=cms.vdouble(
        0.0006,
        0.0006,
        0.0006
      ),
      type=cms.string('CA')

    ),
    pluginPatternRecognitionByCLUE3D=cms.PSet(
          algo_verbosity=cms.int32(0),
          criticalDensity=cms.vdouble(v[0], v[1], v[2]),
          criticalEtaPhiDistance=cms.vdouble(0.025, 0.025, 0.025),
          criticalSelfDensity=cms.vdouble(v[3], v[4], v[5]),
          criticalXYDistance=cms.vdouble(v[6], v[7], v[8]),
          criticalZDistanceLyr=cms.vint32(int(v[9]), int(v[10]), int(v[11])),
          cutHadProb=cms.double(0.5),
          densityEtaPhiDistanceSqr=cms.vdouble(0.0008, 0.0008, 0.0008),
          densityOnSameLayer=cms.bool(False),
          densitySiblingLayers=cms.vint32(int(v[12]), int(v[13]), int(v[14])),
          densityXYDistanceSqr=cms.vdouble(v[15], v[16], v[17]),
          doPidCut=cms.bool(False),
          eid_input_name=cms.string('input'),
          eid_min_cluster_energy=cms.double(1),
          eid_n_clusters=cms.int32(10),
          eid_n_layers=cms.int32(50),
          eid_output_name_energy=cms.string('output/regressed_energy'),
          eid_output_name_id=cms.string('output/id_probabilities'),
          kernelDensityFactor=cms.vdouble(v[18], v[19], v[20]),
          minNumLayerCluster=cms.vint32(int(v[21]), int(v[22]), int(v[23])),
          nearestHigherOnSameLayer=cms.bool(False),
          outlierMultiplier=cms.vdouble(v[24], v[25], v[26]),
          rescaleDensityByZ=cms.bool(False),
          type=cms.string('CLUE3D'),
          useAbsoluteProjectiveScale=cms.bool(True),
          useClusterDimensionXY=cms.bool(False)
      ),
    pluginPatternRecognitionByFastJet=cms.PSet(
      algo_verbosity=cms.int32(0),
      antikt_radius=cms.double(0.09),
      minNumLayerCluster=cms.int32(5),
      eid_input_name=cms.string('input'),
      eid_output_name_energy=cms.string('output/regressed_energy'),
      eid_output_name_id=cms.string('output/id_probabilities'),
      eid_min_cluster_energy=cms.double(1),
      eid_n_layers=cms.int32(50),
      eid_n_clusters=cms.int32(10),
      computeLocalTime=cms.bool(False),
      type=cms.string('FastJet')

    ),
    pluginPatternRecognitionByPassthrough=cms.PSet(
      algo_verbosity=cms.int32(0),
      type=cms.string('Passthrough')

    ),
    mightGet=cms.optional.untracked.vstring
    )
  )
    setattr(process, 'filteredLayerClustersPassThrough' + str(i), cms.EDProducer('FilteredLayerClustersProducer',
      LayerClusters = cms.InputTag("hgcalMergeLayerClusters"),
      LayerClustersInputMask = cms.InputTag("ticlTrackstersCLUE3D" + str(i)),
      algo_number = cms.vint32(6, 7, 8),
      clusterFilter = cms.string('ClusterFilterBySize'),
      iteration_label = cms.string('Passthrough'),
      max_cluster_size = cms.int32(9999),
      max_layerId = cms.int32(9999),
      mightGet = cms.optional.untracked.vstring,
      min_cluster_size = cms.int32(2),
      min_layerId = cms.int32(0)
     )	
    )
    setattr(process, 'ticlTrackstersPassThrough' + str(i), cms.EDProducer('TrackstersProducer',
      detector = cms.string('HGCAL'),
      layer_clusters = cms.InputTag('hgcalMergeLayerClusters'),
      filtered_mask = cms.InputTag('filteredLayerClustersPassThrough'+str(i), 'Passthrough'),
      original_mask = cms.InputTag('ticlTrackstersCLUE3D'+str(i)),
      time_layerclusters = cms.InputTag('hgcalMergeLayerClusters', 'timeLayerCluster'),
      layer_clusters_tiles = cms.InputTag('ticlLayerTileProducer'),
      layer_clusters_hfnose_tiles = cms.InputTag('ticlLayerTileHFNose'),
      seeding_regions = cms.InputTag('ticlSeedingGlobal'),
      patternRecognitionBy = cms.string('Passthrough'),
      itername = cms.string('unknown'),
      tfDnnLabel = cms.string('tracksterSelectionTf'),
      pluginPatternRecognitionByCA = cms.PSet(
        algo_verbosity = cms.int32(0),
        oneTracksterPerTrackSeed = cms.bool(False),
        promoteEmptyRegionToTrackster = cms.bool(False),
        out_in_dfs = cms.bool(True),
        max_out_in_hops = cms.int32(10),
        min_cos_theta = cms.double(0.915),
        min_cos_pointing = cms.double(-1),
        root_doublet_max_distance_from_seed_squared = cms.double(9999),
        etaLimitIncreaseWindow = cms.double(2.1),
        skip_layers = cms.int32(0),
        max_missing_layers_in_trackster = cms.int32(9999),
        shower_start_max_layer = cms.int32(9999),
        min_layers_per_trackster = cms.int32(10),
        filter_on_categories = cms.vint32(0),
        pid_threshold = cms.double(0),
        energy_em_over_total_threshold = cms.double(-1),
        max_longitudinal_sigmaPCA = cms.double(9999),
        max_delta_time = cms.double(3),
        eid_input_name = cms.string('input'),
        eid_output_name_energy = cms.string('output/regressed_energy'),
        eid_output_name_id = cms.string('output/id_probabilities'),
        eid_min_cluster_energy = cms.double(1),
        eid_n_layers = cms.int32(50),
        eid_n_clusters = cms.int32(10),
        computeLocalTime = cms.bool(False),
        siblings_maxRSquared = cms.vdouble(
          0.0006,
          0.0006,
          0.0006
        ),
        type = cms.string('CA')
      ),
      pluginPatternRecognitionByCLUE3D = cms.PSet(
            algo_verbosity = cms.int32(0),
            criticalDensity = cms.vdouble(0.6, v[0], v[1]),
            criticalEtaPhiDistance = cms.vdouble(0.025, 0.025, 0.025),
            criticalSelfDensity = cms.vdouble(0.15, v[2], v[3]),
            criticalXYDistance = cms.vdouble(1.8, v[4], v[5]),
            criticalZDistanceLyr = cms.vint32(5, int(v[6]), int(v[7])),
            cutHadProb = cms.double(0.5),
            densityEtaPhiDistanceSqr = cms.vdouble(0.0008, 0.0008, 0.0008),
            densityOnSameLayer = cms.bool(False),
            densitySiblingLayers = cms.vint32(3, int(v[8]), int(v[9])),
            densityXYDistanceSqr = cms.vdouble(3.24, v[10], v[11]),
            doPidCut = cms.bool(False),
            eid_input_name = cms.string('input'),
            eid_min_cluster_energy = cms.double(1),
            eid_n_clusters = cms.int32(10),
            eid_n_layers = cms.int32(50),
            eid_output_name_energy = cms.string('output/regressed_energy'),
            eid_output_name_id = cms.string('output/id_probabilities'),
            kernelDensityFactor = cms.vdouble(0.2, v[12], v[13]),
            minNumLayerCluster = cms.vint32(2, int(v[14]), int(v[15])),
            nearestHigherOnSameLayer = cms.bool(False),
            outlierMultiplier = cms.vdouble(2, v[16], v[17]),
            rescaleDensityByZ = cms.bool(False),
            type = cms.string('CLUE3D'),
            useAbsoluteProjectiveScale = cms.bool(True),
            useClusterDimensionXY = cms.bool(False)
        ),
      pluginPatternRecognitionByFastJet = cms.PSet(
        algo_verbosity = cms.int32(0),
        antikt_radius = cms.double(0.09),
        minNumLayerCluster = cms.int32(5),
        eid_input_name = cms.string('input'),
        eid_output_name_energy = cms.string('output/regressed_energy'),
        eid_output_name_id = cms.string('output/id_probabilities'),
        eid_min_cluster_energy = cms.double(1),
        eid_n_layers = cms.int32(50),
        eid_n_clusters = cms.int32(10),
        computeLocalTime = cms.bool(False),
        type = cms.string('FastJet')
      
      ),
      pluginPatternRecognitionByPassthrough = cms.PSet(
        algo_verbosity = cms.int32(0),
        type = cms.string('Passthrough')
      
      ),
      mightGet = cms.optional.untracked.vstring
    )
    )
    setattr(process, 'mergedTrackstersProducer' + str(i), cms.EDProducer("TracksterLinksProducer",
          linkingPSet = cms.PSet(
            track_time_quality_threshold = cms.double(0.5),
            wind = cms.double(4),
            min_num_lcs = cms.uint32(7),
            min_trackster_energy = cms.double(10),
            pca_quality_th = cms.double(0.85),
            dot_prod_th = cms.double(0.97),
            max_distance_projective_sqr = cms.vdouble(
              9,
              40
            ),
            min_distance_z = cms.vdouble(
              30,
              30
            ),
            max_distance_projective_sqr_closest_points = cms.vdouble(
              9,
              40
            ),
            max_z_distance_closest_points = cms.vdouble(
              35,
              35
            ),
            cylinder_radius_sqr = cms.vdouble(
              4,
              9
            ),
            algo_verbosity = cms.int32(0),
            type = cms.string('Skeletons')
          
          ),
          tracksters_collections = cms.VInputTag('ticlTrackstersCLUE3D' + str(i), 'ticlTrackstersPassThrough' + str(i)),
          original_masks = cms.VInputTag('hgcalMergeLayerClusters:InitialLayerClustersMask'),
          layer_clusters = cms.InputTag('hgcalMergeLayerClusters'),
          layer_clustersTime = cms.InputTag('hgcalMergeLayerClusters', 'timeLayerCluster'),
          regressionAndPid = cms.bool(False),
          tfDnnLabel = cms.string('tracksterSelectionTf'),
          eid_input_name = cms.string('input'),
          eid_output_name_energy = cms.string('output/regressed_energy'),
          eid_output_name_id = cms.string('output/id_probabilities'),
          eid_min_cluster_energy = cms.double(2.5),
          eid_n_layers = cms.int32(50),
          eid_n_clusters = cms.int32(10),
          detector = cms.string('HGCAL'),
          propagator = cms.string('PropagatorWithMaterial'),
          mightGet = cms.optional.untracked.vstring
       )
      )
    setattr(process, "layerClusterToSimTracksterAssociation", cms.EDProducer("LCToTSAssociatorProducer",
        layer_clusters = cms.InputTag("hgcalMergeLayerClusters"),
        mightGet = cms.optional.untracked.vstring,
        tracksters = cms.InputTag("ticlSimTracksters")
        )
    )
    setattr(process, "layerClusterToTracksterMergeAssociation" + str(i), cms.EDProducer("LCToTSAssociatorProducer",
        layer_clusters = cms.InputTag("hgcalMergeLayerClusters"),
        mightGet = cms.optional.untracked.vstring,
        tracksters = cms.InputTag("mergedTrackstersProducer" + str(i))
        )
    )

    setattr(process, "layerClusterToCLUE3DTracksterAssociation" + str(i), cms.EDProducer("LCToTSAssociatorProducer",
        layer_clusters = cms.InputTag("hgcalMergeLayerClusters"),
        mightGet = cms.optional.untracked.vstring,
        tracksters = cms.InputTag("ticlTrackstersCLUE3D" + str(i))
        )
    )

    setattr(process, "tracksterSimTracksterFromCPsAssociationPR" + str(i), cms.EDProducer("TracksterToSimTracksterAssociatorProducer",
      	layerClusters = cms.InputTag("hgcalMergeLayerClusters"),
      	mightGet = cms.optional.untracked.vstring,
      	simTracksterMap = cms.InputTag("layerClusterToSimTracksterAssociation"),
      	simTracksters = cms.InputTag("ticlSimTracksters"),
      	tracksterMap = cms.InputTag("layerClusterToTracksterMergeAssociation" + str(i)),
      	tracksters = cms.InputTag("ticlTrackstersCLUE3D" + str(i))
      	)
      )
    setattr(process, "tracksterSimTracksterFromCPsAssociationLinking" + str(i), cms.EDProducer("TracksterToSimTracksterAssociatorProducer",
      	layerClusters = cms.InputTag("hgcalMergeLayerClusters"),
      	mightGet = cms.optional.untracked.vstring,
      	simTracksterMap = cms.InputTag("layerClusterToSimTracksterAssociation"),
      	simTracksters = cms.InputTag("ticlSimTracksters"),
      	tracksterMap = cms.InputTag("layerClusterToTracksterMergeAssociation" + str(i)),
      	tracksters = cms.InputTag("mergedTrackstersProducer" + str(i))
      	)
      )
    setattr(process, "simpleValidation" + str(i), cms.EDAnalyzer('SimpleValidation',
      trackstersclue3d = cms.InputTag('ticlTrackstersCLUE3D' + str(i)),
      trackstersMerged = cms.InputTag('mergedTrackstersProducer' + str(i)),
      simtrackstersCP = cms.InputTag('ticlSimTracksters', 'fromCPs'),
      caloParticles = cms.InputTag('mix', 'MergedCaloTruth'),
      layerClusters = cms.InputTag('hgcalMergeLayerClusters'),
      recoToSimAssociatorCP = cms.InputTag('tracksterSimTracksterFromCPsAssociationPR'+str(i), 'tracksterToSimTracksterMap'),
      simToRecoAssociatorCP = cms.InputTag('tracksterSimTracksterFromCPsAssociationPR'+str(i), 'simTracksterToTracksterMap'),
      MergerecoToSimAssociatorCP = cms.InputTag('tracksterSimTracksterFromCPsAssociationLinking'+str(i), 'tracksterToSimTracksterMap'),
      MergesimToRecoAssociatorCP = cms.InputTag('tracksterSimTracksterFromCPsAssociationLinking'+str(i), 'simTracksterToTracksterMap'),
      mightGet = cms.optional.untracked.vstring
    )
    )

    
taskListTrackstersCLUE3D = [] 
taskListTrackstersCLUE3D.extend([getattr(process, 'ticlTrackstersCLUE3D' + str(i)) for i in range(len(params))])
taskListTrackstersCLUE3D.extend([getattr(process, 'filteredLayerClustersPassThrough' + str(i)) for i in range(len(params))])
taskListTrackstersCLUE3D.extend([getattr(process, 'ticlTrackstersPassThrough' + str(i)) for i in range(len(params))])

mergedTrackstersProducerTasks = []
mergedTrackstersProducerTasks.extend([getattr(process, 'mergedTrackstersProducer' + str(i)) for i in range(len(params))])

TaskAssociations = [getattr(process, 'layerClusterToSimTracksterAssociation')]
TaskAssociations.extend([getattr(process, 'tracksterSimTracksterFromCPsAssociationLinking' + str(i)) for i in range(len(params))])
TaskAssociations.extend([getattr(process, 'layerClusterToCLUE3DTracksterAssociation' + str(i)) for i in range(len(params))])
TaskAssociations.extend([getattr(process, 'layerClusterToTracksterMergeAssociation' + str(i)) for i in range(len(params))])
TaskAssociations.extend([getattr(process, 'tracksterSimTracksterFromCPsAssociationPR' + str(i)) for i in range(len(params))])
taskSimpleValidation = [getattr(process, 'simpleValidation' + str(i)) for i in range(len(params))]

process.TFESSource = cms.Task(process.trackdnn_source)
process.hgcalLayerClustersTask = cms.Task(process.hgcalLayerClustersEE,
                                          process.hgcalLayerClustersHSi,
                                          process.hgcalLayerClustersHSci,
                                          process.hgcalMergeLayerClusters)

process.trackstersProducersTask = cms.Task(process.ticlSeedingRegionProducer, process.filteredLayerClustersCLUE3DHigh, *taskListTrackstersCLUE3D, *mergedTrackstersProducerTasks)
#process.Tracer = cms.Service('Tracer')
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

process.schedule = cms.Schedule(process.TICL, process.TICLValidation,  process.consume_step)

process.options.wantSummary = True
process.options.numberOfThreads =  16 
process.options.numberOfStreams = 0 

process=setCrossingFrameOn(process)


# Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
process=customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
process=customiseEarlyDelete(process)
# End adding early deletion
