import numpy as np
import uproot

# calculate the metrics from validation results
def get_metrics(uproot_file, id):
    tree = uproot_file['ticlDumperSimple' + str(id)]
    associations = tree['associations'].arrays(tree['associations'].keys())
    tracksterLinks = tree['ticlTracksterLinks'].arrays(tree['ticlTracksterLinks'].keys())
    simTracksters = tree['simtrackstersCP'].arrays(tree['simtrackstersCP'].keys())

    # Extracting relevant association maps
    sim_to_reco_index = associations['ticlTracksterLinks_simToReco_CP']
    sim_to_reco_score = associations['ticlTracksterLinks_simToReco_CP_score']
    sim_to_reco_sharedE = associations['ticlTracksterLinks_simToReco_CP_sharedE']
    reco_to_sim_index = associations['ticlTracksterLinks_recoToSim_CP']
    reco_to_sim_score = associations['ticlTracksterLinks_recoToSim_CP_score']
    reco_to_sim_sharedE = associations['ticlTracksterLinks_recoToSim_CP_sharedE']

    num_events = len(sim_to_reco_index)
    efficient_sim_tracksters = []
    fake_reco_tracksters = []

    for event in range(num_events):
        simT = simTracksters[event]
        tLinks = tracksterLinks[event]
        sim_tracksters = sim_to_reco_index[event]
        sim_scores = sim_to_reco_score[event]
        sim_sharedE = sim_to_reco_sharedE[event]
        reco_tracksters = reco_to_sim_index[event]
        reco_scores = reco_to_sim_score[event]
        reco_sharedE = reco_to_sim_sharedE[event]

        # Efficiency Calculation
        for sim_idx in range(len(sim_tracksters)):
            if len(sim_scores[sim_idx]) > 0:
                maxSE = np.max(sim_sharedE[sim_idx])
                norm_sharedE = maxSE / np.sum(sim_sharedE[sim_idx])

                if np.sum(norm_sharedE > 0.5) > 0:
                    print(norm_sharedE, maxSE, np.sum(sim_sharedE[sim_idx])) 
                    efficient_sim_tracksters.append(sim_idx)

        # Fake Calculation
        for reco_idx in range(len(reco_tracksters)):
            sim_associations = reco_to_sim_score[reco_idx]

            # Check if this reco_idx appears in any SimToReco map with more than 1 score < 0.2
            low_simToReco_count = 0
            for sim_idx, reco_list in enumerate(sim_to_reco_index[event]):
                if reco_idx in reco_list:
                    reco_positions = np.where(reco_list == reco_idx)[0]  # Get positions
                    low_simToReco_count += np.sum(np.array(sim_to_reco_score[event][sim_idx])[reco_positions] < 0.2)

            has_low_simToReco = low_simToReco_count >= 2  # Must have at least 2 such scores

            # Check if this reco_trackster has any RecoToSim score < 0.6
            has_low_recoToSim = np.any(np.array(sim_associations) < 0.6)

            if not has_low_simToReco and not has_low_recoToSim:
                fake_reco_tracksters.append(reco_idx)

    efficiency = len(efficient_sim_tracksters) / len(sim_to_reco_index) if len(sim_to_reco_index) > 0 else 0
    fake_rate = len(fake_reco_tracksters) / len(reco_to_sim_index) if len(reco_to_sim_index) > 0 else 0

    return {
        'efficiency': efficiency,
        'fake_rate': fake_rate,
        'efficient_sim_tracksters': efficient_sim_tracksters,
        'fake_reco_tracksters': fake_reco_tracksters
    }

# read a csv file, return a matrix
def read_csv(filename):
    matrix = np.genfromtxt(filename, delimiter=",", dtype=float)
    if matrix.ndim == 2:
        return np.genfromtxt(filename, delimiter=",", dtype=float)
    return np.array([matrix])
    
# write a matrix to a csv file
def write_csv(filename, matrix):
    np.savetxt(filename, matrix, fmt='%.18f', delimiter=',')



fileIn = "test.root"
fileUproot = uproot.open(fileIn)
for i in range(0, 4):
    print(get_metrics(fileUproot, i))
