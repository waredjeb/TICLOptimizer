import numpy as np
import uproot

# calculate the metrics from validation results
import numpy as np


def get_metrics(uproot_file, id, num_bins=10):
    fileName = f"ticlDumperSimple{id}"
    tree = uproot_file[fileName]
    associations = tree['associations'].arrays(tree['associations'].keys())
    tracksterLinks = tree['ticlTracksterLinks'].arrays(
        tree['ticlTracksterLinks'].keys())
    simTracksters = tree['simtrackstersSC'].arrays(
        tree['simtrackstersSC'].keys())

    # Extracting relevant association maps
    sim_to_reco_index = associations['ticlTracksterLinks_simToReco_SC']
    sim_to_reco_score = associations['ticlTracksterLinks_simToReco_SC_score']
    sim_to_reco_sharedE = associations['ticlTracksterLinks_simToReco_SC_sharedE']
    reco_to_sim_index = associations['ticlTracksterLinks_recoToSim_SC']
    reco_to_sim_score = associations['ticlTracksterLinks_recoToSim_SC_score']
    reco_to_sim_sharedE = associations['ticlTracksterLinks_recoToSim_SC_sharedE']

    # Define energy bins
    sim_energies = np.concatenate(
        simTracksters['raw_energy'])  # Flatten all energies
    min_energy, max_energy = np.min(sim_energies), np.max(sim_energies)
    bin_edges = np.linspace(min_energy, max_energy,
                            num_bins + 1)  # Equal-width bins

    def bin_indices(e): return np.clip(
        np.digitize(e, bin_edges) - 1, 0, num_bins - 1)

    # Storage for binned metrics
    efficiency_bins = np.zeros(num_bins)
    fake_rate_bins = np.zeros(num_bins)
    sim_bin_counts = np.zeros(num_bins)  # Total sim tracksters per bin
    reco_bin_counts = np.zeros(num_bins)  # Total reco tracksters per bin

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
        # Get bin indices for SimTracksters
        sim_energy_bins = bin_indices(simT['raw_energy'])

        # Efficiency Calculation
        for sim_idx in range(len(sim_tracksters)):
            # Get energy bin for this SimTrackster
            sim_bin = sim_energy_bins[sim_idx]

            sim_bin_counts[sim_bin] += 1  # Count total SimTracksters per bin

            if len(sim_scores[sim_idx]) > 0:
                maxSE = np.max(sim_sharedE[sim_idx])
                norm_sharedE = maxSE / np.sum(sim_sharedE[sim_idx])

                if norm_sharedE > 0.7:
                    efficient_sim_tracksters.append(sim_idx)
                    # Count efficient trackster in bin
                    efficiency_bins[sim_bin] += 1

        # Get bin indices for RecoTracksters
        reco_energy_bins = bin_indices(tLinks['raw_energy'])

    # Fake Calculation
    for reco_idx in range(len(reco_tracksters)):
        # Ensure reco_bin is within valid range
        reco_bin = reco_energy_bins[reco_idx]
        if reco_bin < 0 or reco_bin >= num_bins:
            continue  # Skip invalid bin index
    
        # Count total RecoTracksters per bin
        reco_bin_counts[reco_bin] += 1
    
        # Get boolean array for fake & merge conditions
        fake_sim_list = reco_scores[reco_idx] < 0.6  # NumPy boolean array
        merge_sim_list = reco_scores[reco_idx] < 0.2  # Different threshold for merging?
    
        # Check fake tracksters (no valid RecoToSim association)
        if (np.sum(fake_sim_list) == 0 or np.sum(merge_sim_list) > 1):  # Count "True" values
#            print(f"Event {event} reco_idx {reco_idx} energy {tLinks['raw_energy'][reco_idx] }  {list(reco_scores[reco_idx][:10])}")
#            print(np.sum(fake_sim_list), np.sum(merge_sim_list))
            fake_rate_bins[reco_bin] += 1
            fake_reco_tracksters.append(reco_idx)
            continue  # Skip further checks for this reco_trackster
    
    # Compute efficiency and fake rate per bin
    efficiency_bins = np.divide(
        efficiency_bins, sim_bin_counts, where=sim_bin_counts > 0)
    fake_rate_bins = np.divide(
        fake_rate_bins, reco_bin_counts, where=reco_bin_counts > 0)
    #print(efficiency_bins)
    #print(fake_rate_bins)

    # Overall Efficiency and Fake Rate
    overall_efficiency = len(efficient_sim_tracksters) / \
        len(sim_to_reco_index) if len(sim_to_reco_index) > 0 else 0
    overall_fake_rate = len(fake_reco_tracksters) / \
        len(reco_to_sim_index) if len(reco_to_sim_index) > 0 else 0
    # Compute bin weights (higher weights for higher energy bins)
    # Midpoint of each bin
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    # Normalize so sum(weights) = 1
    weights = bin_centers / np.sum(bin_centers)

    # Compute weighted objectives
    weighted_efficiency = np.sum(weights * efficiency_bins)
    weighted_fake_rate = np.sum(weights * fake_rate_bins)

    # Adjust efficiency for minimization (1 - efficiency)
    # Now minimizing this improves efficiency
    objective_efficiency = 1 - weighted_efficiency
    objective_fake_rate = weighted_fake_rate  # Already minimizing fake rate
    print(f"Objective Efficiency {objective_efficiency}, Objective Fake-Merge rate {objective_fake_rate}")
    return objective_efficiency, objective_fake_rate


# read a csv file, return a matrix
def read_csv(filename):
    matrix = np.genfromtxt(filename, delimiter=",", dtype=float)
    if matrix.ndim == 2:
        return np.genfromtxt(filename, delimiter=",", dtype=float)
    return np.array([matrix])

# write a matrix to a csv file


def write_csv(filename, matrix):
    np.savetxt(filename, matrix, fmt='%.18f', delimiter=',')

#fileIn = "PSOTICLv5CLUE3D/simple_validation.root"
#fu = uproot.open(fileIn)
#print(get_metrics(fu, 0))

