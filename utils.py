import numpy as np

# calculate the metrics from validation results
def get_metrics(uproot_file, id):
    tree = uproot_file['simpleValidation' + str(id)]['output']
#    print(tree.keys())
    total_recC3D = tree['total_recoCLUE3D'].array()[0]
    total_recMerged = tree['total_recoMerged'].array()[0]
    number_of_fake = tree['number_of_fake'].array()[0]
    number_of_merge = tree['number_of_merge'].array()[0]
    number_of_sim = tree['number_of_sim'].array()[0]
    number_of_sim_eff = tree['number_of_sim_eff'].array()[0]
    number_of_eff = tree['number_of_eff'].array()[0]
 

    return [total_recC3D/number_of_sim, number_of_fake / total_recMerged, 1 - number_of_eff/number_of_sim_eff] 

# read a csv file, return a matrix
def read_csv(filename):
    matrix = np.genfromtxt(filename, delimiter=",", dtype=float)
    if matrix.ndim == 2:
        return np.genfromtxt(filename, delimiter=",", dtype=float)
    return np.array([matrix])
    
# write a matrix to a csv file
def write_csv(filename, matrix):
    np.savetxt(filename, matrix, fmt='%.18f', delimiter=',')
