import optimizer
import subprocess
import math as m
from utils import get_metrics, write_csv
import numpy as np
import uproot
import argparse
import os

# parsing argument
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--continuing', type=int, action='store')
parser.add_argument('-d', '--default', action='store_true')
parser.add_argument('-p2', '--phase2', action='store_true')
parser.add_argument('-p', '--num_particles', default=5,
                    type=int, action='store')
parser.add_argument('-i', '--num_iterations',
                    default=1, type=int, action='store')
parser.add_argument('-e', '--num_events', default=100,
                    type=int, action='store')
args = parser.parse_args()

num_iterations = args.num_iterations

##EM PARAMS
optimizer.Logger.setLevel('DEBUG')

optimizer.Randomizer.rng = np.random.default_rng(46)
#Critical Density
defaults = []


for p in defaults:
    print(f"{p:.18f}", end=',')
config = 'reconstructionHGCALTICLv5PU75KSkeletonsOnly.py'
input_file = 'step3.root'

# configure parameters
critical_density_lb = [0.4,0.4]
critical_density_ub = [1.2,1.2]
critical_self_density_lb = [0.01, 0.01]
critical_self_density_ub = [0.3, 0.3]
critical_xy_distance_lb = [1.]
critical_xy_distance_ub = [2.5]
critical_z_distance_lb = [4, 4]
critical_z_distance_ub = [7, 7]
density_sibling_lay_lb = [2,2]
density_sibling_lay_ub = [5,5]
density_xydist_lay_lb = [2.0,2.0]
density_xydist_lay_ub = [4.0,4.0]
kernel_density_fac_lb = [0.1, 0.1]
kernel_density_fac_ub = [0.5, 0.5]
outlier_mul_lb = [1., 1.]
outlier_mul_ub = [2., 2.]

lb = []
lb.extend(critical_density_lb)
lb.extend(critical_self_density_lb)
lb.extend(critical_xy_distance_lb)
lb.extend(critical_z_distance_lb)
lb.extend(density_sibling_lay_lb)
lb.extend(density_xydist_lay_lb)
lb.extend(kernel_density_fac_lb)
lb.extend(outlier_mul_lb)

ub = []
ub.extend(critical_density_ub)
ub.extend(critical_self_density_ub)
ub.extend(critical_xy_distance_ub)
ub.extend(critical_z_distance_ub)
ub.extend(density_sibling_lay_ub)
ub.extend(density_xydist_lay_ub)
ub.extend(kernel_density_fac_ub)
ub.extend(outlier_mul_ub)

working_dir = 'PSOTICLv5CLUE3D'
def reco_and_validate(params):
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)
    write_csv(f'{working_dir}/parameters.csv', params)
    validation_result = f'{working_dir}/simple_validation.root'
    subprocess.run(['cmsRun', config, 'nEvents=' + str(args.num_events),
         f'parametersFile={working_dir}/parameters.csv', 'outputFile=' + validation_result]
                     )
    print('cmsRun', config, 'nEvents=' + str(args.num_events),f'parametersFile={working_dir}/parameters.csv', 'outputFile=' + validation_result)
    num_particles = len(params)
    with uproot.open(validation_result) as uproot_file:
        #print(f"Get Metric {get_metrics(uproot_file,0)}")
        population_fitness = np.array(
            [get_metrics(uproot_file, i) for i in range(num_particles)], dtype = float)
#    print(f" Pop fitness {population_fitness}, {params}")
    return population_fitness


# get default metrics
if args.default:
    defaults = [0.6,0.6,0.15,0.15,1.8,5,5,3,3,3.24,3.24,0.2,0.2,2.0,2.0]
    print(f' Len defaults {len(defaults)}')
    default_params = [defaults]
    default_metrics = reco_and_validate(default_params)
    write_csv(f'{working_dir}/default.csv',
              [np.concatenate([default_params[0], default_metrics[0]])])

objective = optimizer.Objective(reco_and_validate, 2)
optimizer.FileManager.working_dir=working_dir
optimizer.FileManager.loading_enabled = False 
optimizer.FileManager.saving_enabled = True

print(f"Len ub {len(ub)}, lb {len(lb)}")
pso = optimizer.MOPSO(objective=objective, lower_bounds=lb, upper_bounds=ub, 
            num_particles=args.num_particles,
            inertia_weight=0.5, cognitive_coefficient=1.5, social_coefficient=1.5)

pso.optimize(num_iterations, max_iterations_without_improvement=10)
