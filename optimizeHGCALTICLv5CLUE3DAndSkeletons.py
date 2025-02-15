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
config = 'reconstructionTICLv5.py'

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
        print(f" Pop fitness {population_fitness}, {params}")
    return population_fitness

cylinder_radius_sqrEM_ub = 12.0 
cylinder_radius_sqrHAD_ub = 15.0
cylinder_radius_sqrEM_lb = 6.0 
cylinder_radius_sqrHAD_lb = 9.0
cylinder_radius_sqr_split_ub = 9.0
cylinder_radius_sqr_split_lb = 6.0 
deltaRxy_ub = 10.0
deltaRxy_lb = 2.0
dot_prod_th_ub = 0.98
dot_prod_th_lb = 0.90
lower_boundaryEM_ub = 30 
lower_boundaryHAD_ub = 30
lower_boundaryEM_lb = 10
lower_boundaryHAD_lb = 10
lower_distance_projective_sqrEM_ub = 40
lower_distance_projective_sqrHAD_ub = 40
lower_distance_projective_sqrEM_lb = 20
lower_distance_projective_sqrHAD_lb = 20
lower_distance_projective_sqr_closest_pointsEM_ub = 60
lower_distance_projective_sqr_closest_pointsHAD_ub = 60
lower_distance_projective_sqr_closest_pointsEM_lb = 20
lower_distance_projective_sqr_closest_pointsHAD_lb = 20
min_num_lcs_ub = 15
min_num_lcs_lb = 7 
min_trackster_energy_ub = 20
min_trackster_energy_lb = 5 
upper_boundaryEM_ub = 200
upper_boundaryHAD_ub = 200 
upper_boundaryEM_lb = 100 
upper_boundaryHAD_lb = 100
upper_distance_projective_sqrEM_ub = 40
upper_distance_projective_sqrHAD_ub = 70
upper_distance_projective_sqrEM_lb = 30
upper_distance_projective_sqrHAD_lb = 30

ub = [cylinder_radius_sqrEM_ub, cylinder_radius_sqrHAD_ub, cylinder_radius_sqr_split_ub, deltaRxy_ub, dot_prod_th_ub, lower_boundaryEM_ub, lower_boundaryHAD_ub, lower_distance_projective_sqrEM_ub, lower_distance_projective_sqrHAD_ub, lower_distance_projective_sqr_closest_pointsEM_ub, lower_distance_projective_sqr_closest_pointsHAD_ub, min_num_lcs_ub, min_trackster_energy_ub, upper_boundaryEM_ub, upper_boundaryHAD_ub, upper_distance_projective_sqrEM_ub, upper_distance_projective_sqrHAD_ub] 
lb = [cylinder_radius_sqrEM_lb, cylinder_radius_sqrHAD_lb, cylinder_radius_sqr_split_lb, deltaRxy_lb, dot_prod_th_lb, lower_boundaryEM_lb, lower_boundaryHAD_lb, lower_distance_projective_sqrEM_lb, lower_distance_projective_sqrHAD_lb, lower_distance_projective_sqr_closest_pointsEM_lb, lower_distance_projective_sqr_closest_pointsHAD_lb, min_num_lcs_lb, min_trackster_energy_lb, upper_boundaryEM_lb, upper_boundaryHAD_lb, upper_distance_projective_sqrEM_lb, upper_distance_projective_sqrHAD_lb ] 


for i in range(len(ub)):
    if(ub[i] <= lb[i]):
        print(ub[i], lb[i], i)
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
