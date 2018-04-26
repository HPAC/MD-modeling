import math
import os.path
import subprocess

import General.Config as Config
import General.utils as utils
from Lammps.PPPMDisp.Auxiliary.KSpaceGrids import get_grid_sizes
from Lammps.PPPMDisp.Auxiliary.KSpaceErrorGridConfig import generate_kspace_error_config_file

#
# Real space error estimates for a single set of parameters
#
def error_real( V, N, C, cutoff, beta ):
    return (C * math.sqrt(math.pi) * beta**5) / (math.sqrt(N * V * cutoff)) * \
            ((6. / (cutoff*beta)**6) + (6. / (cutoff*beta)**4) + \
             (3. / (cutoff*beta)**2) + 1) * \
            math.exp(-(cutoff*beta)**2)

#
# Grids of error estimates
#
def get_error_real_grid( box, N, C ): # domain?
    V = reduce( lambda x,y: x*y, box )
    #V = utils.mult( box )
    grid = []
    for cutoff in Config.cutoff_values:
        grid.append([])
        for ewald in Config.ewald_values:
            grid[-1].append( error_real( V, N, C, cutoff, ewald ) )
    return grid

# Externalized to C + MPI code for performance reasons
def get_error_reciprocal_grid( domain, box, N, C ):
    error_grid_file = Config.get_test_kspace_error_grid_file( Config.test )
    if not os.path.isfile( error_grid_file ):
        print "Error estimates not available, running..."
        config_file = Config.get_test_kspace_error_grid_config_file( Config.test )
        f = open(config_file, "w")
        generate_kspace_error_config_file( Config.test, f )
        f.close()
        subprocess.check_call("mpirun -np %d C/MPI-qopt %s %s > %s" % \
                (Config.test.nprocesses, Config.test.diff_mode, config_file, error_grid_file), shell=True)
        print " Done"
        #raise Exception # abort until done

    grids = get_grid_sizes( N, domain, box )
    data = []
    with open( error_grid_file, "r" ) as fp:
        for P in Config.order_values:
            data.append([])
            for grid in grids:
                data[-1].append([])
                for ewald in Config.ewald_values:
                    line = fp.readline()
                    err = float(line.strip().split()[-1])
                    data[-1][-1].append( err )
    return data

# 
# Mask: [] -> combination of cutoff, P, and grid, results
#             in inaccurate results (according to the error bounds)
#       else -> the estimated accuracy
#
def get_accurate_grid( N, domain, box, _ ):
    c = [2.0 for i in range(N)]
    C = sum([ch**2 for ch in c])
    
    # Either accuracy is set, or the other two are
    accuracy   = Config.test.accuracy
    acc_real   = Config.test.accuracy_real
    acc_kspace = Config.test.accuracy_kspace
    #if not accuracy:
        #acc_real = accuracy
    if not acc_real:
        acc_real = accuracy
    if not acc_kspace:
        acc_kspace = accuracy

    error_real_grid = get_error_real_grid( box, N, C )
    error_reciprocal_grid = get_error_reciprocal_grid( domain, box, N, C )

    # Get min and max possible beta, given the ranges in Config
    # and the desired accuracy
    beta_idx = 0
    error = error_real_grid[-1][beta_idx]
    while error > acc_real and beta_idx < len(Config.ewald_values)-1: #(len(error_real_grid[3][-1])-1):
        beta_idx += 1
        error = error_real_grid[-1][beta_idx]
    if error > acc_real: # increase max_beta and recursively call?
        raise Exception # improve
    min_beta_idx = beta_idx
    min_beta = Config.ewald_values[0] + beta_idx * Config.d_ewald

    max_beta_idx = -1
    for P_idx in range(len(Config.order_values)):
        beta_idx = len(Config.ewald_values)-1 #len(error_reciprocal_grid[0])-1
        error = error_reciprocal_grid[P_idx][-1][beta_idx]
        while error > acc_kspace and beta_idx > 0:
            beta_idx -= 1
            error = error_reciprocal_grid[P_idx][-1][beta_idx]
        if not error > acc_kspace and beta_idx > max_beta_idx:
            max_beta_idx = beta_idx
    if max_beta_idx == -1:
        raise Exception # improve
    max_beta = Config.ewald_values[0] + max_beta_idx * Config.d_ewald

    # Calculate the error estimate within the accurate region
    # The rest is set to 0
    grids = get_grid_sizes( N, domain, box )
    error_grid = [
                   [ 
                     [[] for i in range(len(grids))] for j in range(len(Config.order_values))
                   ] for k in range(len(Config.cutoff_values)) 
                 ]
    # Could iterate over regions with a relaxed accuracy constraint
    # and then check actual accuracy before counting it in for min and max flops
    for beta_idx in range(min_beta_idx, max_beta_idx+1):
        for cutoff_idx, cutoff in enumerate(Config.cutoff_values):
            error_real = error_real_grid[cutoff_idx][beta_idx]
            if error_real > acc_real: # could run backwards and cut as soon as "inaccurate"
                continue
            for P_idx, P in enumerate(Config.order_values):
                for grid_idx, grid in enumerate(grids):
                    error_reciprocal = error_reciprocal_grid[P_idx][grid_idx][beta_idx]
                    if accuracy: # single overall, then norm
                        error = math.sqrt(error_real**2 + error_reciprocal**2)
                        accurate = error < accuracy
                    else: # Otherwise, satisfy error is less than acc_ksapce
                        error = error_reciprocal
                        accurate = error < acc_kspace
                    if accurate:
                        error_grid[cutoff_idx][P_idx][grid_idx].append( (
                                  Config.ewald_values[0] + Config.d_ewald * beta_idx,
                                  error
                                ) )
    return error_grid

if __name__ == "__main__":
    # Isele-Holder et al. J. Chem. Phys. 137, 174107 (2012)
    # Section 3.C (Figs. 2 and 3)

    boxsize = [15.0, 15.0, 15.0]
    npart = 2000

    V = reduce( lambda x,y: x*y, boxsize )
    N = npart
    c = [2.0 for i in range(N)]
    C = sum([ch**2 for ch in c])

    beta = 0.8
    cutoff = 2.0
    print error_real( V, N, C, cutoff, beta )
    cutoff = 3.0
    print error_real( V, N, C, cutoff, beta )
    cutoff = 4.0
    print error_real( V, N, C, cutoff, beta )

    iorder = 5
    h = boxsize[0]/32.
    print error_reciprocal( V, N, C, beta, h, iorder )
