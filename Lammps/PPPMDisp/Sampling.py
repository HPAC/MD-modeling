import sys, os
import subprocess
import math
import textwrap

sys.path.append(os.getcwd())
#sys.path.append(os.path.join(os.getcwd(), "pkg/subprocess32-3.2.6/build/lib.linux-x86_64-2.6/" ))
import subprocess32
import General.Config as Config
from General.utils import mult
import LammpsUtilities.PPPMUtilities as PPPMUtilities
from Lammps.PPPMDisp.Auxiliary.KSpaceGrids import get_grid_sizes
from Lammps.PPPMDisp.Auxiliary.InterfacePoints import get_interface_points
#from Lammps.PPPMDisp.Fitting import real_space_fitting_f, neigh_fitting_f, kspace_fitting_f
from Lammps.PPPMDisp.Auxiliary.ExtractTimings import extract_timings_from_log

real_space_sampling_params = {
    "ewald" : 0.5,
    "order" : 2,
    "grid"  : (1, 1, 1)
}
rsp = real_space_sampling_params

kspace_sampling_params = {
    "ewald"  : 0.5,
    #"cutoff" : 5.3
    "cutoff" : 2.0
}
ksp = kspace_sampling_params

d_sampled_cutoff = 0.1
def get_sampled_cutoffs( npart, domain, box ):
    test = Config.test
    cutoffs_at_interface = [ Config.cutoff_values[p[0]] \
            for p in get_interface_points( test.npart, test.domain, test.box ) ]
    cutoffs_at_interface = sorted(cutoffs_at_interface)
    c_min = cutoffs_at_interface[0]
    c_max = cutoffs_at_interface[-1]
    #c_min = 2.0
    #c_max = 6.1

    cutoffs = []
    c = c_min
    while c <= c_max:
        cutoffs.append(c)
        c += d_sampled_cutoff
    return cutoffs

def real_space_fitting_f( x, c, off ):
    return off + c*x**3

def neigh_fitting_f( x, c, off ):
    return off + c*x**3

def kspace_fitting_f( x, p, g ):
    return p + g*x

def sample_neigh( test ):
    print "** Sampling for Neighbor list"
    samples_dir = Config.get_test_samples_dir(Config.test)
    Config.create_dir( samples_dir )

    lmp_config = PPPMUtilities.LammpsPPPMConfig()
    lmp_config.npart = test.npart
    lmp_config.timesteps = test.timesteps
    lmp_config.ewald = real_space_sampling_params["ewald"]
    # cutoff set in the loop
    lmp_config.order = real_space_sampling_params["order"]
    lmp_config.grid = real_space_sampling_params["grid"]
    lmp_config.mixing = test.mixing
    lmp_config.diff = test.diff_mode
    #
    cutoffs = get_sampled_cutoffs( test.npart, test.domain, test.box )
    timings = [None for _ in cutoffs]

    ewald = rsp["ewald"]
    P     = rsp["order"]
    grid  = rsp["grid"]
    for i,c in enumerate(get_sampled_cutoffs(test.npart, test.domain, test.box)):
        try:
            t = extract_timings_from_log( os.path.join(samples_dir, "Output-%s-%s-%.2f-%.2f-%d-%dx%dx%d" % \
                    (test.diff_mode, test.mixing, ewald, c, P, grid[0], grid[1], grid[2])) )[3]
            timings[i] = t
        except Exception, e:
            #print e
            pass
    dynamic_sampling( lmp_config, "neigh", cutoffs, cutoffs, timings, neigh_fitting_f )

def sample_cutoff( test ):
    print "** Sampling for Real Space"
    samples_dir = Config.get_test_samples_dir(Config.test)
    Config.create_dir( samples_dir )

    lmp_config = PPPMUtilities.LammpsPPPMConfig()
    lmp_config.npart = test.npart
    lmp_config.timesteps = test.timesteps
    lmp_config.ewald = real_space_sampling_params["ewald"]
    # cutoff set in the loop
    lmp_config.order = real_space_sampling_params["order"]
    lmp_config.grid = real_space_sampling_params["grid"]
    lmp_config.mixing = test.mixing
    lmp_config.diff = test.diff_mode
    #
    cutoffs = get_sampled_cutoffs( test.npart, test.domain, test.box )
    timings = [None for _ in cutoffs]

    ewald = rsp["ewald"]
    P     = rsp["order"]
    grid  = rsp["grid"]
    for i,c in enumerate(get_sampled_cutoffs(test.npart, test.domain, test.box)):
        try:
            t = extract_timings_from_log( os.path.join(samples_dir, "Output-%s-%s-%.2f-%.2f-%d-%dx%dx%d" % \
                    (test.diff_mode, test.mixing, ewald, c, P, grid[0], grid[1], grid[2])) )[1]
            timings[i] = t
        except Exception, e:
            pass
    dynamic_sampling( lmp_config, "cutoff", cutoffs, cutoffs, timings, real_space_fitting_f )

def get_sampled_grid_sizes( npart, domain, box ):
    grids = get_grid_sizes( npart, domain, box )
    test = Config.test
    grids_at_interface = [ grids[p[2]] \
            for p in get_interface_points( test.npart, test.domain, test.box ) ]
    grids_at_interface = sorted(list(set(grids_at_interface)))
    return grids_at_interface

def get_sampled_Ps( npart, domain, box ):
    test = Config.test
    Ps_at_interface = [ Config.order_values[p[1]] \
            for p in get_interface_points( test.npart, test.domain, test.box ) ]
    Ps_at_interface = sorted(Ps_at_interface)
    return range(Ps_at_interface[0], Ps_at_interface[-1]+1)

def sample_kspace( test ):
    print "** Sampling for KSpace"
    samples_dir = Config.get_test_samples_dir(Config.test)
    Config.create_dir( samples_dir )

    lmp_config = PPPMUtilities.LammpsPPPMConfig()
    lmp_config.npart = test.npart
    lmp_config.timesteps = test.timesteps
    lmp_config.ewald  = kspace_sampling_params["ewald"]
    lmp_config.mixing = test.mixing
    lmp_config.diff = test.diff_mode
    #
    grids = get_sampled_grid_sizes( test.npart, test.domain, test.box )
    grids_cont = [mult(g) for g in grids]
    Ps = get_sampled_Ps( test.npart, test.domain, test.box )
    #
    cutoffs = get_sampled_cutoffs( test.npart, test.domain, test.box )
    c_min = cutoffs[0]
    c_max = cutoffs[-1]

    ewald  = kspace_sampling_params["ewald"]
    for c in [c_min, c_max]:
        lmp_config.cutoff = c
        for P in Ps:
            lmp_config.order = P
            timings = [None for _ in grids]
            for i, g in enumerate(grids):
                try:
                    t = extract_timings_from_log( os.path.join(samples_dir, "Output-%s-%s-%.2f-%.2f-%d-%dx%dx%d" % \
                            (test.diff_mode, test.mixing, ewald, c, P, g[0], g[1], g[2])) )[2]
                    if not math.isnan(t):
                        timings[i] = t
                except Exception, e:
                    print c, P, g, e
                    pass
            dynamic_sampling( lmp_config, "kspace", grids_cont, grids, timings, kspace_fitting_f )

import numpy as np
import scipy
from scipy.optimize import curve_fit

def fit( x, y, f ):
    popt, pcov = curve_fit(f, x, y)
    return popt

def error( xs, ys, f, *args ):
    max_error = 0
    for x,y in zip(xs,ys):
        y_est = f( x, *args )
        err = abs((y-y_est)/y)
        if err > max_error:
            max_error = err
    return max_error

def split_error( x, y, f, i ):
    args = fit( x[:i], y[:i], f )
    err_left = error( x[:i], y[:i], f, *args)
    args = fit( x[i:], y[i:], f )
    err_right = error( x[i:], y[i:], f, *args)
    return max(err_left, err_right)

def top_down( x, y, f, threshold = 0.05 ):
    def top_down_rec( x, y, f, threshold ):
        # if already good fit, stop
        vars = fit( x, y, f )
        if error( x, y, f, *vars ) < threshold or len(x) <= 3:
            return []

        best_so_far = 2**31 # ~ inf
        breakpoint = 0
        for i in range(2, len(x)-1):
            err = split_error( x, y, f, i )
            if err < best_so_far:
                best_so_far = err
                breakpoint = i
        return top_down_rec( x[:breakpoint], y[:breakpoint], f, threshold ) + \
               [breakpoint] + \
               [breakpoint + b for b in top_down_rec( x[breakpoint:], y[breakpoint:], f, threshold )]
    breakpoints = top_down_rec( x, y, f, threshold )

    # grid size bounds
    bounds = [0] + [x[i] for i in breakpoints] + [2**31]
    # fitting vars per range
    breakpoints = [0] + breakpoints + [len(x)]
    fit_vars = [None]
    for i in range(len(breakpoints)-1):
        fit_vars.append( fit(x[breakpoints[i] : breakpoints[i+1]],
                             y[breakpoints[i] : breakpoints[i+1]], f) )
    return (bounds, fit_vars)

def sample( lmp_config, target, value ):
    if target in ("cutoff", "neigh"):
        lmp_config.cutoff = value
        if target == "cutoff":
            log_pos = 1
        else:
            log_pos = 3
        print "Sampling cutoff:", value
    elif target == "kspace": # gonna need to deal with grid and mult(grid)
        lmp_config.grid = value
        log_pos = 2
        print "Sampling grid:", value
    #return
    # Set paths and write lammps config file
    lmp_config.set_file_paths( Config.get_test_samples_dir(Config.test) )
    lmp_config.write_config_file( lmp_config.infile )
    #
    done = False
    attempt = 1
    while not done and attempt <= 10:
        try:
            os.remove( lmp_config.outfile )
        except OSError:
            pass
        try:
            #subprocess32.check_call("mpirun -np %d $HOME/MD-libs/lammps-22Jan14/src/lmp_openmpi -l /dev/null < %s > %s" % (Config.test.nprocesses, lmp_config.infile, lmp_config.outfile), shell=True, timeout=1500)
            subprocess32.check_call("$MPIEXEC -np %d $HOME/MD-libs/lammps-22Jan14/src/lmp_openmpi -l /dev/null < %s > %s" % (Config.test.nprocesses, lmp_config.infile, lmp_config.outfile), shell=True, timeout=1500)
            t = extract_timings_from_log( lmp_config.outfile )[log_pos]
            if not math.isnan(t):
                done = True
        except subprocess32.TimeoutExpired:
            pass
        attempt += 1
    if done:
        return t
    else:
        raise Exception

# x_cont: g[0] * g[1] * g[2]
# x_disc: (g[0], g[1], g[2])
def dynamic_sampling( lmp_config, target, x_cont, x_disc, y, f, threshold = 0.05 ):
    def dyn_samp_rec( i, j, f, threshold ):
        if (j-i) == 1:
            return 
        # run and parse mid point
        midpoint = (i+j)/2
        if not y[midpoint]:
            y[midpoint] = sample( lmp_config, target, x_disc[midpoint] )
        _x = np.array([x_cont[i], x_cont[midpoint], x_cont[j]])
        _y = np.array([y[i], y[midpoint], y[j]])

        print error(_x, _y, f, *fit(_x, _y, f))
        if error(_x, _y, f, *fit(_x, _y, f)) > threshold:
            print "recursive", i, midpoint, j
            # left sampling
            left = dyn_samp_rec( i , midpoint, f, threshold )
            # right sampling
            right = dyn_samp_rec( midpoint, j, f, threshold )

    if not y[0]:
        y[0]  = sample( lmp_config, target, x_disc[0] )
    if not y[-1]:
        y[-1] = sample( lmp_config, target, x_disc[-1] )
    dyn_samp_rec( 0, len(x_cont)-1, f, threshold )

def run_samples( test, diff_mode ):
    Config.test = test
    Config.test.diff_mode = diff_mode
    #
    sample_cutoff( test )
    sample_neigh( test )
    sample_kspace( test )

if __name__ == "__main__":
    import sys
    import TestSuite

    if len(sys.argv) != 3:
        print >> sys.stderr, "Usage: %s TestName diff" % sys.argv[0]
        sys.exit(-1)

    test = TestSuite.TestSuite[sys.argv[1]]
    Config.test = test
    Config.test.diff_mode = sys.argv[2]

    sample_cutoff( test )
    sample_neigh( test )
    sample_kspace( test )
