import sys
import os
import math

import numpy
import scipy
from scipy.optimize import curve_fit

sys.path.append(os.getcwd())
import General.Config as Config
from General.utils import mult
import Lammps.PPPMDisp.Sampling as Sampling
from Sampling import real_space_sampling_params as rsp, kspace_sampling_params as ksp
from Lammps.PPPMDisp.Sampling import get_sampled_cutoffs, get_sampled_grid_sizes, get_sampled_Ps
from Lammps.PPPMDisp.Auxiliary.ExtractTimings import extract_timings_from_log

import General.topdown as topdown

_real_piecewise = []
_neigh_piecewise = []
_kspace_piecewise = [[],[]]

_c = []

def get_kspace_prediction( cutoff, P_idx, grid_size ):
    if _kspace_piecewise == []:
        raise Exception # improve

    bounds, fit = _kspace_piecewise[0][P_idx]
    i = 0
    while grid_size >= bounds[i]:
        i += 1
    left = kspace_fitting_f( grid_size, *fit[i] )

    bounds, fit = _kspace_piecewise[1][P_idx]
    i = 0
    while grid_size >= bounds[i]:
        i += 1
    right = kspace_fitting_f( grid_size, *fit[i] )

    popt = topdown.fit( numpy.array([_c[0], _c[1]]), numpy.array([left, right]), across_k_fitting_f )
    return across_k_fitting_f( cutoff, *popt )

def get_real_space_prediction( cutoff ):
    if _real_piecewise == []:
        raise Exception # improve
    #
    bounds, c = _real_piecewise
    i = 0
    while cutoff > (bounds[i] - 1e-4):
        i += 1
    return real_space_fitting_f( cutoff, *c[i] )

def get_neigh_prediction( cutoff ):
    if _neigh_piecewise == []:
        raise Exception # improve
    #
    bounds, fit = _neigh_piecewise
    i = 0
    while cutoff > (bounds[i] - 1e-4):
        i += 1
    return neigh_fitting_f( cutoff, *fit[i] )

def fit_cutoff( test ):
    samples_dir = Config.get_test_samples_dir(Config.test)
    cutoffs, timings = [], []

    ewald = rsp["ewald"]
    P     = rsp["order"]
    grid  = rsp["grid"]
    for c in get_sampled_cutoffs(test.npart, test.domain, test.box):
        try:
            t = extract_timings_from_log( os.path.join(samples_dir, "Output-%s-%s-%.2f-%.2f-%d-%dx%dx%d" % \
                    (test.diff_mode, test.mixing, ewald, c, P, grid[0], grid[1], grid[2])) )[1]
            cutoffs.append(c)
            timings.append(t)
        except Exception, e:
            pass
    bounds, fit_vars = topdown.top_down( numpy.array(cutoffs),
                                         numpy.array(timings),
                                         real_space_fitting_f )
    _real_piecewise.extend((bounds, fit_vars))

def fit_neigh( test ):
    samples_dir = Config.get_test_samples_dir(Config.test)
    cutoffs, timings = [], []

    ewald = rsp["ewald"]
    P     = rsp["order"]
    grid  = rsp["grid"]
    for c in get_sampled_cutoffs(test.npart, test.domain, test.box):
        try:
            t = extract_timings_from_log( os.path.join(samples_dir, "Output-%s-%s-%.2f-%.2f-%d-%dx%dx%d" % \
                    (test.diff_mode, test.mixing, ewald, c, P, grid[0], grid[1], grid[2])) )[3]
            cutoffs.append(c)
            timings.append(t)
        except Exception, e:
            pass
    bounds, fit_vars = topdown.top_down( numpy.array(cutoffs),
                                         numpy.array(timings),
                                         neigh_fitting_f )
    _neigh_piecewise.extend((bounds, fit_vars))


def fit_samples( npart, domain, box, mixing, np, nn, ts, diff ):
    cwd = os.getcwd()
    os.chdir(Config.get_test_samples_dir(Config.test))

    fits = []

    fit_cutoff( Config.test )
    fits.append(get_real_space_prediction)

    fit_neigh( Config.test )
    fits.append(get_neigh_prediction)

    # kspace samples 
    ewald  = ksp["ewald"]
    cutoffs = get_sampled_cutoffs( npart, domain, box )
    _c.extend( (cutoffs[0], cutoffs[-1]) )
    Ps = get_sampled_Ps( npart, domain, box )
    grids = get_sampled_grid_sizes( npart, domain, box )
    for i,cutoff in enumerate([cutoffs[0], cutoffs[-1]]):
        _kspace_piecewise[i].extend([None for _ in range(Ps[0] - Config.order_values[0])])
        for P in Ps:
            timings = []
            gs = []
            for grid in grids:
                try:
                    t = extract_timings_from_log( "Output-%s-%s-%.2f-%.2f-%d-%dx%dx%d" % (diff, mixing, ewald, cutoff, P, grid[0], grid[1], grid[2]) )[2]
                    gs.append(grid)
                    timings.append(t)
                except Exception, e:
                    pass
            try:
                x, y = filter_data( [mult(g_) for g_ in gs], timings )
                bounds, fit_vars = topdown.top_down( numpy.array(x), numpy.array(y), kspace_fitting_f )
                _kspace_piecewise[i].append((bounds, fit_vars))
            except ValueError, e:
                print e
                print >> sys.stderr, "Cannot fit. Infs or Nans found."
                (p, g) = (float('NaN'), float('NaN'))
                raise Exception
    fits.append(get_kspace_prediction)

    os.chdir(cwd)

    return fits

def real_space_fitting_f( x, c, off ):
    return off + c*x**3

def neigh_fitting_f( x, c, off ):
    return off + c*x**3

def kspace_fitting_f( x, p, g ):
    return p + g*x

def across_k_fitting_f( x, c, off ):
    return off + c*x**3

# For now, if two equal entries in x,
#  use only the first.
# Assumes x entries are sorted
def filter_data( x, y ):
    _x, _y = [], []
    for i in range( 1, len( x ) ):
        if x[i] != x[i-1]:
            _x.append( x[i] )
            _y.append( y[i] )
    return _x, _y

if __name__ == "__main__":
	import sys
	import TestSuite

	if len(sys.argv) != 2:
		print >> sys.stderr, "Usage: %s TestName" % sys.argv[0]
		sys.exit(-1)

	test = TestSuite.TestSuite[sys.argv[1]]
	Config.test = test
	Config.test.diff_mode = "ad"

	fit_cutoff( test )
