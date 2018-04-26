#!/usr/bin/env python

import argparse
import sys, os, re
import math

import numpy

from General.utils import mult
import General.Config as Config

from Lammps.PPPMDisp.Fitting import fit_samples as pppm_fit_samples
from Lammps.PPPMDisp.Empirical import get_interface_timings
from Lammps.PPPMDisp.Auxiliary.KSpaceGrids import get_grid_sizes
from Lammps.PPPMDisp.Auxiliary.InterfacePoints import get_interface_points

def get_whole_interface_prediction( npart, domain, box, mixing, np, nn, ts, diff, acc ):
	fitting_vars = pppm_fit_samples( npart, domain, box, mixing, np, nn, ts, diff )
	
	grids = get_grid_sizes( npart, domain, box )
	# c estimator
	f_real_space_prediction = fitting_vars[0]
	f_neigh_prediction = fitting_vars[1]
	f_kspace_prediction = fitting_vars[2]
	estimated_timings = []
	for params in get_interface_points( npart, domain, box ):
		(cutoff_idx, P_idx, grid_idx, ewald_error) = params
		## parameters of interest
		cutoff = Config.cutoff_values[cutoff_idx]
		grid = grids[grid_idx]
		n_grid_points = mult(grid)
		## estimated time
		pt = f_real_space_prediction(cutoff)
		nt = f_neigh_prediction(cutoff)
		dt = f_kspace_prediction(cutoff, P_idx, n_grid_points)
		#tt = pt + dt # c*cutoff**3 + p + g*n_grid_points
		tt = pt + nt + dt
		# store prediction
		estimated_timings.append(((tt, pt, dt, nt), params))
	return estimated_timings

def get_min_prediction_from_interface( interface_prediction, npart, domain, box ):
	grids = get_grid_sizes( npart, domain, box )
	# find minimum timing
	min_est_t, min_est_params = min(interface_prediction)
	min_est_ci, min_est_pi, min_est_gi, ewald_error = min_est_params
	# extract parameters of fastest configuration
	cutoff = Config.cutoff_values[min_est_ci]
	P      = Config.order_values[min_est_pi]
	grid   = grids[min_est_gi]

	return (cutoff, P, grid, ewald_error[0][0], min_est_t)

def get_prediction( npart, domain, box, np, nn, ts, acc ):
	estimated_timings = get_whole_interface_prediction( npart, domain, box, np, nn, ts, acc )
	return get_min_prediction_from_interface( estimated_timings, npart, domain, box )

def get_prediction_with_statistics( npart, domain, box, mixing, np, nn, ts, diff, acc ):
	estimated_timings = get_whole_interface_prediction( npart, domain, box, mixing, np, nn, ts, diff, acc )
	empirical_timings = get_interface_timings( npart, domain, box, mixing, np, nn, ts, diff, acc )

	grids = get_grid_sizes( npart, domain, box )
	differences = []
	for i, (est, emp) in enumerate(zip( estimated_timings, empirical_timings )):
		t_est = est[0][0]
		#
		cutoff = Config.cutoff_values[est[1][0]]
		P = Config.order_values[est[1][1]]
		grid = grids[est[1][2]]
		print "%2d | %6.3f %6.3f %6.3f %6.3f | %4.2f %d (%2d, %2d %3d)" % (i, t_est, est[0][1], est[0][3], est[0][2], cutoff, P, grid[0], grid[1], grid[2] )
		t_emp = emp[0][0]
		#t_emp = emp[0][1] + emp[0][2]
		print "   | %6.3f %6.3f %6.3f %6.3f | %4.2f %d (%2d, %2d %3d)" % (t_emp, emp[0][1], emp[0][3], emp[0][2], cutoff, P, grid[0], grid[1], grid[2] )
		#print "   | %6.3f %6.3f %6.3f | %4.2f %d (%2d, %2d %3d)" % (emp[0][0], emp[0][1], emp[0][2], cutoff, P, grid[0], grid[1], grid[2] )
		print
		differences.append( abs(t_emp - t_est) )
	
	min_diff = min(differences)
	max_diff = max(differences)
	mean_diff = numpy.mean(differences)
	median_diff = numpy.median(differences)
	std_diff = numpy.std(differences)

	min_estimated = get_min_prediction_from_interface( estimated_timings, npart, domain, box )
	return ( min_estimated, (min_diff, max_diff, mean_diff, median_diff, std_diff), (estimated_timings, empirical_timings) )
