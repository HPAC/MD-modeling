#!/usr/bin/env python

import math
import itertools

#import os, sys
#sys.path.append(os.getcwd())
from General.utils import mult

def get_range_grid_side( N_side, size_side ):
	grid_min = 2
	grid_max = int(math.ceil( N_side ))
	grid_max *= 2

	min_pow5 = 0
	max_pow5 = int(math.floor(math.log(grid_max, 5)))
	min_pow3 = 0
	max_pow3 = int(math.floor(math.log(grid_max, 3)))
	min_pow2 = 0
	max_pow2 = int(math.floor(math.log(grid_max, 2)))
	grid_sides = []
	for pow5 in range(min_pow5, max_pow5+1):
		for pow3 in range(min_pow3, max_pow3+1):
			for pow2 in range(min_pow2, max_pow2+1):
				grid_sides.append(5**pow5 * 3**pow3 * 2**pow2)

	grid_sides = sorted(set(grid_sides))
	for idx, grid_side in enumerate(grid_sides):
		if grid_side > grid_max:
			idx_beyond_grid_max = idx
			break
	return grid_sides[1:idx_beyond_grid_max] # First value (1) is not valid

def get_grid_sizes( N, domain, box=None ):
	# Assume full domain with same density as in box (if box).
	# Then, assuming regular distribution, take the number of particles
	#   along the largest side and use it to get the range of grid sides
	if box:
		N_extrapolated = N * mult([ d/b for d,b in zip(domain,box) ])
	else:
		N_extrapolated = N
	sorted_sides = sorted(domain)
	largest_side_particles = N_extrapolated / sorted_sides[0] / sorted_sides[1]
	gr = get_range_grid_side( largest_side_particles, sorted_sides[2] )
	# Allow min grid_size (for largest side) to be half of size_side
	min_grid_side = 0
	while gr[min_grid_side] < ( sorted_sides[2] / 2.0): # max sep between grid points 2*sigma
		min_grid_side += 1
	# We want to keep grid spacing even. Thus, for each value of grid side
	# of one dimension, the others have to get values for a similar grid spacing
	max_side = float(max(domain)) # should be float already
	factors = [ side/max_side for side in domain ]
	grids = []
	for size in gr[min_grid_side:]:
		grid = [ f * size for f in factors ]
		# adjust the grid to valid values
		adj_grids = adjust_grid( grid, gr )
		grids.extend(adj_grids)
	#return grids
	# Sort by decreasing spacing, then increasing number of grid points 
	# (from less accurate to more accurate, then from fastest to slowest)
	return sorted(grids, key=lambda grid: (min( grid_side/domain_side for domain_side, grid_side in zip (domain, grid) ), mult(grid)))

# set the grid sides that do not match any in grid_range
# to the one before and the one after
def adjust_grid( grid, grid_range ):
	possible_values = []
	for side in grid:
		side = int(side)
		if side in grid_range:
			possible_values.append([side])
		else:
			i = 0
			while grid_range[i] < side: # should never go beyond the length of the range
				i += 1
			if i == 0: # possible?
				possible_values.append([grid_range[i]])
			else:
				possible_values.append([grid_range[i-1], grid_range[i]])
	return list(itertools.product(*possible_values))

if __name__ == "__main__":
	N = 6000
	domain = (11.01, 11.01, 66.06)
	print get_grid_sizes( N, domain )
	N = 4000
	domain = (11.01, 11.01, 176.16)
	box = (11.01, 11.01, 44.04)
	print get_grid_sizes( N, domain, box )
