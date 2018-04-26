#!/usr/bin/env python

import os

#
# Global var to simplify code transition
#   should cleanup later on
test = None

def get_test_root_dir( test ):
	return os.path.join(results_dir, test.name)
def get_test_p_dir( test ):
	return os.path.join( get_test_root_dir( test ), "%dp" % test.nprocesses)
def get_test_acc_dir( test ):
    if test.accuracy:
    	return os.path.join( get_test_p_dir( test ), "%.3e" % test.accuracy)
    else:
    	return os.path.join( get_test_p_dir( test ), "%.3e" % test.accuracy_kspace)
def get_test_method_dir( test ):
	return os.path.join( get_test_acc_dir( test ), test.diff_mode )
def get_test_samples_dir( test ):
	return os.path.join( get_test_method_dir( test ), "Samples" )
def get_test_timings_dir( test ):
	return os.path.join( get_test_method_dir( test ), "Timings" )

def get_test_kspace_error_grid_config_file( test ):
	return os.path.join( get_test_method_dir( test ), "kspace-error.cfg" )
def get_test_kspace_error_grid_file( test ):
	return os.path.join( get_test_method_dir( test ), "kspace-error.txt" )

# Parameter choices
#
# General
#
# short range cutoff
d_cutoff = 0.1
min_cutoff = 2.0
max_cutoff = 6.0
cutoff_values = [min_cutoff + i*d_cutoff for i in range( int( (max_cutoff - min_cutoff)/d_cutoff + 1 ) ) ]

#
# Algorithm specific
#
# PPPM
d_ewald = 0.01
min_ewald = d_ewald
max_ewald = 1.0
ewald_values = [min_ewald + i*d_ewald for i in range( int( (max_ewald - min_ewald)/d_ewald+ 1 ) ) ]
order_values = [2, 3, 4, 5, 6]
mixing_values = ["geom", "arit", "none"]
diff_values = ["ad", "ik"]

# Directories
modeling_dir = "."

results_dir = os.path.join(modeling_dir, "Results")
lammps_config_dir = os.path.join(modeling_dir, "config/")
jobs_dir = os.path.join(modeling_dir, "jobs/")
rolf_flat_output_dir = os.path.join(modeling_dir, "Results/Rolf-flat/")

# Directory creation and walking
def create_dir( dirname ):
	if os.path.exists(dirname):
		if not os.path.isdir(dirname):
			raise Exception( "% exists and is not a directory. Exiting" % dirname )
	else:
		#os.mkdir( dirname )
		os.makedirs( dirname )

def create_dirs( dirs_list ):
	if dirs_list:
		for d in dirs_list[0]:
			cwd = os.getcwd()
			create_dir( d )
			os.chdir( d )
			create_dirs( dirs_list[1:] )
			os.chdir( cwd )

def walk_dirs( dirs_list, f, data=None ):
	if dirs_list:
		for d in dirs_list[0]:
			cwd = os.getcwd()
			os.chdir( d )
			walk_dirs( dirs_list[1:], f, data )
			os.chdir( cwd )
	else:
		if data is not None:
			f(data)
		else:
			f()
