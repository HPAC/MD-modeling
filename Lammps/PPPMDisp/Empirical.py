import sys, os
import subprocess

import textwrap

sys.path.append(os.getcwd())
import General.Config as Config

from Lammps.PPPMDisp.Auxiliary.KSpaceGrids import get_grid_sizes
from Lammps.PPPMDisp.Auxiliary.InterfacePoints import get_interface_points
from Lammps.PPPMDisp.Auxiliary.ExtractTimings import extract_timings_from_log

import LammpsUtilities.PPPMUtilities as PPPMUtilities
from JobSubmission.LSF import JobScript

def run_interface_timings( npart, domain, box, mixing, np, nn, ts, diffmode, accuracy ):
	grids = get_grid_sizes( npart, domain, box )
	nprocs = np

	#Config.create_dir( timings_dir )
	timings_dir = Config.get_test_timings_dir(Config.test)
	Config.create_dir( timings_dir )

	# LSF Job submission
	lsf_config = JobScript.LSF_Config()
	lsf_config.group = "aices"
	lsf_config.jobname = "PPPM-Timings"
	#lsf_config.time = "%d:00" % (min(120, 2.5*np))
	lsf_config.time = "2:00"
	lsf_config.memory = "1000"
	lsf_config.nthreads = np
	lsf_config.parallelism_type = "openmpi"
	lsf_config.arch_string = "Harpertown"
	# Command
	lsf_config.command = []
	# Lammps Config
	lmp_config = PPPMUtilities.LammpsPPPMConfig()
	lmp_config.npart = npart
	lmp_config.timesteps = ts
	lmp_config.mixing = mixing
	lmp_config.diff = diffmode
	for (cutoff_idx, P_idx, grid_idx, ewald_error) in get_interface_points( npart, domain, box ):
		lmp_config.ewald = ewald_error[0][0]
		lmp_config.cutoff = Config.cutoff_values[cutoff_idx]
		lmp_config.order = Config.order_values[P_idx]
		lmp_config.grid = grids[grid_idx]

		print lmp_config.ewald
		print lmp_config.cutoff, "(", cutoff_idx, ")" 
		print lmp_config.order, "(", P_idx,  ")" 
		print lmp_config.grid, "(", grid_idx,  ")" 

		lmp_infile = os.path.join( 
			timings_dir, 
			"Lammps-%s-%s-%.2f-%.2f-%d-%dx%dx%d" % (
				lmp_config.diff, 
				lmp_config.mixing, 
				lmp_config.ewald, 
				lmp_config.cutoff, 
				lmp_config.order, 
				lmp_config.grid[0], lmp_config.grid[1], lmp_config.grid[2]
			)
		)
		lmp_config.write_config_file( lmp_infile )
		# LSF
		outputfile = os.path.join( 
			timings_dir, 
			"Output-%s-%s-%.2f-%.2f-%d-%dx%dx%d" % (
				lmp_config.diff, 
				lmp_config.mixing, 
				lmp_config.ewald, 
				lmp_config.cutoff, 
				lmp_config.order, 
				lmp_config.grid[0], lmp_config.grid[1], lmp_config.grid[2]
			)
		)
		lsf_config.command.append(textwrap.dedent("""
			export OMP_NUM_THREADS=1
			mpirun -np %d $HOME/MD-libs/lammps-22Jan14/src/lmp_openmpi -l /dev/null < %s > %s
		""" % (np, lmp_infile, outputfile)))

	# Script
	lsf_config.command = "\n".join(lsf_config.command)
	jobscript = os.path.join( 
		timings_dir, 
		"Jobscript-Interface"
	)
	lsf_config.write_jobscript( jobscript )

	print "Submitting job (Interface timing)... ",
	subprocess.check_call("bsub < %s" % jobscript, shell=True)
	print "Done"

# nn not used so far
def get_interface_timings( npart, domain, box, mixing, np, nn, ts, diff, accuracy ):
	diff = diff
	mixing = mixing
	grids = get_grid_sizes( npart, domain, box )
	timings = []
	for (cutoff_idx, P_idx, grid_idx, ewald_error) in get_interface_points( npart, domain, box ):
		# extract parameters
		cutoff = Config.cutoff_values[cutoff_idx]
		P = Config.order_values[P_idx]
		grid = grids[grid_idx]
		ewald = ewald_error[0][0]
		# load data
		output_dir = Config.get_test_timings_dir( Config.test )
		fp = os.path.join( output_dir, "Output-%s-%s-%.2f-%.2f-%d-%dx%dx%d" % (diff, mixing, ewald, cutoff, P, grid[0], grid[1], grid[2]))
		#t = extract_timings_from_log_full( fp ) # tuple (total, real, kspace)
		t = extract_timings_from_log( fp ) # tuple (total, real, kspace)
		# store
		timings.append( (t, (cutoff_idx, P_idx, grid_idx, ewald_error)) )
		#print timings[-1]
	return timings

if __name__ == "__main__":
	npart = 512000
	box = (88.08, 88.08, 88.08)
	np = 96
	nn = 12
	ts = 1000
