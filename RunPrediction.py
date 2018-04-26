import sys

import General.Config as Config
from TestSuite import TestSuite
from Lammps.PPPMDisp.Sampling import run_samples
from Lammps.PPPMDisp.Empirical import run_interface_timings
from Lammps.PPPMDisp.Prediction import get_prediction_with_statistics


if len(sys.argv) != 4:
	print >> sys.stderr, "Usage: %s TestName [ad|ik] action" % sys.argv[0]
	sys.exit(-1)
	
# Test info
test = TestSuite[sys.argv[1]]
Config.test = test
test.diff_mode = sys.argv[2]


domain = test.domain
box    = test.box
npart  = test.npart
mixing = test.mixing

diff   = test.diff_mode
acc    = test.accuracy

nnodes = test.nnodes
nprocs = test.nprocesses
ts     = test.timesteps

if sys.argv[3] == "--run-samples":
	# run samples for performance prediction
	run_samples( test, diff )

elif sys.argv[3] == "--run-interface":
	# run experiments for the points in the interface
	# these are empirical values to measure the quality
	# of the prediction
	run_interface_timings( npart, domain, box, mixing, nprocs, nnodes, ts, diff, acc )

elif sys.argv[3] == "--predict":
	# give me a prediction (with a statistical measure of quality)
	# est -> estimated/predicted
	# emp -> empirical timings, i.e., actually measured for quality control purposes
	( (cutoff, P, grid, ewald, t), stats, (est_tims, emp_tims) ) = \
		get_prediction_with_statistics( npart, domain, box, mixing, nprocs, nnodes, ts, diff, acc )

	print "Prediction"
	print "----------"
	print "Ewald:  %.2f" % ewald
	print "Cutoff: %.2f" % cutoff
	print "Order:  %d" % P
	print "Grid:  ", grid
	print "Estimated time: ", t
	print

	#print "Stats"
	#print "-----"
	#print "Min difference:       " , stats[0]
	#print "Max difference:       " , stats[1]
	#print "Mean difference:      " , stats[2]
	#print "Median of differences:" , stats[3]
	#print "STD differences:      " , stats[4]
	#print

	print "Ranking"
	print "-------"
	est_tims_sorted = sorted( zip( est_tims, range(len(est_tims)) ), key=lambda t: t[0][0][0] ) # t[timing[times[totaltime]]]
	emp_tims_sorted = sorted( zip( emp_tims, range(len(emp_tims)) ), key=lambda t: t[0][0][0] ) # t[timing[times[totaltime]]]
	#est_tims_sorted = sorted( zip( est_tims, range(len(est_tims)) ), key=lambda t: t[0][0][1]+t[0][0][2] ) # t[timing[times[totaltime]]]
	#emp_tims_sorted = sorted( zip( emp_tims, range(len(emp_tims)) ), key=lambda t: t[0][0][1]+t[0][0][2] ) # t[timing[times[totaltime]]]
	print " Rank   Empir. Time | Rank   Estim. Time"
	print "----------------------------------------"
	for est, emp in zip(est_tims_sorted, emp_tims_sorted)[:20]:
		t_emp_r = emp[0][0][0]
		#t_emp_r = emp[0][0][1] + emp[0][0][2]
		t_emp = emp_tims[est[1]][0][0]
		#t_emp = emp_tims[est[1]][0][1] + emp_tims[est[1]][0][2]
		t_est = est[0][0][0]
		#t_est = est[0][0][1] + est[0][0][2]
		#print " %3d %14.3f | %3d %14.3f | (%.2f%%)" % (emp[1], emp[0][0][0], est[1], est[0][0][0], emp_tims[est[1]][0][0]/est[0][0][0] )
		print " %3d %14.3f | %3d %14.3f | (%+.2f%%) -- (%+6.2f%%, %+6.2f%%, %+6.2f%%)" % (emp[1], t_emp_r, est[1], t_est, \
				(1 - t_emp/t_est) * 100, \
				(1 - emp_tims[est[1]][0][1]/est[0][0][1]) * 100, \
				(1 - emp_tims[est[1]][0][3]/est[0][0][3]) * 100, \
				(1 - emp_tims[est[1]][0][2]/est[0][0][2]) * 100 )
	# Total loop time as reported by LAMMPS (not really, there is still some averaging)
	#loop_times_sorted = sorted( zip( emp_tims, range(len(emp_tims)) ), key=lambda t: t[0][0][6] )
	#for lt in loop_times_sorted:
		#print " %3d %14.3f " % (lt[1], lt[0][0][6])

else:
	print "Valid options: --run-samples, --run-interface, and --predict"
