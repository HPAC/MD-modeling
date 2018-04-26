import sys
import os
import re

# Parsing timings
#
pppmpair  = re.compile("\[\s*(\d+)\] total time in PairLJPPPM Compute: ([0-9eE.-]+)")
pppmdisp  = re.compile("\[\s*(\d+)\] total time in PPPMDisp Compute: ([0-9eE.-]+) secs")
pppmmap   = re.compile("\[\s*(\d+)\] total time in PPPMDisp Mapping: ([0-9eE.-]+) secs")
pppmmap2  = re.compile("\[\s*(\d+)\] total time in PPPMDisp Mapping2: ([0-9eE.-]+) secs")
pppmcomm  = re.compile("\[\s*(\d+)\] total time in PPPMDisp Comm:    ([0-9eE.-]+) secs")
pppmcomm2 = re.compile("\[\s*(\d+)\] total time in PPPMDisp Comm2:    ([0-9eE.-]+) secs")
pppmfft1  = re.compile("\[\s*(\d+)\] total time in PPPMDisp FFT1:     ([0-9eE.-]+) secs")
pppmfft2  = re.compile("\[\s*(\d+)\] total time in PPPMDisp FFT2:     ([0-9eE.-]+) secs")
pppmfft2_map  = re.compile("\[\s*(\d+)\] total time in PPPMDisp FFT2_Map:     ([0-9eE.-]+) secs")
pppmfft_other = re.compile("\[\s*(\d+)\] total time in PPPMDisp FFT_Other:     ([0-9eE.-]+) secs")
pppmother = re.compile("\[\s*(\d+)\] total time in PPPMDisp Other:   ([0-9eE.-]+) secs")

pppm_neigh = re.compile("Neigh time \(%\) = ([0-9eE.-]+).*")

pppmlooptime = re.compile("Loop time of ([0-9eE.-]+)")

def extract_timings_from_log( fpath ):
	with open(fpath) as f:
		output = f.read()
		pairtimes  = sorted([(int(i), float(t)) for i,t in pppmpair.findall(output)])
		disptimes  = sorted([(int(i), float(t)) for i,t in pppmdisp.findall(output)])
		neightimes = [float(t) for t in pppm_neigh.findall(output)]
		try:
			max_idx = max(pairtimes, key=lambda t: t[1])[0]
			ptime = pairtimes[max_idx][1]
			ntime = neightimes[0]
			dtime = disptimes[max_idx][1]
			#ttime = ptime + dtime
			ttime = ptime + ntime + dtime
		except ValueError, e: # something went wrong with the test, job Exited
			print e
			print >> sys.stderr, "[Warning] Log file %s corrupt. Setting timings to NaN." % fpath
			ptime = float('NaN')
			ntime = float('NaN')
			dtime = float('NaN')
			ttime = float('NaN')
	return (ttime, ptime, dtime, ntime)

def extract_timings_from_log_full( fpath ):
	#print fpath
	with open(fpath) as f:
		output = f.read()

		pairtimes  = sorted([(int(i), float(t)) for i,t in pppmpair.findall(output)])
		disptimes  = sorted([(int(i), float(t)) for i,t in pppmdisp.findall(output)])
		maptimes   = sorted([(int(i), float(t)) for i,t in pppmmap.findall(output)])
		maptimes2  = sorted([(int(i), float(t)) for i,t in pppmmap2.findall(output)])
		commtimes  = sorted([(int(i), float(t)) for i,t in pppmcomm.findall(output)])
		commtimes2 = sorted([(int(i), float(t)) for i,t in pppmcomm2.findall(output)])
		ffttimes1  = sorted([(int(i), float(t)) for i,t in pppmfft1.findall(output)])
		ffttimes2  = sorted([(int(i), float(t)) for i,t in pppmfft2.findall(output)])
		ffttimes2M = sorted([(int(i), float(t)) for i,t in pppmfft2_map.findall(output)])
		ffttimesO  = sorted([(int(i), float(t)) for i,t in pppmfft_other.findall(output)])
		othertimes = sorted([(int(i), float(t)) for i,t in pppmother.findall(output)])

		looptime = float(pppmlooptime.findall(output)[0][0])

		try:
			max_idx = max(pairtimes, key=lambda t: t[1])[0]
			ptime = pairtimes[max_idx][1]
			dtime = disptimes[max_idx][1]
			ttime = ptime + dtime
			mtime = [m1+m2 for m1,m2 in zip(maptimes, maptimes2)][max_idx][1]
			ctime = [c1+c2 for c1,c2 in zip(commtimes, commtimes2)][max_idx][1]
			ftime = [f1+f2+f3+f4 for f1,f2,f3,f4 in zip(ffttimes1, ffttimes2,ffttimes2M, ffttimesO)][max_idx][1]
		except ValueError: # something went wrong with the test, job Exited
			print >> sys.stderr, "[Warning] Log file %s corrupt. Setting timings to NaN." % fpath
			ptime = float('NaN')
			dtime = float('NaN')
			ttime = float('NaN')
	return (ttime, ptime, dtime, mtime, ctime, ftime, looptime)


# Just for debugging and modeling purposes
def extract_timings_minmax( fpath ):
	with open(fpath) as f:
		output = f.read()
		pairtimes  = sorted([(int(i), float(t)) for i,t in pppmpair.findall(output)])
		disptimes  = sorted([(int(i), float(t)) for i,t in pppmdisp.findall(output)])
		maptimes   = sorted([(int(i), float(t)) for i,t in pppmmap.findall(output)])
		maptimes2  = sorted([(int(i), float(t)) for i,t in pppmmap2.findall(output)])
		commtimes  = sorted([(int(i), float(t)) for i,t in pppmcomm.findall(output)])
		commtimes2 = sorted([(int(i), float(t)) for i,t in pppmcomm2.findall(output)])
		ffttimes1  = sorted([(int(i), float(t)) for i,t in pppmfft1.findall(output)])
		ffttimes2  = sorted([(int(i), float(t)) for i,t in pppmfft2.findall(output)])
		ffttimes2M = sorted([(int(i), float(t)) for i,t in pppmfft2_map.findall(output)])
		ffttimesO  = sorted([(int(i), float(t)) for i,t in pppmfft_other.findall(output)])
		othertimes = sorted([(int(i), float(t)) for i,t in pppmother.findall(output)])
		
		othertimes = [t1+t2 for t1,t2 in zip(ffttimesO, othertimes)]
		timestable = [
			[p[1]+d[1] for p,d in zip(pairtimes,disptimes)],
			[t[1] for t in pairtimes],
			[t[1] for t in disptimes],
			[t[1] for t in maptimes],
			[t[1] for t in commtimes],
			[t[1] for t in ffttimes1],
			[t[1] for t in ffttimes2],
			#[t[1] for t in ffttimes2M],
			#[t[1] for t in ffttimesO],
			[t[1] for t in commtimes2],
			[t[1] for t in maptimes2]#,
			#[t[1] for t in othertimes]
		]
		print
		print "        | Total |      Pair      |  Disp |     Mapping    |      Comm      |      FFT1      |      FFT2      |      Comm2     |    Mapping2    "
		print "-----------------------------------------------------------------------------------------------------------------------------------------------"
		for npr in range(len(timestable[0])):
			cumm = 0
			print "Proc %2d:" % npr,
			for tim in range(len(timestable)):
				if tim not in (0,2):
					cumm += timestable[tim][npr]
					print " %.5f (%.3f)" % (timestable[tim][npr], cumm),
				else:
					print " %.3f " % timestable[tim][npr],
			print
		print
		try:
			min_idx, max_idx = min(pairtimes, key=lambda t: t[1])[0], max(pairtimes, key=lambda t: t[1])[0]
			min_times = (pairtimes[min_idx][1] + disptimes[min_idx][1],
					     pairtimes[min_idx][1], disptimes[min_idx][1],
						 maptimes[min_idx][1], commtimes[min_idx][1], ffttimes1[min_idx][1]) ########
			max_times = (pairtimes[max_idx][1] + disptimes[max_idx][1],
					     pairtimes[max_idx][1], disptimes[max_idx][1],
						 maptimes[max_idx][1], commtimes[max_idx][1], ffttimes1[max_idx][1])
		except ValueError: # something went wrong with the test, job Exited
			print >> sys.stderr, "[Warning] Log file %s corrupt. Setting timings to NaN." % fpath
			ptime = float('NaN')
			dtime = float('NaN')
			ttime = float('NaN')
	return min_times, max_times
