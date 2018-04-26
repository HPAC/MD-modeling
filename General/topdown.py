import numpy as np
import scipy
from scipy.optimize import curve_fit

def linear( x, a, b ):
	return a*x + b

def fit( x, y, f ):
	popt, pcov = curve_fit(f, x, y)
	return popt

def error( xs, ys, f, *args ):
	max_error = 0
	for x,y in zip(xs,ys):
		y_est = f( x, *args )
		err = abs((y-y_est)/y)
		#print "*** Inv error:", err
		if err > max_error:
			#print err
			max_error = err
	return max_error

def split_error( x, y, f, i ):
	args = fit( x[:i], y[:i], f )
	err_left = error( x[:i], y[:i], f, *args)
	args = fit( x[i:], y[i:], f )
	err_right = error( x[i:], y[i:], f, *args)
	#print "* Error: (%d)" % len(x), err_left, err_right
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
		#print breakpoint, best_so_far
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


if __name__ == "__main__":
	with open("kspace-times.txt", "r") as f:
		x, y = [], []
		for line in f.readlines():
			gx,gy,gz, t = line.strip().split()
			gx = int(gx)
			gy = int(gy)
			gz = int(gz)
			t  = float(t)
			x.append(gx*gy*gz)
			y.append(t)

	np_x = np.array(x)
	np_y = np.array(y)
	print top_down( np_x, np_y, linear )
