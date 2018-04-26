import os

class LammpsPPPMConfig:
	def __init__( self ):
		# simulation
		self.npart = None
		self.timesteps = None
		# PPPM
		self.ewald = None
		self.cutoff = None
		self.order = None
		self.grid = None
		self.mixing = None
		self.diff = None
	
	def set_file_paths( self, dir ):
		self.infile = os.path.join( 
			dir, 
			"Lammps-%s-%s-%.2f-%.2f-%d-%dx%dx%d" % (
				self.diff, 
				self.mixing, 
				self.ewald, 
				self.cutoff, 
				self.order, 
				self.grid[0], self.grid[1], self.grid[2]
			)
		)
		self.outfile = os.path.join( 
			dir, 
			"Output-%s-%s-%.2f-%.2f-%d-%dx%dx%d" % (
				self.diff, 
				self.mixing, 
				self.ewald, 
				self.cutoff, 
				self.order, 
				self.grid[0], self.grid[1], self.grid[2]
			)
		)

	def write_config_file( self, output_path ):
		if not self.npart or \
		   not self.timesteps or \
		   not self.ewald or \
		   not self.cutoff or \
		   not self.order or \
		   not self.grid or \
		   not self.mixing or \
		   not self.diff:
			raise Exception
		
		lines = []
		lines.append("echo            screen")
		lines.append("units           lj")
		lines.append("atom_style      atomic")
		lines.append("")
		lines.append("read_restart        config/e0_%d.rs" % self.npart) ######
		lines.append("reset_timestep      0") ###
		lines.append("")
		lines.append("pair_style     lj/long/coul/long long off %.2f" % self.cutoff)
		lines.append("pair_coeff     1 1 1.0 1.0")
		lines.append("kspace_style   pppm/disp  1e-3")    ####### Do something with this?
		lines.append("kspace_modify  mesh/disp %d %d %d" % self.grid)
		lines.append("kspace_modify  mix/disp %s" % self.mixing)
		lines.append("kspace_modify  gewald/disp %.2f" % self.ewald)
		lines.append("kspace_modify  order/disp %d" % self.order)
		lines.append("kspace_modify  diff %s" % self.diff)
		lines.append("")
		lines.append("fix             nhvt all nvt temp 0.7 0.7 10")
		lines.append("thermo          100")
		lines.append("fix             rcent all recenter NULL NULL 0.5 units fraction")
		lines.append("")
		lines.append("run %d" % self.timesteps)
		lines.append("")

		with open( output_path, "w" ) as f:
			f.write( "\n".join(lines) )
