class LSF_Config:
	def __init__(self):
		self.outputfile = ""
		self.email = ""
		self.group = ""
		self.jobname = ""
		self.time = ""    # hh:mm
		self.memory = ""  # MBs
		self.nthreads = ""
		self.parallelism_type = ""
		self.arch_string = ""
		self.command = ""

	def __repr__(self):
		#not self.outputfile or \
		if not self.time or \
		   not self.memory or \
		   not self.nthreads or \
		   not self.parallelism_type or \
		   not self.command:
			raise Exception # improve

		lines = []
		if self.outputfile:
			lines.append("#BSUB -o %s" % self.outputfile)
		else:
			lines.append("#BSUB -o /dev/null/")
		if self.email:
			lines.append("#BSUB -B")
			lines.append("#BSUB -N")
			lines.append("#BSUB -u %s" % self.email)
		if self.group:
			lines.append("#BSUB -P %s" % self.group)
		if self.jobname:
			lines.append("#BSUB -J %s" % self.jobname)
		lines.append("#BSUB -W %s" % self.time)
		lines.append("#BSUB -M %s" % self.memory)
		lines.append("#BSUB -n %s" % self.nthreads)
		lines.append("#BSUB -a \"%s\"" % self.parallelism_type)
		if self.arch_string:
			lines.append("#BSUB -R \"select[model==%s]\"" % self.arch_string)
		lines.append("#BSUB -x")
		lines.append("")
		lines.append(self.command)
		lines.append("")
		lines.append("echo 'Done'")
		lines.append("")

		return "\n".join(lines)

	def write_jobscript(self, filepath):
		with open(filepath, "w") as f:
			f.write( str(self) )

if __name__ == "__main__":
	config = LSF_Config()
	config.outputfile = "LSF_test"
	config.group = "AICES"
	config.jobname = "PPPM"
	config.time = "1:00"
	config.memory = "8000"
	config.nthreads = 8
	config.parallelism_type = "openmpi"
	config.arch_string = "Harpertown"
	config.command = "command"

	print config

