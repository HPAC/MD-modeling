class Simulation:
    def __init__( self, name, npart, domain, box, mixing, restart, \
                  timesteps, nnodes, nprocesses, \
                  accuracy, accuracy_real = None, accuracy_kspace = None ):
        self.name = name
        self.npart = npart
        self.domain = domain
        self.box = box
        self.mixing = mixing
        #
        self.restart = restart
        #
        self.timesteps = timesteps
        self.nnodes = nnodes
        self.nprocesses = nprocesses
        self.diff_mode = None
        self.accuracy = accuracy
        self.accuracy_real   = accuracy_real
        self.accuracy_kspace = accuracy_kspace
