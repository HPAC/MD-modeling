import sys

import General.Config as Config
from Lammps.PPPMDisp.Auxiliary.KSpaceGrids import get_grid_sizes

def generate_kspace_error_config_file( test, output=sys.stdout ):
    print >> output, test.domain[0], test.domain[1], test.domain[2]
    print >> output, test.npart
    if test.accuracy:
        print >> output, test.accuracy
    else:
        print >> output, test.accuracy_kspace
    print >> output, Config.min_ewald, Config.max_ewald+Config.d_ewald/10, Config.d_ewald
    print >> output, Config.order_values[0], Config.order_values[-1]
    grids = get_grid_sizes( test.npart, test.domain, test.box )
    print >> output, len(grids)
    for g in grids:
        print >> output, "%4d %4d %4d" % g
