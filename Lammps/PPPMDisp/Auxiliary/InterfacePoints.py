import os
import re

import General.Config as Config
from Lammps.PPPMDisp.Auxiliary.KSpaceGrids import get_grid_sizes
from Lammps.PPPMDisp.ErrorBounds import get_accurate_grid

def get_interface_points( npart, domain, box, width=1 ):
    error_grid = get_accurate_grid( npart, domain, box, None )
    grids = get_grid_sizes( npart, domain, box )

    for cutoff_idx, cutoff in enumerate(Config.cutoff_values):
        for P_idx, P in enumerate(Config.order_values):
            for grid_idx, grid_side in enumerate(grids):
                if error_grid[cutoff_idx][P_idx][grid_idx] and ( \
                        (cutoff_idx < width or not error_grid[cutoff_idx-width][P_idx][grid_idx]) and \
                        (P_idx < width      or not error_grid[cutoff_idx][P_idx-width][grid_idx]) and \
                        (grid_idx < width   or not error_grid[cutoff_idx][P_idx][grid_idx-width])
                    ):
                    yield (cutoff_idx, P_idx, grid_idx, error_grid[cutoff_idx][P_idx][grid_idx])
