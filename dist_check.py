# -*- coding: utf-8 -*-
# This script check the minimum distance of all peptide pairs in the system.
# If any pair is too close (< 10 angstrom) exit with error code.

import itertools
import sys
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import pandas as pd


TOP_FILENAME = 'initial_wat_ion.pdb'
TRAJ_FILENAME = 'initial_wat_ion.pdb'


universe = mda.Universe(TOP_FILENAME, TRAJ_FILENAME)


# ref: reference, conf: configure
segment_ref = universe.select_atoms(f'segid 0')
segment_conf = universe.select_atoms(f'segid 1')
minimum_dist = distance_array(segment_ref.atoms.positions, segment_conf.atoms.positions).min()
print(minimum_dist)

if (minimum_dist < 10):
    sys.exit(1)
