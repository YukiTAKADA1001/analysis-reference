# -*- coding: utf-8 -*-
import random
import warnings

import MDAnalysis as mda
import MDAnalysis.analysis.align

warnings.simplefilter('ignore')

BASEDIR = '/save/users/bip/simulation6/main/acd'

TOP_NAME = f'{BASEDIR}/repo/monomer.parm7'
TRAJ_NAME = f'{BASEDIR}/repo/mdcrd.310.nc'

universe = mda.Universe(TOP_NAME, TRAJ_NAME, topology_format='PARM7')
n_frames = universe.trajectory.n_frames

for i in range(2):
    frame = random.randint(0, n_frames - 1)
    with mda.Writer(f"sample{i}.pdb") as pdb:
        for ts in universe.trajectory[frame:frame+1]:
            pdb.write(universe)
