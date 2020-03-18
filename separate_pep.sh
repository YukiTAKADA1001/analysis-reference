#/bin/bash



readonly BASEDIR=/save/users/bip/simulation6/extension/acd

cd $BASEDIR

sim_dirs=$(ls -d 20*)

for dir in $sim_dirs
do
    cd $dir

    if ! test -e extension_pep1.ncdf
    then
        python -c "
import MDAnalysis as mda

BASE_DIR = '/save/users/bip/simulation6/'
PATH_TOP = BASE_DIR + 'main/acd/analysis/repo/stripped.initial_wat_ion.parm7'

universe = mda.Universe(PATH_TOP, 'extension_dry.nc')


protein_pep1 = universe.select_atoms('resid 1:42')

with mda.Writer('extension_pep1.ncdf', protein_pep1.n_atoms) as W:
    for ts in universe.trajectory:
        W.write(protein_pep1)


protein_pep2 = universe.select_atoms('resid 43:84')

with mda.Writer('extension_pep2.ncdf', protein_pep2.n_atoms) as W:
    for ts in universe.trajectory:
        W.write(protein_pep2)
"
        echo $(pwd): OK
    fi

    cd $BASEDIR
done
