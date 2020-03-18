#!/bin/bash

sim_dirs=$(ls -d 20*)

for dir in $sim_dirs
do
    cd $dir

    if test -e extension.nc
    then
        cpptraj -p initial_wat_ion.parm7 <<EOF
trajin extension.nc
strip :WAT,Cl-,Na+
trajout extension_dry.nc
run
exit
EOF
    fi

    cd ../
done
