#/bin/bash

sim_dirs=$(ls -d 20*)

for sim_dir in $sim_dirs
do
    cd $sim_dir
    if ! test -e extension.nc
    then
        \cp ../job.sh ./
        jsub -q PN job.sh
    fi

    cd ../
done
