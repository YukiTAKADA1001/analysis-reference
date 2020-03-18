#/bin/bash

if test $(pwd) == /save/users/bip/simulation6/extension/acd
then
    echo Not proper working directory.
    echo Exit.
    exit 1
fi

mkdir hoge
mv initial_wat_ion.parm7 production.ncrst hoge
mv hoge ../
\rm *
mv ../hoge ./
mv hoge/* ./
rmdir hoge
cp ../job.sh ./
jsub -q PN job.sh
