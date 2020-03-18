#!/bin/bash

readonly PARALLEL=15
readonly LOOP=10

steps=$(for j in $(seq 1 $LOOP); do echo job.sh; done)

for i in $(seq 1 $PARALLEL)
do
    jsub -q PN --stepany $(echo $steps)
done
