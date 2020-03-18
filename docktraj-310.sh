#!/bin/bash

# This script extract frames only at 310K from T-REMD trajectories.

cpptraj << EOF

parm ../prmtop

trajin ../mdcrd.rep.001 remdtraj remdtrajtemp 310.0
trajout ./mdcrd.310

EOF
