#!/bin/bash

cpptraj << EOF

readdata ../rem.log
remlog rem.log stats

EOF
