#!/bin/bash
#
# Install Julia
curl -fsSL https://install.julialang.org | sh
# Run the tests
julia runtests.jl water_box.inp \
                  ieee_signaling.inp 

# check if output files are properly generated in a failed run
./test_failed.sh water_box_failed.inp packmol.log "FORCED" 
