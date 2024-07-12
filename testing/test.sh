#!/bin/bash
#
# Install Julia
curl -fsSL https://install.julialang.org | sh
# Run the tests
julia runtests.jl ./input_files/water_box.inp \
                  ./input_files/ieee_signaling.inp \
                  ./input_files/mixture.inp \
                  ./input_files/spherical.inp \
                  ./input_files/bilayer.inp

# check if output files are properly generated in a failed run
./test_failed.sh ./input_files/water_box_failed.inp packmol.log "FORCED" 
