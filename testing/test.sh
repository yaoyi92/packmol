#!/bin/bash
#
# Install Julia
curl -fsSL https://install.julialang.org | sh
# Run the tests
julia runtests.jl water_box.inp \
                  ieee_signaling.inp