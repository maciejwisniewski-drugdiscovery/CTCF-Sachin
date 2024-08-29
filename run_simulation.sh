#!/bin/bash


python run_simulation.py  --pdbid 5T00 \
                          --outdir /Users/maciejwisniewski/data/SachinCTCF/dynamics/5T00_1 \
                          --complex_file /Users/maciejwisniewski/data/SachinCTCF/rawPDB/5T00_1.pdb \
                          --platform CPU \
                          --add_solvate \
                          --padding 1.0 \
                          --ionic_strength 0.15 \
                          --friction_value 1.0 \
                          --timestep_value 1.0 \
                          --initial_temp_value 50.0 \
                          --heating_step 100 \
                          --heating_write_interval 100 \
                          --temp_value 300 \
                          --production_write_interval 1000 \
                          --max_steps 100000000 \

# 100ns simulation with interval every 10 ps

# 100 000 000 / 10 000 => 10000timeshots