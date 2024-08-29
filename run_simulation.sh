#!/bin/bash


python run_simulation.py  --pdbid 10gs \
                          --outdir /Users/maciejwisniewski/data/QBioLipDynamics/results \
                          --protein_files /Users/maciejwisniewski/data/QBioLipDynamics/10gs_1.pdb \
                          --ligand_files /Users/maciejwisniewski/data/QBioLipDynamics/10gs_1_VWW_E.pdb \
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
                          --production_write_interval 10000 \
                          --max_steps 100000000