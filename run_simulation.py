import os
import argparse
from dynamic.dynamics import OpenMMSimulation


def run(config):
    '''
    Run the Molecular Dynamics Simulation Pipeline
    :param config: Dictionary of Parameters for the MD Simulation.
    '''

    os.makedirs(os.path.join(config['outdir']), exist_ok=True)
    openmm = OpenMMSimulation(config)


    return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Molecular Dynamic Simulation with OpenMM")
    parser.add_argument('--pdbid', type=str, required=True, help='Specific PDB ID... please... ;-;')
    parser.add_argument('--outdir', type=str, required=True, help='Path to the Output Directory.')
    parser.add_argument('--complex_file',nargs='+', type=str, required=True, help='Complex PDB Structure Filepath.')
    parser.add_argument('--platform',type=str,required=False,default='CUDA',help='Define if you want the OpenCL or CUDA support')
    # Molecular Dynamics Simulation Forcefields

    parser.add_argument('--pff',type=str,required=False,default='amber14-all.xml',help='Define the type of Protein Force Field, eg. \'amber14-all.xml\'')
    parser.add_argument('--wff',type=str,required=False,default='amber14/tip3pfb.xml',help='Define the type of Water Force Field, eg. \'amber14/tip3pfb.xml\'')
    parser.add_argument('--nff',type=str,required=False,default='amber14/DNA.OL15.xml',help='Define the type of Water Force Field, eg. \'amber14/tip3pfb.xml\'')
    parser.add_argument('--mff',type=str,required=False,default='gaff-2.11',help='Define the type of Small Molecule Force Field, eg. \'gaff-2.11\'')
    parser.add_argument("--water_model", default="tip3p",choices=["tip3p", "spce", "tip4pew", "tip5p", "swm4ndp"],help="Water model for solvation")
    # Molecular Dynamics Box Preparation
    parser.add_argument('--use_modeller', action='store_true',help='Define if you want to rebuild the Protein Structures with the modeller.')
    parser.add_argument('--add_solvate', action='store_true', help="Add Solvent Molecules to the Complex Model")
    parser.add_argument('--add_prot_const', action='store_true', help="Add Protein Constraints to the Complex Model")
    parser.add_argument("--padding", type=float, default=1.0, help="Padding for solvent box [nm]")
    parser.add_argument("--ionic_strength", type=float, default=0.15, help="Ionic strength for solvation [molar]")
    # Molecular Dynamics Simulation Parameters
    parser.add_argument("--friction_value", type=float, default=1.0,
                        help="Friction value [1/ps]")
    parser.add_argument("--timestep_value", type=float, default=1.0,
                        help="TimeStep value during simulation [fs].\nEach frame corresponds to a change after specified <timestep_value>.")
    parser.add_argument("--heating_step", type=int, default=100,
                        help="How many steps to WarmUp Simulation by 10K")
    parser.add_argument("--initial_temp_value", type=float, default=50.0,
                        help="Initial Temperature Value during WarmUp Simulation Run [K]")
    parser.add_argument("--heating_write_interval", type=int, default=100,
                        help="Define how often information should be saved in the WarmUp Data Reporter in terms of the number of frames. ")
    parser.add_argument("--temp_value", type=float, default=300.0,
                        help="Temperature Value during Production Simulation Run [K]")
    parser.add_argument("--production_write_interval", type=int, default=10000,
                        help="Define how often information should be saved in the Production Data Reporter in terms of the number of frames. ")
    parser.add_argument("--max_steps", type=int, default=100000000,
                        help="Define how many steps during Production Simulation [fs] [1000fs = 1ps] [1000ps = 1ns]")
    parser.add_argument('--non_bonded_reporter', action='store_true',
                        help="Add Non-Bonded Forces Reporter to the Simulation. Actually only works for singe ligand simulation, okey?")

    # Parse the command line arguments
    args = parser.parse_args()

    # Asserts
    assert os.path.exists(args.complex_file), f'{args.complex_file}. Complex File does not exist!'

    config = {
        'pdbid': args.pdbid,
        'outdir': args.outdir,
        'complexFile': args.complex_file,
        'platformName': args.platform,
        # Force Fields Parameters
        'nucleicFF': args.nff,
        'proteinFF': args.pff,
        'waterFF': args.wff,
        'ligandFF': args.mff,
        'waterModel': args.water_model,
        # MD Box Parameters
        'modellerRebuilder': args.use_modeller,
        'addSolvent': args.add_solvate,
        'proteinConstraint': args.add_prot_const,
        'paddingValue': args.padding,
        'ionicStrengthValue': args.ionic_strength,
        # Simulation Parameters
        'frictionValue': args.friction_value,
        'timeStepValue': args.timestep_value,
        'NonBondedReporter': args.non_bonded_reporter,
        # WarmUp Parameters
        'initialTemperatureValue': args.initial_temp_value,
        'heatingStep': args.heating_step,
        'heatingWriteInterval': args.heating_write_interval,
        # Production Parameters
        'temperatureValue': args.temp_value,
        'maxSteps': args.max_steps,
        'productionWriteInterval': args.production_write_interval,
    }
    run(config)