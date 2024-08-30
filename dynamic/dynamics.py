import os
import pandas as pd
import numpy as np
from rdkit import Chem
from openbabel import pybel
from openff.toolkit.topology import Molecule
import pdbfixer
from openmm import app, unit
from openmm import Platform, XmlSerializer, LangevinIntegrator, NonbondedForce
from openmmforcefields.generators import SystemGenerator

class ForceReporter(object):
    def __init__(self, file, reportInterval):
        if not os.path.exists(file):
            self._out = open(file, 'w')
            self._out.write('Step\tCoulomb Interaction Energy\tLJ Interaction Energy\n')
        else:
            self._out = open(file,'a')

        self._reportInterval = reportInterval
    def __del__(self):
        self._out.close()
    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, False, False, True, False, None)
    def energy(self,simulation,protein_coulomb_scale,protein_lj_scale,ligand_coulomb_scale,ligand_lj_scale,solvent_coulomb_scale,solvent_lj_scale):
        simulation.context.setParameter("protein_coulomb_scale", protein_coulomb_scale)
        simulation.context.setParameter("protein_lj_scale", protein_lj_scale)
        simulation.context.setParameter("ligand_coulomb_scale", ligand_coulomb_scale)
        simulation.context.setParameter("ligand_lj_scale", ligand_lj_scale)
        simulation.context.setParameter("solvent_coulomb_scale", solvent_coulomb_scale)
        simulation.context.setParameter("solvent_lj_scale", solvent_lj_scale)
        return simulation.context.getState(getEnergy=True, groups={0}).getPotentialEnergy()
    def report(self, simulation, state):
        # Coulomb Energy
        protein_ligand_coulomb_energy = self.energy(simulation,1,0,1,0,0,0)
        protein_coulomb_energy = self.energy(simulation,1,0,0,0,0,0)
        ligand_coulomb_energy = self.energy(simulation,0,0,1,0,0,0)

        # Lennard-Jones Energy
        protein_ligand_lj_energy = self.energy(simulation,0,1,0,1,0,0)
        protein_lj_energy = self.energy(simulation,0,1,0,0,0,0)
        ligand_lj_energy = self.energy(simulation,0,0,0,1,0,0)

        # Interaction

        coulomb_energy = protein_ligand_coulomb_energy - protein_coulomb_energy - ligand_coulomb_energy
        coulomb_energy = coulomb_energy.value_in_unit(unit.kilocalories_per_mole)

        lj_energy =  protein_ligand_lj_energy - protein_lj_energy - ligand_lj_energy
        lj_energy = lj_energy.value_in_unit(unit.kilocalories_per_mole)

        self._out.write(f'{str(simulation.currentStep)}\t{str(coulomb_energy)}\t{str(lj_energy)}\n')
        self._out.flush()
class OpenMMSimulation():
    def __init__(self,config):
        self.homedir = os.getcwd()
        self.pdbid = config['pdbid']

        assert os.path.exists(config['outdir']), "Output Directory does not exist!"
        self.outdir = config['outdir']

        ###############
        # INPUT FILES #
        ###############

        # Specify protein files
        self.complex_file = config['complexFile']
        assert len(self.complex_file) != 0, 'There are no protein files in the OutDir!'

        ################
        # OUTPUT FILES #
        ################

        # Topology Files
        # Initial Topology File
        self.InitialTopologyFile = os.path.join(self.outdir, 'InitialTopology.pdb')
        # Topology after WarmUp File
        self.WarmedTopologyFile = os.path.join(self.outdir, 'WarmedTopology.pdb')

        # Systems
        self.InitialSystemFile = os.path.join(self.outdir, 'InitialSystem.xml')
        self.WarmUpSystemFile = os.path.join(self.outdir, 'WarmUpSystem.xml')
        self.ProductionSystemFile = os.path.join(self.outdir, 'ProductionSystem.xml')

        # Output Trajectories
        self.WarmUpTrajectoryFile = os.path.join(self.outdir, 'WarmUpTrajectory.xtc')
        self.ProductionTrajectoryFile = os.path.join(self.outdir, 'ProductionTrajectory.xtc')

        # Output State Data Reports
        self.WarmUpStateDataReporterFile = os.path.join(self.outdir, 'WarmUpStateDataReport.txt')
        self.ProductionStateDataReporterFile = os.path.join(self.outdir, 'ProductionStateDataReport.txt')

        # Output Interaction Force Data Reporter
        self.WarmupInteractionForceDataReporterFile = os.path.join(self.outdir, 'WarmUpInteractionForceDataReport.txt')
        self.ProductionInteractionForceDataReporterFile = os.path.join(self.outdir,
                                                                       'ProductionInteractionForceDataReport.txt')

        # Output CheckPoints
        self.WarmedCheckpointFile = os.path.join(self.outdir, 'WarmedCheckpoint.chk')
        self.WarmedStateFile = os.path.join(self.outdir, 'WarmedState.state')
        self.ProductionCheckpointFile = os.path.join(self.outdir, 'ProductionCheckpoint.chk')
        self.ProductionStateFile = os.path.join(self.outdir, 'ProductionState.state')

        ##############
        # Parameters #
        ##############
        self.config = config

        # Number of WarmUp Steps
        self.heating_steps = (self.config['temperatureValue'] - self.config['initialTemperatureValue']) * self.config['heatingStep']
        # Number of Prtoduction Steps
        self.steps_to_run = self.how_many_steps()

        #############
        # Variables #
        #############

        # List of proteins, read by OpenMM and PDBFixer
        self.complex = None
        # Molecular Dynamic Simulation Box
        self.complex_model = None
        # Molecular Dynamic Simulation System
        self.system = None
        # Molecular Dynamic Simulation ForceFields
        self.forcefield = None
        # Hardware Platform
        self.platform = self.get_platform()

    def how_many_steps(self):
        ''' Get number of steps left to the End of the Simulation '''
        if os.path.exists(self.ProductionCheckpointFile):
            df = pd.read_csv(self.ProductionStateDataReporterFile,sep='\t')
            last_step = df['Step'].tolist()[-1]
            steps_to_run = int(self.config['maxSteps']) - (int(last_step) - self.heating_steps)
        else:
            steps_to_run = self.config['maxSteps']
        return steps_to_run
    def get_platform(self):

        platform = Platform.getPlatformByName(self.config['platformName'])

        print('Using platform', platform.getName())

        # if it's GPU platform set the precision to mixed
        if platform.getName() == 'CUDA' or platform.getName() == 'OpenCL':
            platform.setPropertyDefaultValue('Precision', 'mixed')
            platform.setPropertyDefaultValue('CudaDeviceIndex','0') # one GPU, ok?
            print('Set precision for platform', platform.getName(), 'to mixed')

        return platform
    def parse_complex(self, ignore_missing_residues=True, ignore_terminal_missing_residues=True):

        # Initialize PDBFixer Class
        fixer = pdbfixer.PDBFixer(self.complex_file)
        # Co-Crystallized Ligands are unknown for PDBFixer
        fixer.removeHeterogens()
        # Identify Missing Residues, needed later for identification of Missing Atoms
        fixer.findMissingResidues()
        # If Missing Terminal Residues shall be ignored, remove them from the dictionary
        if ignore_terminal_missing_residues:
            # Iterate over the Missing Residues
            chains = list(fixer.topology.chains())
            keys = fixer.missingResidues.keys()
            for key in list(keys):
                chain = chains[key[0]]
                if key[1] == 0 or key[1] == len(list(chain.residues())):
                    del fixer.missingResidues[key]

        # If All Missing Residues shall be ignored, clear the dictionary.
        if ignore_missing_residues:
            fixer.missingResidues = {}

        # Find Non-Standard Residues
        fixer.findNonstandardResidues()
        # Replace Non-Standard Residues with Standard Ones
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        print("Adding missing atoms...")
        fixer.addMissingAtoms()
        print("Adding missing hydrogens...")
        fixer.addMissingHydrogens(7)
        # Add Protein PDBFixer Class to MD Simulation Protein List
        self.complex = fixer
        return self
    def create_complex(self):
        '''
        Create an Initial Topology Complex File based on already loaded Proteins and Ligands.
        '''
        complex_model = app.Modeller(self.complex.topology,self.complex.positions)
        self.complex_model = complex_model
        return self
    def generate_initial_system(self):
        ''' Generate a System '''

        forcefield_kwargs = {
            'constraints': None,
            'rigidWater': True,
            'removeCMMotion': False,
            'hydrogenMass': 4*unit.amu
        }

        print('\t[**] Prepare Protein and Generate OpenMM Proteins...')
        self. parse_complex()
        print('\t[**] Merge OpenMM Proteins and Molecules...')
        self.create_complex()
        print('\t[**] Generate the system of Protein - Ligand Complex')

        system_generator = SystemGenerator(
            forcefields=[self.config['proteinFF'], self.config['nucleicFF'],self.config['waterFF']],
            small_molecule_forcefield=self.config['ligandFF'],
            molecules=None,
            forcefield_kwargs=forcefield_kwargs
        )

        if self.config['addSolvent']:
            print('\t[**] Adding solvent')
            self.complex_model.addSolvent(
                system_generator.forcefield,
                model = self.config['waterModel'],
                padding = self.config['paddingValue'] * unit.nanometers,
                ionicStrength = self.config['ionicStrengthValue'] * unit.molar,
            )
            print('\tSystem has %d atoms' % self.complex_model.topology.getNumAtoms())
        if self.config['proteinConstraint']:
            print('\t[**] Constrain Protein Movement...')
            for atom in self.complex_model.topology.atoms():
                if atom.residue.name not in ('HOH', 'Cl', 'Na', 'UNK'):
                    # By setting Protein Particle Mass to 0 we constrain That Particle Movement
                    self.system.setParticleMass(atom.index, 0 * unit.amu)

        # Save Generated Topology
        print('\t[**] Save Generated Topology')
        with open(self.InitialTopologyFile, 'w') as outfile:
            app.PDBFile.writeFile(self.complex_model.topology, self.complex_model.positions, outfile)

        # Save Generated System
        self.system = system_generator.create_system(self.complex_model.topology, molecules=None)
        print('\t[**] Save Initial System')
        with open(self.InitialSystemFile, 'w') as outfile:
            outfile.write(XmlSerializer.serialize(self.system))

        # Define Atoms of Solvent, Ligand and Protein for Calculating L-J and Electriostatics Forces.

        self.solvent_atoms = set([a.index for a in self.complex_model.topology.atoms() if a.residue.name in ('HOH', 'Cl', 'Na')])
        self.ligand_atoms = set([a.index for a in self.complex_model.topology.atoms() if a.residue.name in ('UNK')])
        self.protein_atoms = set([a.index for a in self.complex_model.topology.atoms() if a.residue.name not in ('UNK', 'HOH', 'Cl', 'Na')])
    def load_model(self):
        '''
        Load an Initial Topology Complex File.
        '''
        temp = app.pdbfile.PDBFile(self.InitialTopologyFile)
        self.complex_model = app.Modeller(temp.topology, temp.positions)

        # Define Atoms of Solvent, Ligand and Protein for Calculating L-J and Electriostatics Forces.

        self.solvent_atoms = set([a.index for a in self.complex_model.topology.atoms() if a.residue.name in ('HOH', 'Cl', 'Na')])
        self.ligand_atoms = set([a.index for a in self.complex_model.topology.atoms() if a.residue.name in ('UNK')])
        self.protein_atoms = set([a.index for a in self.complex_model.topology.atoms() if a.residue.name not in ('UNK', 'HOH', 'Cl', 'Na')])
    def load_system(self):
        '''Load the most current System File.'''
        if os.path.exists(self.InitialSystemFile):
            with open(self.InitialSystemFile) as input:
                self.system = XmlSerializer.deserialize(input.read())
        else:
            print('There are no System Files')
            self.generate_initial_system()
    def warmup(self):
        ''' MD Simulation to WarmUp Studied Complex from Initial Temperature to 300K. '''

        # Define Integrator for Simulation
        integrator = LangevinIntegrator(self.config['initialTemperatureValue'] * unit.kelvin,
                                        self.config['frictionValue'] / unit.picoseconds,
                                        self.config['timeStepValue'] * unit.femtoseconds)

        # Load Model if not Loaded (?)
        if not self.complex_model:
            self.load_model()

        # Load System if not Loaded (?)
        if not self.system:
            self.load_system()

        # Initialize NonBonded Forces Reporter
        if self.config['NonBondedReporter']:
            for force in self.system.getForces():
                if isinstance(force, NonbondedForce):
                    force.setForceGroup(0)
                    force.addGlobalParameter("protein_coulomb_scale", 1)
                    force.addGlobalParameter("protein_lj_scale", 1)
                    force.addGlobalParameter("ligand_coulomb_scale", 1)
                    force.addGlobalParameter("ligand_lj_scale", 1)
                    force.addGlobalParameter("solvent_coulomb_scale", 1)
                    force.addGlobalParameter("solvent_lj_scale", 1)
                    for i in range(force.getNumParticles()):
                        charge, sigma, epsilon = force.getParticleParameters(i)

                        force.setParticleParameters(i, 0, 0, 0)
                        if i in self.protein_atoms:
                            force.addParticleParameterOffset("protein_coulomb_scale", i, charge, 0, 0)
                            force.addParticleParameterOffset("protein_lj_scale", i, 0, sigma, epsilon)
                        elif i in self.ligand_atoms:
                            force.addParticleParameterOffset("ligand_coulomb_scale", i, charge, 0, 0)
                            force.addParticleParameterOffset("ligand_lj_scale", i, 0, sigma, epsilon)
                        else:
                            force.addParticleParameterOffset("solvent_coulomb_scale", i, charge, 0, 0)
                            force.addParticleParameterOffset("solvent_lj_scale", i, 0, sigma, epsilon)
                    for i in range(force.getNumExceptions()):
                        p1, p2, chargeProd, sigma, epsilon = force.getExceptionParameters(i)
                        force.setExceptionParameters(i, p1, p2, 0, 0, 0)
                else:
                    force.setForceGroup(2)

        warmup_simulation = app.Simulation(self.complex_model.topology, self.system, integrator, platform=self.platform)
        warmup_simulation.context.setPositions(self.complex_model.positions)

        # Initial Minimization
        print('\t[**] WarmUp Energy Minimization')
        warmup_simulation.minimizeEnergy()

        # Add Trajectory Reporter
        warmup_simulation.reporters.append(
            app.xtcreporter.XTCReporter(
                file=self.WarmUpTrajectoryFile,
                reportInterval=self.config['heatingWriteInterval'],
                enforcePeriodicBox=False,
                append=os.path.exists(self.WarmUpTrajectoryFile)
            )
        )

        # Add State Reporter
        warmup_simulation.reporters.append(
            app.StateDataReporter(
                file=self.WarmUpStateDataReporterFile,
                reportInterval=self.config['heatingWriteInterval'],
                step=True,
                potentialEnergy=True,
                kineticEnergy=True,
                temperature=True,
                volume=True,
                progress=True,
                remainingTime=True,
                speed=True,
                totalSteps=self.heating_steps,
                separator="\t",
                append=os.path.exists(self.WarmUpStateDataReporterFile)
            )
        )

        # Add Non-Bonded Forces Reporter such as L-J and vdW
        if self.config['NonBondedReporter']:
            warmup_simulation.reporters.append(
                ForceReporter(
                    file=self.WarmupInteractionForceDataReporterFile,
                    reportInterval=self.config['heatingWriteInterval']
                )
            )

        # Set the Velocities (Speed of Particles) based on Initial Temperature

        warmup_simulation.context.setVelocitiesToTemperature(self.config['initialTemperatureValue'] * unit.kelvin)

        # Warming Up!
        print('\t[**] Warming Up...')
        for i in range(int(self.config['temperatureValue'] - self.config['initialTemperatureValue'])):
            print(f"\t\t[***] Actual Temperature: {str(self.config['temperatureValue'] - self.config['initialTemperatureValue'])}")
            warmup_simulation.step(self.config['heatingStep'])
            current_temperature = (self.config['initialTemperatureValue'] + int(i)) * unit.kelvin
            integrator.setTemperature(current_temperature)

        # Save State and Checkpoint of Warmed Simulation
        # The Temperature should be sth around 300K
        warmup_simulation.saveState(self.WarmedStateFile)
        warmup_simulation.saveCheckpoint(self.WarmedCheckpointFile)

        # Save Warmed Topology
        with open(self.WarmedTopologyFile, "w") as pdb_file:
            app.PDBFile.writeFile(
                warmup_simulation.topology,
                warmup_simulation.context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(),
                file=pdb_file,
                keepIds=True,
            )

        # NVT Equilibration
        # Production Run
    def production(self):

        # Define Integrator for Simulation
        integrator = LangevinIntegrator(self.config['temperatureValue'] * unit.kelvin,
                                        self.config['frictionValue'] / unit.picoseconds,
                                        self.config['timeStepValue'] * unit.femtoseconds)
        # Load Model if not Loaded (?)
        if not self.complex_model:
            self.load_model()

        # Load System if not Loaded (?)
        if not self.system:
            self.load_system()

        # Initialize NonBonded Forces Reporter
        if self.config['NonBondedReporter']:
            for force in self.system.getForces():
                if isinstance(force, NonbondedForce):
                    force.setForceGroup(0)
                    force.addGlobalParameter("protein_coulomb_scale", 1)
                    force.addGlobalParameter("protein_lj_scale", 1)
                    force.addGlobalParameter("ligand_coulomb_scale", 1)
                    force.addGlobalParameter("ligand_lj_scale", 1)
                    force.addGlobalParameter("solvent_coulomb_scale", 1)
                    force.addGlobalParameter("solvent_lj_scale", 1)
                    for i in range(force.getNumParticles()):
                        charge, sigma, epsilon = force.getParticleParameters(i)

                        force.setParticleParameters(i, 0, 0, 0)
                        if i in self.protein_atoms:
                            force.addParticleParameterOffset("protein_coulomb_scale", i, charge, 0, 0)
                            force.addParticleParameterOffset("protein_lj_scale", i, 0, sigma, epsilon)
                        elif i in self.ligand_atoms:
                            force.addParticleParameterOffset("ligand_coulomb_scale", i, charge, 0, 0)
                            force.addParticleParameterOffset("ligand_lj_scale", i, 0, sigma, epsilon)
                        else:
                            force.addParticleParameterOffset("solvent_coulomb_scale", i, charge, 0, 0)
                            force.addParticleParameterOffset("solvent_lj_scale", i, 0, sigma, epsilon)
                    for i in range(force.getNumExceptions()):
                        p1, p2, chargeProd, sigma, epsilon = force.getExceptionParameters(i)
                        force.setExceptionParameters(i, p1, p2, 0, 0, 0)
                else:
                    force.setForceGroup(2)

        # Initialize Production Simulation
        production_simulation = app.Simulation(self.complex_model.topology, self.system, integrator, platform=self.platform)

        # Check if we have any checkpoint: Production or Warmed

        if os.path.exists(self.ProductionCheckpointFile):
            # Load Production Checkpoint if possible
            production_simulation.loadCheckpoint(self.ProductionCheckpointFile)
        elif os.path.exists(self.WarmedCheckpointFile):
            # Load Warmed Up Checkpoint if possible
            production_simulation.loadCheckpoint(self.WarmedCheckpointFile)

        # Set CheckPoint Reporter
        production_simulation.reporters.append(
            app.checkpointreporter.CheckpointReporter(
                file = self.ProductionCheckpointFile,
                reportInterval = self.config['productionWriteInterval'])
        )
        # Set Trajectory Reporter
        production_simulation.reporters.append(
            app.xtcreporter.XTCReporter(file = self.ProductionTrajectoryFile,reportInterval = self.config['productionWriteInterval'],enforcePeriodicBox = False,append = os.path.exists(self.ProductionTrajectoryFile))
        )
        # Set State Data Reporter
        production_simulation.reporters.append(
            app.StateDataReporter(
                file = self.ProductionStateDataReporterFile,
                reportInterval = self.config['productionWriteInterval'],
                step = True,
                potentialEnergy = True,
                kineticEnergy = True,
                temperature = True,
                volume = True,
                progress = True,
                remainingTime = True,
                speed = True,
                totalSteps = self.config['maxSteps'] + self.heating_steps,
                separator = "\t",
                append = os.path.exists(self.ProductionStateDataReporterFile),
            )
        )
        # Set Non-Bonded Forces Reporter such as L-J and vdW
        if self.config['NonBondedReporter']:
            production_simulation.reporters.append(
                ForceReporter(
                    file=self.WarmupInteractionForceDataReporterFile,
                    reportInterval=self.config['productionWriteInterval']
                )
            )

        state = production_simulation.context.getState(getVelocities=True, getPositions=True)
        positions = state.getPositions()
        velocities = state.getVelocities()
        production_simulation.context.setPositions(positions)
        production_simulation.context.setVelocities(velocities)
        production_simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
        production_simulation.step(self.steps_to_run)
