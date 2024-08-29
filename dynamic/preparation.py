import os
import requests
from openbabel import pybel
from Bio import SeqIO, PDB
from Bio.Align import PairwiseAligner
import modeller
import modeller.automodel
from dynamic import constant


class ProteinPreparation():
    def __init__(self, config):
        self.homedir = os.getcwd()
        assert isinstance(config['pdbid'], str), "pdb_id variable must be string!"
        self.pdbid = config['pdbid']

        assert os.path.exists(config['outdir']), "Output Directory does not exist!"
        self.outdir = config['outdir']

        #############
        # RAW FILES #
        #############

        self.raw_protein_pdb_files = [filepath for filepath in config['proteinFiles'] if os.path.exists(filepath)]
        self.protein_pdb_filenames = [os.path.basename(x).split('.')[0] for x in self.raw_protein_pdb_files]
        assert len(self.raw_protein_pdb_files) != 0, "Protein PDB files do not exist."

        #######################
        # PREPROCESSING FILES #
        #######################

        # List of files with protein fasta sequences for MODELLER reconstruction
        self.rcsb_protein_fasta_files = [os.path.join(self.outdir, x + '.fasta') for x in self.protein_pdb_filenames]
        # List of files with protein alignment sequences for MODELLER reconstruction
        self.alignment_modeller_files = [os.path.join(self.outdir, x + '.ali') for x in self.protein_pdb_filenames]
        # List of files with protein structure from MODELLER reconstruction
        self.modeller_protein_pdb_files = [os.path.join(self.outdir, x + '_protein.modeller.B99990001.pdb') for x in self.protein_pdb_filenames]
        # xChimera Dock Prep Script File
        self.dockprep_file = os.path.join(self.outdir, 'dockprep.cxc')
        # List of files with protein structure after xChimera Dock Preparation
        self.prepared_protein_pdb_files = [os.path.join(self.outdir, x + '_protein.prepared.pdb') for x in self.protein_pdb_filenames]

    def fetch_rcsb_pdb(self):
        ''' Fetch RCSB PDB complex sequence FASTA file. '''
        for rcsb_protein_fasta_file in self.rcsb_protein_fasta_files:
            url = 'https://www.rcsb.org/fasta/entry/' + str(self.pdbid.upper())
            response = requests.get(url)
            if response.status_code == 200:
                with open(rcsb_protein_fasta_file, 'wb') as f:
                    f.write(response.content)
                print("FASTA File downloaded successfully.")
            else:
                print("An error occurred while downloading the file.")
                print(response.text)
    def find_most_similar_chain(self, raw_seq: str, full_seqs: str):
        '''
        A function that compares the raw sequence from a PDB file with the full
        sequence from the RCSB Protein Data Bank to find the most similar string.

        :param raw_seq: Raw Sequence from PDB file.
        :param full_seqs: Complete Sequence from RCSB Protein Data Bank.
        :return:
        '''
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 1
        aligner.mismatch_score = 0
        aligner.open_gap_score = -1
        aligner.extend_gap_score = -0.5
        best_smiliarity = None
        best_sequence = None
        for i, full_seq in enumerate(full_seqs):
            print(i)
            print('full')
            print(full_seq)
            print('raw')
            print(raw_seq)
            alignment = aligner.align(raw_seq, full_seq)
            score = alignment.score
            if best_smiliarity is None or score > best_smiliarity:
                best_smiliarity = score
                best_sequence = full_seq
        return best_sequence, best_smiliarity
    def modeller_run(self):
        '''
        Run MODELLER for protein structure reconstruction.
        Remember to export MODELLER Key!
        Steps:
            * Compare Full Sequence from RSCB Protein Data Bank with our PDB Structure with missing residues.
            * Align both sequences.
            * Rebuild missing residues in raw structure.
        '''
        # Go to the Output Directory
        os.chdir(self.outdir)
        # Create Full Sequence Structure via MODELLER Loop
        for raw_protein_pdb_file, rcsb_protein_fasta_file, alignment_modeller_file, protein_pdb_filename in zip(self.raw_protein_pdb_files,self.rcsb_protein_fasta_files, self.alignment_modeller_files,self.protein_pdb_filenames):
            structure_index = self.protein_pdb_filenames.index(protein_pdb_filename)

            # Load Structure from Raw PDB File
            parser = PDB.PDBParser()
            structure = parser.get_structure(self.pdbid, file=raw_protein_pdb_file)

            # Read Sequence from RCSB
            records = list(SeqIO.parse(rcsb_protein_fasta_file, "fasta"))
            full_sequences = [str(x.seq) for x in records]

            # Read Sequence from Raw Protein PDB File
            structure_sequences = []
            for model in structure:
                no_protein_res = 0
                for chain in model:
                    chain_seq = []
                    for residue in chain:
                        try:
                            chain_seq.append(constant.aa_to_aa_id[residue.resname])
                        except:
                            no_protein_res = no_protein_res + 1
                    structure_sequences.append(chain_seq)
            structure_sequences = [''.join(x) for x in structure_sequences]
            structure_sequences = list(filter(None, structure_sequences))

            # Find Matching Sequences for Raw Sequences
            chosen_sequences = []
            for structure_sequence in structure_sequences:
                best_sequence, best_smiliarity = self.find_most_similar_chain(structure_sequence, full_sequences)
                chosen_sequences.append(best_sequence)

            # Fill Raw Sequences with dots if Non-Protein Residues are present in Full Sequence
            structure_sequences = '/'.join(structure_sequences) + no_protein_res * '.'
            chosen_sequences = '/'.join(chosen_sequences) + no_protein_res * '.'

            # Alignment of FASTA sequences
            print('\t[**]Modeller - Alignment with FASTA sequence')

            env = modeller.Environ()        # Define Environment class
            env.io.hetatm = True            # Allow to add HETATM into Output File
            aln = modeller.Alignment(env)   # Define Alignment class
            mdl = modeller.Model(env)       # Define Model class

            # Add Raw Protein PDB File to our Model
            mdl.read(file=raw_protein_pdb_file)

            #############
            # ALIGNMENT #
            #############

            # Add Raw Protein PDB File to our Aligner
            aln.append_model(mdl, align_codes=raw_protein_pdb_file, atom_files=raw_protein_pdb_file)

            # Add Full Sequence to add missing residues in Raw PDB Structure
            aln.append_sequence(chosen_sequences)

            # Rename Alignment instances
            aln[1].code = self.pdbid + '_seq'

            # Define alignment penalties
            aln.malign(gap_penalties_1d=(-500, -300))

            # Save Alignment file
            aln.write(file=alignment_modeller_file)

            #############
            # MODELLING #
            #############

            try:
                # Try to model missing loops in
                modelling = modeller.automodel.LoopModel(env, alnfile=alignment_modeller_file,
                                                         knowns=raw_protein_pdb_file,
                                                         sequence=self.pdbid + '_seq',
                                                         root_name=protein_pdb_filename + '_protein.modeller')
                modelling.starting_model = 1
                modelling.ending_model = 1
                modelling.loop.starting_model = 1
                modelling.loop.ending_model = 2
                modelling.loop.md_level = modeller.automodel.refine.fast
                modelling.make()

            except:
                if not os.path.exists(self.modeller_protein_pdb_files[structure_index]):
                    self.modeller_protein_pdb_files[structure_index] = raw_protein_pdb_file

        # Back to Home Directory
        os.chdir(self.homedir)
    def save_perpared_file(self):
        for modeller_protein_pdb_file, prepared_protein_pdb_file in zip(self.modeller_protein_pdb_files,
                                                                        self.prepared_protein_pdb_files):
            with open(modeller_protein_pdb_file, 'r') as infile, open(prepared_protein_pdb_file, 'w') as outfile:
                # Read the content of the input file
                content = infile.read()
                # Write the content to the output file
                outfile.write(content)
    def save_unprepared_file(self):
        for raw_protein_pdb_file, prepared_protein_pdb_file in zip(self.raw_protein_pdb_files,
                                                                   self.prepared_protein_pdb_files):
            with open(raw_protein_pdb_file, 'r') as infile, open(prepared_protein_pdb_file, 'w') as outfile:
                # Read the content of the input file
                content = infile.read()
                # Write the content to the output file
                outfile.write(content)

class LigandPreparation():
    def __init__(self, config):
        self.raw_ligand_pdb_files = config['ligandFiles']
        self.new_sdf_ligand_files = [filepath.replace('pdb','sdf') for filepath in config['ligandFiles']]
        self.new_mol2_ligand_files = [filepath.replace('pdb','mol2') for filepath in config['ligandFiles']]
    def pdb_to_mol2(self):
        for raw_file, mol2_file in zip(self.raw_ligand_pdb_files, self.new_mol2_ligand_files):
            pybel_mol = next(pybel.readfile('pdb', raw_file))
            pybel_mol.write('mol2',mol2_file,overwrite=True)
    def pdb_to_sdf(self):
        for raw_file, sdf_file in zip(self.raw_ligand_pdb_files, self.new_sdf_ligand_files):
            pybel_mol = next(pybel.readfile('pdb', raw_file))
            pybel_mol.write('sdf',sdf_file,overwrite=True)