import os
import pandas as pd
import numpy as np
import re
from Bio import PDB
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline


PDB_NUCLEOTIDES_DICT = {
    '5CM':'C',
    'DC':'C',
    'DA':'A',
    'DT':'T',
    'DG':'G',
}

def identify_chain_type(seq):
    alphabets = {'dna': re.compile('^[acgtn]*$', re.I),
                 'protein': re.compile('^[acdefghiklmnpqrstvwy]*$', re.I)}

    if alphabets['dna'].search(str(seq)) is not None:
        return 'dna'
    elif alphabets['protein'].search(str(seq)) is not None:
        return 'protein'
    else:
        return None
def extract_chains(record):
    chain_ids = record.split('|')[1]
    chain_ids = chain_ids.replace('Chains ','').replace('Chain','')
    chain_ids = re.sub(r'([A-Za-z0-9]+)\[auth\s([A-Za-z0-9]+)\]', r'\2', chain_ids)
    chain_ids = re.sub(r'\s+', ' ', chain_ids).strip()
    chain_ids = re.sub(r'\s+', '', chain_ids)
    chain_ids = chain_ids.split(',')


    return chain_ids

class Assembly:
    def __init__(self,structure_filepath,
                 list_of_fasta_filepaths):

        self.assembly = os.path.splitext(os.path.basename(structure_filepath))[0].split('.')[0]
        self.entry = self.assembly.split('_')[0]
        self.structure_filepath = structure_filepath
        self.fasta_filepath = [fasta_filepath for fasta_filepath in list_of_fasta_filepaths if self.entry in fasta_filepath][0]
        self.chains_types = self.identify_chains_types()
        self.DNA_FASTA_sequence = self.get_dna_fasta_sequence()
        self.DNA_PDB_sequence = self.get_dna_pdb_sequence()
        self.methylation, self.methylation_indexes = self.is_methylated()
    def identify_chains_types(self):
        chains_dict = {}
        # PDB Parser
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("structure", self.structure_filepath)
        pdb_chains_ids = []
        for model in structure:
            for chain in model:
                pdb_chains_ids.append(chain.id)
        fasta_chains_ids = []
        fasta_chains_types = []
        for record in SeqIO.parse(self.fasta_filepath, "fasta"):
            fasta_chain_id = extract_chains(record.description)
            fasta_chain_type = identify_chain_type(record.seq)
            fasta_chains_ids.append(fasta_chain_id)
            fasta_chains_types.append(fasta_chain_type)

        for fasta_chains_id, fasta_chain_type in zip(fasta_chains_ids, fasta_chains_types):
            for fasta_id in fasta_chains_id:
                if fasta_id in pdb_chains_ids:
                    chains_dict[fasta_id] = fasta_chain_type
        return chains_dict
    def get_dna_fasta_sequence(self):
        seq_dict = {}
        for chain_id in self.chains_types:
            if self.chains_types[chain_id] == 'dna':
                for record in SeqIO.parse(self.fasta_filepath, "fasta"):
                    fasta_chain_ids = extract_chains(record.description)
                    if chain_id in fasta_chain_ids:
                        seq = str(record.seq)
                        seq_dict[chain_id] = seq
        return seq_dict
    def get_dna_pdb_sequence(self):
        seq_dict = {}
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("structure", self.structure_filepath)
        for chain_id in self.chains_types:
            if self.chains_types[chain_id] == 'dna':
                for model in structure:
                    for chain in model:
                        if chain.id == chain_id:
                            seq=''
                            for residue in chain:
                                seq += PDB_NUCLEOTIDES_DICT[residue.resname] if residue.resname in PDB_NUCLEOTIDES_DICT.keys() else ''
                            seq_dict[chain_id] = seq
        return seq_dict
    def get_assembly_dict(self):
        assembly_dict = {}
        assembly_dict["assembly"] = self.assembly
        assembly_dict["entry"] = self.entry
        assembly_dict["structure_filepath"] = self.structure_filepath
        assembly_dict["fasta_filepath"] = self.fasta_filepath
        assembly_dict["chains_types"] = self.chains_types
        assembly_dict['DNA_FASTA_sequence'] = self.DNA_FASTA_sequence
        assembly_dict['DNA_PDB_sequence'] = self.DNA_PDB_sequence
        assembly_dict['is_methylated'] = self.methylation
        assembly_dict['methylated_indexes'] = self.methylation_indexes
        return assembly_dict
    def is_methylated(self):
        meth_dict = {}
        meth_indexes_dict = {}
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("structure", self.structure_filepath)
        for chain_id in self.chains_types:
            if self.chains_types[chain_id] == 'dna':
                for model in structure:
                    for chain in model:
                        if chain.id == chain_id:
                            meth_indexes = []
                            is_meth = False
                            for res_index, residue in enumerate(chain):
                                if residue.resname == '5CM':
                                    meth_indexes.append(res_index+1)
                            if len(meth_indexes) > 0:
                                is_meth = True
                            meth_dict[chain_id] = is_meth
                            meth_indexes_dict[chain_id] = meth_indexes

        return meth_dict, meth_indexes_dict

def fasta_dna_sequence_similarity_analysis(df,WORKDIR):
    # Get all Strand One Sequences
    seqs = {}
    for _, row in df.iterrows():
        chain = next(iter(row['DNA_PDB_sequence']))
        seq = row['DNA_FASTA_sequence'][chain]
        seqs[row['assembly']+'_'+chain] = seq

    # Align Sequences
    temp_fasta_file_input = os.path.join(WORKDIR,'temp','temp_fasta_strand_1_dna_sequences.fasta')
    temp_fasta_file_output = os.path.join(WORKDIR,'temp','temp_aligned_fasta_strand_1_dna_sequences.fasta')
    clustalo
    with open(temp_fasta_file_input,'w') as f:
        for key, seq in seqs.items():
            f.write(f'>{key}\n{seq}\n')
    if
    clustalomega_cline = ClustalOmegaCommandline(infile=temp_fasta_file_input,
                                                 outfile=temp_fasta_file_output,
                                                 verbose=True,
                                                 auto=True)
    clustalomega_cline()
    alignment = AlignIO.read(temp_fasta_file_output, "fasta")
    for record in alignment:
        print(record.id, record.seq)

    a=1



def run(WORKDIR):
    list_of_assemblies = sorted(os.listdir(os.path.join(WORKDIR,'rawPDB')), key=lambda s: s[:4])
    list_of_pdb_filepaths = [os.path.join(WORKDIR,'rawPDB',x) for x in list_of_assemblies]
    list_of_fasta_filepaths = [os.path.join(WORKDIR,'rawSeq',x) for x in os.listdir(os.path.join(WORKDIR,'rawSeq'))]


    DATAFRAME_FILEPATH = os.path.join(WORKDIR,'dataframe.pkl')

    if not os.path.exists(DATAFRAME_FILEPATH):
        list_of_assemblies = []
        for structure_filepath in list_of_pdb_filepaths:
            assembly = Assembly(structure_filepath, list_of_fasta_filepaths)
            assembly_dict = assembly.get_assembly_dict()
            list_of_assemblies.append(assembly_dict)
        df = pd.DataFrame.from_dict(list_of_assemblies)
        df.to_pickle(DATAFRAME_FILEPATH)
    else:
        df = pd.read_pickle(DATAFRAME_FILEPATH)
    # Alignment i Podobieństwo Sekwencji DNA z FASTA
    fasta_dna_sequence_similarity_analysis(df,WORKDIR=WORKDIR)
    # Alignment i Podobieństwo Sekwencji DNA z PDB
    # Podobieństwo strukturalne łańcuchów CTCF
if __name__ == '__main__':
    WORKDIR = '/Users/maciejwisniewski/data/SachinCTCF'
    run(WORKDIR)
