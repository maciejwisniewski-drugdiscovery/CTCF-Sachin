import os
import argparse
import tempfile

import pandas as pd
import numpy as np
import re
from Bio import PDB
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
import seaborn as sns
sns.set_theme()
import matplotlib.pyplot as plt



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
        self.PROTEIN_FASTA_sequence = self.get_protein_fasta_sequence()
        self.PROTEIN_PDB_sequence = self.get_protein_pdb_sequence()
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
                                if residue.resname.strip() in list(PDB_NUCLEOTIDES_DICT.keys()):
                                    seq += PDB_NUCLEOTIDES_DICT[residue.resname.strip()]
                            seq_dict[chain_id] = seq
        return seq_dict
    def get_protein_fasta_sequence(self):
        seq_dict = {}
        for chain_id in self.chains_types:
            if self.chains_types[chain_id] == 'protein':
                for record in SeqIO.parse(self.fasta_filepath, "fasta"):
                    fasta_chain_ids = extract_chains(record.description)
                    if chain_id in fasta_chain_ids:
                        seq = str(record.seq)
                        seq_dict[chain_id] = seq
        return seq_dict
    def get_protein_pdb_sequence(self):
        seq_dict = {}
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("structure", self.structure_filepath)
        for chain_id in self.chains_types:
            if self.chains_types[chain_id] == 'protein':
                for model in structure:
                    for chain in model:
                        seq_records = []
                        if chain.id == chain_id:
                            ppb = PDB.PPBuilder()
                            peptides = ppb.build_peptides(chain)
                            for peptide in peptides:
                                seq_records.append(peptide.get_sequence())
                            seq_dict[chain_id] = str(seq_records[0])
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
        assembly_dict['PROTEIN_FASTA_sequence'] = self.PROTEIN_FASTA_sequence
        assembly_dict['PROTEIN_PDB_sequence'] = self.PROTEIN_PDB_sequence
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

    temp_fasta_file_input = os.path.join(WORKDIR,'temp','temp_fasta_strand_1_dna_sequences.fasta')
    temp_fasta_file_output = os.path.join(WORKDIR,'temp','temp_aligned_fasta_strand_1_dna_sequences.aln')
    similarity_matrix_filepath = os.path.join(WORKDIR,'similarity','fasta_dna_sequence_clustalo_similairty_matrix.npy')
    similarity_heatmap_filepath = os.path.join(WORKDIR,'similarity','fasta_dna_sequence_clustalo_similairty_heatmap.png')


    if not os.path.exists(similarity_matrix_filepath):
        # Get all Strand One Sequences
        seqs = {}

        for _, row in df.iterrows():
            chain = next(iter(row['DNA_PDB_sequence']))
            seq = row['DNA_FASTA_sequence'][chain]
            seqs[row['assembly']+'_'+chain] = seq

        with open(temp_fasta_file_input,'w') as f:
            for key, seq in seqs.items():
                f.write(f'>{key}\n{seq}\n')

        clustalomega_cline = ClustalOmegaCommandline(infile=temp_fasta_file_input,
                                                     outfile=temp_fasta_file_output,
                                                     verbose=True,
                                                     auto=True)
        clustalomega_cline()
        alignment = AlignIO.read(temp_fasta_file_output, "fasta")

        num_sequences = len(alignment)
        similarity_matrix = np.zeros((num_sequences, num_sequences))
        #   Obliczanie macierzy podobieństwa
        for i in range(num_sequences):
            for j in range(num_sequences):
                matches = sum(res1 == res2 for res1, res2 in zip(alignment[i], alignment[j]))
                similarity_matrix[i, j] = matches / len(alignment[0])

        np.save(similarity_matrix_filepath, similarity_matrix)
    else:
        similarity_matrix = np.load(similarity_matrix_filepath)

    import matplotlib.pyplot as plt

    # Utwórz mapowanie kolorów dla każdej grupy w kolumnie 'Entry'
    unique_entries = df['entry'].unique()
    colors = sns.color_palette("dark", len(unique_entries))
    entry_color_map = dict(zip(unique_entries, colors))

    # Stwórz listę kolorów odpowiadających każdemu assembly
    assembly_colors = df['entry'].map(entry_color_map)

    # Stwórz heatmapę z dodaniem kolorowych pasków (color bars)
    plt.figure(figsize=(12, 10))

    # Dodaj heatmapę
    sns.heatmap(similarity_matrix, annot=False, cmap="flare",
                xticklabels=df['assembly'], yticklabels=df['assembly'],
                cbar_kws={'label': 'Similarity'}, linewidths=.5, linecolor='lightgrey')

    # Odległość etykiet od heatmapy
    plt.gca().tick_params(axis='x', labelsize=12, pad=10)  # Większy pad (odstęp) na osi X
    plt.gca().tick_params(axis='y', labelsize=12, pad=20)  # Większy pad (odstęp) na osi Y

    # Dodaj paski kolorów na zewnątrz heatmapy
    for (i, color) in enumerate(assembly_colors):
        # plt.gca().add_patch(plt.Rectangle((i, -0.5), 1, 0.5, color=color, transform=plt.gca().transData, clip_on=True))
        plt.gca().add_patch(
            plt.Rectangle((-1, i + 0.05), 0.5, 0.9, color=color, transform=plt.gca().transData, clip_on=False))

    # Dodaj tytuł
    plt.text(x=-0.75, y=-0.5, s='Entries', fontsize=10, ha='center', rotation='vertical')
    plt.text(x=15, y=-0.5, s='HeatMap', fontsize=10, ha='center')
    plt.title(label="Similarity HeatMap of Fasta DNA Sequences with Colored Labels for each Entry", fontsize=12, loc='center', x=0.5, y=1.05)
    plt.ylabel('Assembly')

    plt.savefig(similarity_heatmap_filepath)
def pdb_dna_sequence_similarity_analysis(df,WORKDIR):

    temp_fasta_file_input = os.path.join(WORKDIR,'temp','temp_pdb_strand_1_dna_sequences.fasta')
    temp_fasta_file_output = os.path.join(WORKDIR,'temp','temp_aligned_pdb_strand_1_dna_sequences.aln')
    similarity_matrix_filepath = os.path.join(WORKDIR,'similarity','pdb_dna_sequence_clustalo_similairty_matrix.npy')
    similarity_heatmap_filepath = os.path.join(WORKDIR,'similarity','pdb_dna_sequence_clustalo_similairty_heatmap.png')

    if not os.path.exists(similarity_matrix_filepath):
        # Get all Strand One Sequences
        seqs = {}

        for _, row in df.iterrows():
            chain = next(iter(row['DNA_PDB_sequence']))
            seq = row['DNA_PDB_sequence'][chain]
            seqs[row['assembly']+'_'+chain] = seq

        with open(temp_fasta_file_input,'w') as f:
            for key, seq in seqs.items():
                f.write(f'>{key}\n{seq}\n')

        clustalomega_cline = ClustalOmegaCommandline(infile=temp_fasta_file_input,
                                                     outfile=temp_fasta_file_output,
                                                     verbose=True,
                                                     auto=True)
        clustalomega_cline()
        alignment = AlignIO.read(temp_fasta_file_output, "fasta")

        num_sequences = len(alignment)
        similarity_matrix = np.zeros((num_sequences, num_sequences))
        #   Obliczanie macierzy podobieństwa
        for i in range(num_sequences):
            for j in range(num_sequences):
                matches = sum(res1 == res2 for res1, res2 in zip(alignment[i], alignment[j]))
                similarity_matrix[i, j] = matches / len(alignment[0])

        np.save(similarity_matrix_filepath, similarity_matrix)
    else:
        similarity_matrix = np.load(similarity_matrix_filepath)

    import matplotlib.pyplot as plt

    # Utwórz mapowanie kolorów dla każdej grupy w kolumnie 'Entry'
    unique_entries = df['entry'].unique()
    colors = sns.color_palette("dark", len(unique_entries))
    entry_color_map = dict(zip(unique_entries, colors))

    # Stwórz listę kolorów odpowiadających każdemu assembly
    assembly_colors = df['entry'].map(entry_color_map)

    # Stwórz heatmapę z dodaniem kolorowych pasków (color bars)
    plt.figure(figsize=(12, 10))

    # Dodaj heatmapę
    sns.heatmap(similarity_matrix, annot=False, cmap="flare",
                xticklabels=df['assembly'], yticklabels=df['assembly'],
                cbar_kws={'label': 'Similarity'}, linewidths=.5, linecolor='lightgrey')

    # Odległość etykiet od heatmapy
    plt.gca().tick_params(axis='x', labelsize=12, pad=10)  # Większy pad (odstęp) na osi X
    plt.gca().tick_params(axis='y', labelsize=12, pad=20)  # Większy pad (odstęp) na osi Y

    # Dodaj paski kolorów na zewnątrz heatmapy
    for (i, color) in enumerate(assembly_colors):
        # plt.gca().add_patch(plt.Rectangle((i, -0.5), 1, 0.5, color=color, transform=plt.gca().transData, clip_on=True))
        plt.gca().add_patch(
            plt.Rectangle((-1, i + 0.05), 0.5, 0.9, color=color, transform=plt.gca().transData, clip_on=False))

    # Dodaj tytuł
    plt.text(x=-0.75, y=-0.5, s='Entries', fontsize=10, ha='center', rotation='vertical')
    plt.text(x=15, y=-0.5, s='HeatMap', fontsize=10, ha='center')
    plt.title(label="Similarity HeatMap of PDB Crystal DNA Sequences with Colored Labels for each Entry", fontsize=12, loc='center', x=0.5, y=1.05)
    plt.ylabel('Assembly')

    plt.savefig(similarity_heatmap_filepath)
def fasta_protein_sequence_similarity_analysis(df,WORKDIR):

    # Specify INPUT, TEMP, OUTPUT FILES
    temp_fasta_file_input = os.path.join(WORKDIR,'temp','temp_fasta_protein_sequences.fasta')
    temp_fasta_file_output = os.path.join(WORKDIR,'temp','temp_aligned_fasta_protein_sequences.aln')
    similarity_matrix_filepath = os.path.join(WORKDIR,'similarity','fasta_protein_sequence_clustalo_similairty_matrix.npy')
    similarity_heatmap_filepath = os.path.join(WORKDIR,'similarity','fasta_protein_sequence_clustalo_similairty_heatmap.png')

    # CREATE SIMILARITY MATRIX FOR ALL CHAINS !
    if not os.path.exists(similarity_matrix_filepath):
        # List for all protein chain structures
        all_chains_sequences = []
        # Iterate over our DataFrame
        for _, row in df.iterrows():
            # Iterate over Protein Chains in Specific Assembly
            for chain_id, chain_sequence in row['PROTEIN_FASTA_sequence'].items():
                # Create Protein Chain Dict
                seqdict = {}
                seqdict['structure_id'] = row['assembly']+'_'+chain_id
                seqdict['assembly'] = row['assembly']
                seqdict['entry'] = row['entry']
                seqdict['protein_sequence'] = chain_sequence
                all_chains_sequences.append(seqdict)

        with open(temp_fasta_file_input,'w') as f:
            # Create FASTA file with all proteins
            for sequence in all_chains_sequences:
                f.write(f'>{sequence["structure_id"]}\n{sequence["protein_sequence"]}\n')

        # CLustal Omega Alignment
        clustalomega_cline = ClustalOmegaCommandline(infile=temp_fasta_file_input,
                                                     outfile=temp_fasta_file_output,
                                                     verbose=True,
                                                     auto=True)
        clustalomega_cline()
        alignment = AlignIO.read(temp_fasta_file_output, "fasta")

        num_sequences = len(alignment)
        similarity_matrix = np.zeros((num_sequences, num_sequences))
        #   Obliczanie macierzy podobieństwa
        for i in range(num_sequences):
            for j in range(num_sequences):
                matches = sum(res1 == res2 for res1, res2 in zip(alignment[i], alignment[j]))
                similarity_matrix[i, j] = matches / len(alignment[0])

        np.save(similarity_matrix_filepath, similarity_matrix)
    else:
        similarity_matrix = np.load(similarity_matrix_filepath)

    new_df = pd.DataFrame(all_chains_sequences)

    # Utwórz mapowanie kolorów dla każdej grupy w kolumnie 'Entry' i 'Assembly'
    unique_entries = new_df['entry'].unique()
    unique_assemblies = new_df['assembly'].unique()

    entries_colors = sns.color_palette("dark", len(unique_entries))
    assemblies_colors = sns.color_palette("dark", len(unique_entries))

    entry_color_map = dict(zip(unique_entries, entries_colors))
    assembly_color_map = dict(zip(unique_assemblies, assemblies_colors))

    # Stwórz listę kolorów odpowiadających każdemu assembly i entry
    entry_colors = new_df['entry'].map(entry_color_map)
    assembly_colors = new_df['assembly'].map(assembly_color_map)

    # Stwórz heatmapę z dodaniem kolorowych pasków (color bars)
    plt.figure(figsize=(12, 10))

    # Dodaj heatmapę
    sns.heatmap(similarity_matrix, annot=False, cmap="flare",
                xticklabels=new_df['structure_id'], yticklabels=new_df['structure_id'],
                cbar_kws={'label': 'Similarity'}, linewidths=.5, linecolor='lightgrey')

    # Odległość etykiet od heatmapy
    plt.gca().tick_params(axis='x', labelsize=12, pad=10)  # Większy pad (odstęp) na osi X
    plt.gca().tick_params(axis='y', labelsize=12, pad=20)  # Większy pad (odstęp) na osi Y

    # Dodaj paski kolorów na zewnątrz heatmapy
    for (i, color) in enumerate(assembly_colors):
        # plt.gca().add_patch(plt.Rectangle((i, -0.5), 1, 0.5, color=color, transform=plt.gca().transData, clip_on=True))
        plt.gca().add_patch(
            plt.Rectangle((-1, i + 0.05), 0.5, 0.9, color=color, transform=plt.gca().transData, clip_on=False))
    for (i, color) in enumerate(entry_colors):
        # plt.gca().add_patch(plt.Rectangle((i, -0.5), 1, 0.5, color=color, transform=plt.gca().transData, clip_on=True))
        plt.gca().add_patch(
            plt.Rectangle((-1.5, i + 0.05), 0.5, 0.9, color=color, transform=plt.gca().transData, clip_on=False))
    # Dodaj tytuł
    plt.text(x=-0.75, y=-0.5, s='Assemblies', fontsize=10, ha='center', rotation='vertical')
    plt.text(x=-1.25, y=-0.5, s='Entries', fontsize=10, ha='center', rotation='vertical')

    plt.text(x=15, y=-0.5, s='HeatMap', fontsize=10, ha='center')
    plt.title(label="Similarity HeatMap of Fasta Protein Sequences with Colored Labels for each Entry", fontsize=12, loc='center', x=0.5, y=1.05)
    plt.ylabel('Structures [Entry_Assembly_Chain]')
    plt.savefig(similarity_heatmap_filepath)
def pdb_protein_sequence_similarity_analysis(df,WORKDIR):

    temp_fasta_file_input = os.path.join(WORKDIR,'temp','temp_pdb_protein_sequences.fasta')
    temp_fasta_file_output = os.path.join(WORKDIR,'temp','temp_aligned_pdb_protein_sequences.aln')
    similarity_matrix_filepath = os.path.join(WORKDIR,'similarity','pdb_protein_sequence_clustalo_similairty_matrix.npy')
    similarity_heatmap_filepath = os.path.join(WORKDIR,'similarity','pdb_protein_sequence_clustalo_similairty_heatmap.png')

    # CREATE SIMILARITY MATRIX FOR ALL CHAINS !
    if not os.path.exists(similarity_matrix_filepath):
        all_chains_sequences = []
        for _, row in df.iterrows():
            for chain_id, chain_sequence in row['PROTEIN_PDB_sequence'].items():
                seqdict = {}
                seqdict['structure_id'] = row['assembly']+'_'+chain_id
                seqdict['assembly'] = row['assembly']
                seqdict['entry'] = row['entry']
                seqdict['protein_sequence'] = chain_sequence
                all_chains_sequences.append(seqdict)

        with open(temp_fasta_file_input,'w') as f:
            for sequence in all_chains_sequences:
                f.write(f'>{sequence["structure_id"]}\n{sequence["protein_sequence"]}\n')

        clustalomega_cline = ClustalOmegaCommandline(infile=temp_fasta_file_input,
                                                     outfile=temp_fasta_file_output,
                                                     verbose=True,
                                                     auto=True)
        clustalomega_cline()
        alignment = AlignIO.read(temp_fasta_file_output, "fasta")

        num_sequences = len(alignment)
        similarity_matrix = np.zeros((num_sequences, num_sequences))
        #   Obliczanie macierzy podobieństwa
        for i in range(num_sequences):
            for j in range(num_sequences):
                matches = sum(res1 == res2 for res1, res2 in zip(alignment[i], alignment[j]))
                similarity_matrix[i, j] = matches / len(alignment[0])

        np.save(similarity_matrix_filepath, similarity_matrix)
    else:
        similarity_matrix = np.load(similarity_matrix_filepath)

    new_df = pd.DataFrame(all_chains_sequences)

    # Utwórz mapowanie kolorów dla każdej grupy w kolumnie 'Entry' i 'Assembly'
    unique_entries = new_df['entry'].unique()
    unique_assemblies = new_df['assembly'].unique()

    entries_colors = sns.color_palette("dark", len(unique_entries))
    assemblies_colors = sns.color_palette("dark", len(unique_entries))

    entry_color_map = dict(zip(unique_entries, entries_colors))
    assembly_color_map = dict(zip(unique_assemblies, assemblies_colors))

    # Stwórz listę kolorów odpowiadających każdemu assembly i entry
    entry_colors = new_df['entry'].map(entry_color_map)
    assembly_colors = new_df['assembly'].map(assembly_color_map)

    # Stwórz heatmapę z dodaniem kolorowych pasków (color bars)
    plt.figure(figsize=(12, 10))

    # Dodaj heatmapę
    sns.heatmap(similarity_matrix, annot=False, cmap="flare",
                xticklabels=new_df['structure_id'], yticklabels=new_df['structure_id'],
                cbar_kws={'label': 'Similarity'}, linewidths=.5, linecolor='lightgrey')

    # Odległość etykiet od heatmapy
    plt.gca().tick_params(axis='x', labelsize=12, pad=10)  # Większy pad (odstęp) na osi X
    plt.gca().tick_params(axis='y', labelsize=12, pad=20)  # Większy pad (odstęp) na osi Y

    # Dodaj paski kolorów na zewnątrz heatmapy
    for (i, color) in enumerate(assembly_colors):
        # plt.gca().add_patch(plt.Rectangle((i, -0.5), 1, 0.5, color=color, transform=plt.gca().transData, clip_on=True))
        plt.gca().add_patch(
            plt.Rectangle((-1, i + 0.05), 0.5, 0.9, color=color, transform=plt.gca().transData, clip_on=False))
    for (i, color) in enumerate(entry_colors):
        # plt.gca().add_patch(plt.Rectangle((i, -0.5), 1, 0.5, color=color, transform=plt.gca().transData, clip_on=True))
        plt.gca().add_patch(
            plt.Rectangle((-1.5, i + 0.05), 0.5, 0.9, color=color, transform=plt.gca().transData, clip_on=False))
    # Dodaj tytuł
    plt.text(x=-0.75, y=-0.5, s='Assemblies', fontsize=10, ha='center', rotation='vertical')
    plt.text(x=-1.25, y=-0.5, s='Entries', fontsize=10, ha='center', rotation='vertical')

    plt.text(x=15, y=-0.5, s='HeatMap', fontsize=10, ha='center')
    plt.title(label="Similarity HeatMap of PDB Protein Sequences with Colored Labels for each Entry", fontsize=12,
              loc='center', x=0.5, y=1.05)
    plt.ylabel('Structures [Entry_Assembly_Chain]')
    plt.savefig(similarity_heatmap_filepath)
def pdb_protein_tmscore_analysis(df,WORKDIR):
    tmscore_matrix_filepath = os.path.join(WORKDIR, 'similarity',
                                              'pdb_protein_structure_tmscore_matrix.npy')
    tmscore_heatmap_filepath = os.path.join(WORKDIR, 'similarity',
                                               'pdb_protein_structure_tmscore_heatmap.png')
    protein_chains_folder = os.path.join(WORKDIR, 'temp','protein_chains')
    a=1
    structures = []
    # Generate list of structures
    for _, assembly in df.iterrows():
        for chain_id, chain_type in assembly['chains_types'].items():
            if chain_type == 'protein':
                assembly_dict = {}
                assembly_dict['structure_id'] = assembly['assembly']+'_'+chain_id
                assembly_dict['assembly'] = assembly['assembly']
                assembly_dict['entry'] = assembly['entry']
                assembly_dict['chain_id'] = chain_id
                pdb_chain_filepath = os.path.join(protein_chains_folder,assembly_dict['structure_id'])
                parser = PDB.PDBParser(QUIET=True)

                structure = parser.get_structure(assembly_dict['structure_id'],assembly['structure_filepath'])
                for model in structure:
                    for chain in model:
                        if chain.id == chain_id:
                            new_structure = PDB.Structure.Structure('new_structure')
                            new_model = PDB.Model.Model(0)
                            new_model.add(chain.copy())  # Skopiuj łańcuch do nowego modelu
                            new_structure.add(new_model)
                            with open(pdb_chain_filepath,'w') as f:
                                io = PDB.PDBIO()
                                io.set_structure(new_structure)
                                # Zapisz nowy plik PDB
                                io.save(pdb_chain_filepath)
                            if os.path.exists(pdb_chain_filepath):
                                assembly_dict['pdb_structure_filepath'] = pdb_chain_filepath
                structures.append(assembly_dict)






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

    os.makedirs(os.path.join(WORKDIR,'similarity'), exist_ok=True)
    os.makedirs(os.path.join(WORKDIR, 'temp'), exist_ok=True)
    os.makedirs(os.path.join(WORKDIR, 'temp','protein_chains'), exist_ok=True)
    # Alignment i Podobieństwo Sekwencji DNA z FASTA
    fasta_dna_sequence_similarity_analysis(df,WORKDIR=WORKDIR)
    # Alignment i Podobieństwo Sekwencji DNA z PDB
    pdb_dna_sequence_similarity_analysis(df,WORKDIR=WORKDIR)
    # Alignment i Podobieństwo Sekwencji Białka z FASTA
    fasta_protein_sequence_similarity_analysis(df,WORKDIR=WORKDIR)
    # Alignment i Podobieństwo Sekwencji DNA z PDB
    pdb_protein_sequence_similarity_analysis(df,WORKDIR=WORKDIR)
    # Podobieństwo strukturalne TM-Score łańcuchów CTCF
    #pdb_protein_tmscore_analysis(df,WORKDIR=WORKDIR)
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--workdir', type=str, default='/Users/maciejwisniewski/data/SachinCTCF',
                        help='Path to the working directory (default: current directory)')
    args = parser.parse_args()

    run(args.workdir)
