import os
import argparse
from rcsbsearchapi.search import AttributeQuery
import urllib.request

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--workdir', type=str, default=os.getcwd(),
                        help='Path to the working directory (default: current directory)')
    args = parser.parse_args()

    WORKDIR = args.workdir
    PDB_STRUCTURES_DIRECTORY = os.path.join(WORKDIR, 'rawPDB')
    FASTA_SEQUENCE_DIRECTORY = os.path.join(WORKDIR, 'rawSeq')

    os.makedirs(WORKDIR, exist_ok=True)
    os.makedirs(PDB_STRUCTURES_DIRECTORY, exist_ok=True)
    os.makedirs(FASTA_SEQUENCE_DIRECTORY, exist_ok=True)

    uniprot_query = AttributeQuery(
        attribute="rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
        operator="exact_match", value="P49711")
    organism_query = AttributeQuery(attribute="rcsb_entity_source_organism.taxonomy_lineage.name",
                                    operator="exact_match", value="Homo sapiens")
    dna_query = AttributeQuery(attribute="entity_poly.rcsb_entity_polymer_type", operator="exact_match", value="DNA")

    query = uniprot_query & organism_query & dna_query

    pdb_basic_url = 'https://files.rcsb.org/download/'
    fasta_basic_url = 'https://www.rcsb.org/fasta/entry/'

    for assembly in query("assembly"):
        print(assembly)
        assembly = assembly.split('-')
        assembly_pdb_url = pdb_basic_url + assembly[0] + '.pdb' + assembly[1]
        try:
            urllib.request.urlretrieve(assembly_pdb_url,
                                       os.path.join(PDB_STRUCTURES_DIRECTORY, assembly[0] + '_' + assembly[1] + '.pdb'))
        except:
            print('PDB problem')

    for entry in query("entry"):
        print(entry)
        assembly_fasta_url = fasta_basic_url + entry + '/download'
        try:
            urllib.request.urlretrieve(assembly_fasta_url, os.path.join(FASTA_SEQUENCE_DIRECTORY, entry + '.fasta'))
        except:
            print('Seq problem')
