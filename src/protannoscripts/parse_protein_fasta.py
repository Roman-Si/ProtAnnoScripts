import gzip
import pandas as pd
from Bio import SeqIO
from pyopenms import ProteaseDigestion, AASequence

def parse_proteome_to_df(protein_fasta, mRNA_prefix) -> pd.DataFrame:
    """
    Take a protein FASTA and the mRNA prefix (useful to identify isoforms) and create a dataframe with one row per protein.  
    The dataframe has number of theoretical tryptic peptides as well (useful for proteomics).


    Parameters:
    protein_fasta (str): Path to protein FASTA, can be gzipped as well.
    mRNA_prefix (str): Prefix separating geneID from mRNA identifier. For BRAKER is ".t", for MAKER usually "-R".

    Usage:
    df = parse_proteome_to_df("../proteins.fasta", "-R")

    Returns:
    pandas.DataFrame: Dataframe with gebeId, proteinId, protein_length, theoretical_tryptic_peptides
    """

    data = {
        'proteinId': [],
        'protein_seq': [],
        'protein_length': [],
        'theoretical_tryptic_peptides': [],
    }

    # Open function based on file extension
    open_func = gzip.open if protein_fasta.endswith(".gz") else open

    with open_func(protein_fasta, "rt") as fall:
        records = SeqIO.parse(fall, 'fasta')
        for record in records:
            data['proteinId'].append(record.id)
            data['protein_seq'].append(str(record.seq))
            data['protein_length'].append(len(str(record.seq)))
            # Theoretical trypsin peptides
            result = list()
            dig = ProteaseDigestion()
            dig.setEnzyme("Trypsin/P")
            prot_sequence = AASequence.fromString(str(record.seq))
            dig.digest(prot_sequence, result, 6, 40)
            data['theoretical_tryptic_peptides'].append(len(set(result)))
    
    df = pd.DataFrame(data)
    df["geneId"] = df['proteinId'].str.split(mRNA_prefix).str[0]

    return df
