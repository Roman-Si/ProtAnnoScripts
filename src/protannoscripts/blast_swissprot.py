import re
import gzip
import pandas as pd
from Bio import SeqIO


def parse_sprot_fasta(sprot_fasta_file: str) -> pd.DataFrame:
    """
    Give a SwissProt FASTA file and return a DataFrame with protein accession, protein name, gene name, and taxonomy.
    Latest version can be downloaded with "wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"

    Example usage:
    >>> df = parse_sprot_fasta('uniprot_sprot.fasta.gz')

    Since SwissProt contains million entries you can save the dataframe to a csv file for further usage.
    """

    proteinIDs = []
    protein_names = []
    genes = []
    taxIDs = []

    # Check whether the fasta file is gzipped
    opener = gzip.open if sprot_fasta_file.endswith(".gz") else open

    with opener(sprot_fasta_file, "rt") as fall:
        records = SeqIO.parse(fall, "fasta")
        for record in records:
            proteinIDs.append(
                record.id.split("|")[1] if "|" in record.id else record.id
            )

            description = record.description.replace(record.id, "")

            protein_names.append(re.match(r"(.*) OS\=.*", description).group(1))
            genes.append(
                re.match(r".* GN\=(.*) PE", description).group(1)
                if " GN=" in description
                else None
            )

            taxIDs.append(
                re.match(r".* OX\=([0-9]*) ", description).group(1)
                if " OX=" in description
                else None
            )

    df = pd.DataFrame(
        {
            "sprotId": proteinIDs,
            "Protein name": protein_names,
            "Gene": genes,
            "Taxon Id": taxIDs,
        }
    )
    return df


def blast_sprot(
    blast_output: str, sprot_df: pd.DataFrame, pident=70, qcovs=70, length_difference=20
) -> pd.DataFrame:
    """
    Steps:
    Parse the output of blastp to SwissProt.
    Filter for pident, qcovs and length difference. The default pident and qcovs values are based on 10.1080/21501203.2011.606851 for fungal genome annotation.
    Keep alignment with best bitscore for each query.
    Return Datraframe .

    Parameters:
    blast_output (str): Tabular output of blast to SwissProt in outfmt 6 format. Needs the columns qacc, sacc, pident, qlen,slen, qcovs, bitscore.
    sprot_df (dataframe): Dataframe with the annotations of Swissprot records, created with parse_sprot_fasta function.
    pident (float, optional): Threshold of percentage identity, default 70%
    qcovs (float, optional): Threshold of query coverage, default 70%
    length_difference (float, optional): Threshold for max allowed length difference of query and subject sequences, default 20%

    Example blast command for creating blast_output in command line:
    ```
    echo "qacc sacc pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore qcovs stitle" | sed -E 's/ /\\t/g' > $output
    $BLASTdir/$program -db $blast_sprot_database -query $query_fasta -evalue $eval -num_threads $threads -outfmt "6 qacc sacc pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore qcovs stitle" >> $output
    ```

    Example usage:
    >>> blast_result = blast_sprot('sprot_2022-05.blastp.gz', sprot_df, pident=70, qcovs=70, length_difference=20)


    Returns:
    pandas.DataFrame: DataFrame with best alignment and its annotations.
    """

    df = pd.read_csv(blast_output, sep="\t", dtype={"qacc": object})
    # parse uniprot accession
    df["sacc"] = df["sacc"].apply(lambda x: x.split("|")[1] if "|" in x else x)
    # calculate subject coverage
    df["scovs"] = ((df["send"] - df["sstart"]) / df["slen"]) * 100
    # filter
    df = df.loc[
        (df["pident"] > pident)
        & (df["qcovs"] > qcovs)
        & (abs(df["qlen"] - df["slen"]) / df["qlen"] < length_difference)
    ]
    # keep best hit based on bitscore
    df = (
        df.sort_values(by=["qacc", "bitscore"], ascending=[True, False])
        .groupby("qacc")
        .head(1)
    )
    # merge with sprot_df to get protein annotations
    df = pd.merge(df, sprot_df, left_on="sacc", right_on="sprotId", how="left")
    df = df[
        [
            "qacc",
            "sprotId",
            "Protein name",
            "Gene",
            "Taxon Id",
            "pident",
            "qcovs",
            "scovs",
            "evalue",
            "bitscore",
        ]
    ]
    df["proteinId"] = df["qacc"]
    return df
