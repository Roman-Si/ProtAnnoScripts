from collections import defaultdict
import gzip
import pandas as pd
from Bio import SearchIO


def annotate_transporters(tcdb_blast_path, tcdb_fams_path, tcdb_subs_path, pident = 30, qcovs=70, length_difference=20, hmm_bitscore = 100, tcdb_hmm_path=None) -> pd.DataFrame:
    """
    Annotate transporters using TCDB database.

    Steps
    1.  Parse the output of blastp to TCDB. 
        Filter for pident, qcovs and scovs. The default pident is 30% and coverage 70%. Keep alignment with best bitscore for each query.
    2.  Parse tblout output of hmmsearch to TCDB (optional)
        Filter for bitscore 100.
    3.  Merge annotations and prefer BLAST. Take family description and substrate annotation if available from the files tcdb_fams and tcdb_subs downloaded from TCDB. 
        They are available in data/

    Parameters:
    tcdb_blast_path (str): Tabular output of blast to TCDB in outfmt 6 format. Needs the columns qacc, sacc, pident, qlen,slen, qcovs, bitscore.
    tcdb_hmm (str, optional): Output of hmmsearch to TCDB hmm profiles in --tblout format.
    tcdb_fams_path (str): Path to file with tcdb familie descriptions. Available in data/tcdb_fams.tab
    tcdb_subs_path (str): Path to file with tcdb substrate annotations. Available in data/tcdb_subs.tab
    pident (float, optional): Threshold of BLAST percentage identity, default 30%
    qcovs (float, optional): Threshold of qcovs and scovs coverage, default 70%
    length_difference (float, optional): Threshold for max allowed length difference of query and subject sequences, default 20%    
    hmm_bitscore (float, optional): Threshold for hmm bitscore, default 100
        
    Usage:
    >>> df = annotate_transporters('tcdb.blastp.gz', "path/to/tcdb_fams.tab", "path/to/tcdb_substrates.tab", tcdb_fams_path = 'tcdb.hmm.gz')
    Returns:
    pandas.DataFrame: Dataframe with transporters and TCDB annotations.
    """
    
    ### 1. Read and filter BLAST results
    tcdb_blast = pd.read_csv(tcdb_blast_path, sep="\t")
    tcdb_blast["scovs"] = ((tcdb_blast["send"] - tcdb_blast["sstart"]) / tcdb_blast["slen"]) * 100
    tcdb_blast = tcdb_blast.loc[
        (tcdb_blast["pident"] > pident) &
        (tcdb_blast["qcovs"] > qcovs) &
        (abs(tcdb_blast["qlen"] - tcdb_blast["slen"]) / tcdb_blast["qlen"] < length_difference)]
    tcdb_blast = tcdb_blast.sort_values(by=["qacc", "bitscore"], ascending=[True, False]).groupby("qacc").head(1)
    tcdb_blast['proteinId'] = tcdb_blast['qacc']
    tcdb_blast['tcdb_blast'] = tcdb_blast["sacc"].apply(lambda x: x.split("|")[-1])
    tcdb_blast['tcdb_blast_family'] = tcdb_blast['tcdb_blast'].apply(lambda x: '.'.join(x.split('.')[:3]))
    tcdb = tcdb_blast[['proteinId', 'tcdb_blast', 'tcdb_blast_family', 'pident', 'qcovs']]


    ### 2. Read and filter HMM results
    if tcdb_hmm_path:
        open_func = gzip.open if tcdb_hmm_path.endswith(".gz") else open

        hmm_attribs = ['accession', 'bias', 'bitscore', 'description', 'cluster_num', 'domain_exp_num', 'domain_included_num', 'domain_obs_num', 'domain_reported_num', 'env_num', 'evalue', 'id', 'overlap_num', 'region_num']
        hits = defaultdict(list)

        with open_func(tcdb_hmm_path, "rt") as handle:
            for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
                for hit in queryresult.hits:
                    hits['domain'].append(queryresult.accession)
                    for attrib in hmm_attribs:
                        hits[attrib].append(getattr(hit, attrib))

        tcdb_hmm = pd.DataFrame.from_dict(hits)
        tcdb_hmm = tcdb_hmm.loc[tcdb_hmm['bitscore'] > hmm_bitscore]
        tcdb_hmm = tcdb_hmm.sort_values(by=["id", "bitscore"], ascending=[True, False]).groupby("id").head(1)
        tcdb_hmm['proteinId'] = tcdb_hmm['id']
        tcdb_hmm['tcdb_hmm'] = tcdb_hmm['domain'].str.split(".TD").str[0]
        tcdb_hmm['hmm_bitscore'] = tcdb_hmm['bitscore']

        # take annotation first from BLAST and then from hmmsearch
        tcdb = tcdb.merge(tcdb_hmm[['proteinId', 'tcdb_hmm', 'hmm_bitscore']], on="proteinId", how="left")
        tcdb['tcdb_annotation'] = tcdb['tcdb_blast']
        tcdb['tcdb_family'] = tcdb['tcdb_blast_family']
        tcdb.loc[tcdb['tcdb_annotation'].isnull(), 'tcdb_annotation'] = tcdb['tcdb_hmm']
        tcdb.loc[tcdb['tcdb_family'].isnull(), 'tcdb_family'] = tcdb['tcdb_hmm']

    # if there is no hmm take annotations only from BLAST
    else:
        tcdb['tcdb_annotation'] = tcdb['tcdb_blast']
        tcdb['tcdb_family'] = tcdb['tcdb_blast_family']

    # Merge with family descriptions
    tcdb_fams = pd.read_csv(tcdb_fams_path, sep="\t", header=None, names=['tcdb_family', 'tcdb_family_def'])
    tcdb = tcdb.merge(tcdb_fams, on='tcdb_family', how="left")

    # Add class definitions
    tcdb['tcdb_class'] = tcdb['tcdb_family'].str[0]
    data = {
        'tcdb_class': ['1', '2', '3', '4', '5', '8', '9'],
        'tcdb_class_def': [
            'Channels/Pores', 'Electrochemical Potential-driven Transporters', 'Primary Active Transporters',
            'Group Translocators', 'Transmembrane Electron Carriers', 'Accessory Factors Involved in Transport',
            'Incompletely Characterized Transport Systems'
        ]
    }
    tcdb_clas = pd.DataFrame(data)
    tcdb = tcdb.merge(tcdb_clas, on='tcdb_class', how="left")

    # Merge with substrate data
    tcdb_subs = pd.read_csv(tcdb_subs_path, sep="\t", header=None, names=['tcdb_annotation', 'tcdb_substrate'])
    tcdb = tcdb.merge(tcdb_subs, on='tcdb_annotation', how="left")

    return tcdb
