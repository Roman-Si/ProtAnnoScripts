import pandas as pd


def deeptmhmm_to_df(deeptmhmm_output_file: str) -> pd.DataFrame:
    """
    Parses the DeepTMHMM annotations from the tabular .gff3 output
    
    Parameters:
    deeptmhmm_output_file (str): Path to the deeptmhmm gff3 output.
    
    Returns:
    pandas.DataFrame: DataFrame containing protein IDs, the annotations and their their start and stop positions.
    """
    df = pd.read_csv(deeptmhmm_output_file, sep='\t', header=None, comment='#', names=['proteinId', 'Type', 'Start', 'End'], dtype={'proteinId': object})
    df.dropna(inplace=True)
    df['proteinId'] = df['proteinId'].str.split(' ').str[0]
    df['Start'] = df['Start'].astype(int).astype(str)
    df['End'] = df['End'].astype(int).astype(str)
    df = df.groupby('proteinId').agg({'Start': ';'.join, 'End': ';'.join, 'Type': '+'.join})
    df.columns = ['DeepTMHMM_start', 'DeepTMHMM_end', 'DeepTMHMM']
    return df.reset_index()

def tmhmm_to_df(tmhmm_output_file: str) -> pd.DataFrame:
    """
    Parses the TMHMM annotations from the tabular .gff3 output
    
    Parameters:
    tmhmm_output_file (str): Path to the tmhmm tabular output.
    
    Returns:
    pandas.DataFrame: DataFrame containing protein IDs, the annotations and their their start and stop positions.
    """
    df = pd.read_csv(tmhmm_output_file, sep='\t', header=None,  names=['proteinId', 'len',  'ExpAA','First60', 'PredHel', 'Topology'], dtype={'proteinId': object})
    # every column value contains the column name followed by =, remove them
    for col in df.columns:
        if col != 'proteinId':
            df[col] = df[col].str.replace(f'{col}=', '')

    # Converting data types if necessary
    df['len'] = df['len'].astype(int)
    df['ExpAA'] = df['ExpAA'].astype(float)
    df['First60'] = df['First60'].astype(float)
    df['PredHel'] = df['PredHel'].astype(int)
    return df