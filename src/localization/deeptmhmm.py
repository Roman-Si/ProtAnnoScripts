import pandas as pd
import os
import sys


def deeptmhmm(deeptmhmm_output_file):
    """
    Parses the DeepTMHMM annotations from the tabular .gff3 output
    
    Parameters:
    deeptmhmm_output_file (str): Path to the deeptmhmm gff3 output.
    
    Returns:
    pandas.DataFrame: DataFrame containing protein IDs, the annotations and their their start and stop positions.
    """ 
    df = pd.read_csv(deeptmhmm_output_file, sep='\t', header=None, names=['proteinId', 'Type', 'Start', 'End'], dtype={'proteinId': object})
    df.dropna(inplace=True)
    df['proteinId'] = df['proteinId'].str.split(' ').str[0]
    df['Start'] = df['Start'].astype(int).astype(str)
    df['End'] = df['End'].astype(int).astype(str)
    df = df.groupby('proteinId').agg({'Start': ';'.join, 'End': ';'.join, 'Type': '+'.join})
    df.columns = ['DeepTMHMM_start', 'DeepTMHMM_end', 'DeepTMHMM']
    return df.reset_index()