import pandas as pd
import numpy as np
import os
import sys


def deeploc2(deeploc2_output_file):
    """
    Parses the DeepLoc2 csv output and returns the classifications of these proteins as a dataframe.
    
    Parameters:
    deeploc2_output_file (str): Path to the DeepLoc2 csv output file.
    
    Returns:
    pandas.DataFrame: DataFrame containing protein IDs, a column with the prediction, a column with the signals and then probabilities for all compartments
    """ 
    df = pd.read_csv(deeploc2_output_file, dtype={'Protein_ID': object})
    df['DeepLoc2'] = df['Localizations']
    df['proteinId'] = df['Protein_ID']
    df['Signals_DeepLoc2'] = df['Signals']
    return(df[['proteinId', 'DeepLoc2', 'Signals_DeepLoc2', 'Cytoplasm', 'Nucleus',
       'Extracellular', 'Cell membrane', 'Mitochondrion', 'Plastid',
       'Endoplasmic reticulum', 'Lysosome/Vacuole', 'Golgi apparatus',
       'Peroxisome']])