import pandas as pd
import numpy as np
import os
import sys


def signalp6(signalp6_output_file):
    """
    Filters proteins from SignalP6 output to keep proteins with signal peptides.
    
    Parameters:
    signalp6_output_file (str): Path to the SignalP6 tabular output prediction_results.txt.
    
    Returns:
    pandas.DataFrame: DataFrame containing protein IDs, a column indicating the presence of signal peptides and a column with the cleavage positions.
    """ 
    df = pd.read_csv(signalp6_output_file, sep = '\t', dtype={'# ID': object}, skiprows=1)
    df = df.loc[df['Prediction'] == 'SP']
    df['proteinId'] = df['# ID'].str.split(' ').str[0]
    df['SP'] = 'Y'
    # Now take CS positions
    # In case of low SP probability there is no CS position given and the column is empty
    df['SP_cs'] = df['CS Position'].apply(lambda x: x.split(': ')[1].split('.')[0] if isinstance(x, str) else np.nan)
    return(df[['proteinId', 'SP', 'SP_cs']])
