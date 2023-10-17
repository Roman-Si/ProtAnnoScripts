import pandas as pd
import os
import sys


def netgpi(netgpi_output_file):
    """
    Takes NetGPI v1.1 output. Returns proteins with GPI-anchors.
    
    Parameters:
    netgpi_output_file (str): Path to the NetGPI tabular output
    
    Returns:
    pandas.DataFrame: DataFrame containing protein IDs and a column indicating the presence of a GPI-Anchor
    """   
    df = pd.read_csv(netgpi_output_file, skiprows=1, sep='\t', dtype={'# ID': object})
    df['proteinId'] = df['# ID']
    df['NetGPI'] = 'Y'
    df = df.loc[df['Pred. GPI-Anchored'] == 'GPI-Anchored'][['proteinId', 'NetGPI']]
    return df