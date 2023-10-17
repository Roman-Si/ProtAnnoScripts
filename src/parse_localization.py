import pandas as pd
import numpy as np

import os
import sys

from src.localization.deeptmhmm import deeptmhmm
from src.localization.signalp6 import signalp6
from src.localization.deeploc2 import deeploc2


def location( signalp_output_file, deeptmhmm_output_file, deeploc2_output_file, output_file_path):
    """
    Give outputs of SignalP6, DeepTMHMM and DeepLoc2. This function will classify these proteins with the following operations:
    Extracellular: DeepLoc2 | (SignalP6 & no TM ) 
    Transmembrane: DeepTMHMM
    Intracellular: All the rest
    
    Keep in mind 
    
    Returns the dataframe and saves it in output_file_path
    """  
    # parse the outputs
    signalp_data = signalp6(signalp_output_file)
    deeptmhmm_data = deeptmhmm(deeptmhmm_output_file)
    deeploc2_data = deeploc2(deeploc2_output_file)

    df = pd.merge(deeploc2_data[['proteinId', 'DeepLoc2', 'Signals_DeepLoc2']], signalp_data, on = 'proteinId', how='outer')
    df = pd.merge(df, deeptmhmm_data, on = 'proteinId', how='outer')
    
    # if you do not run DeepTMHMM for all proteins then there will be NaN values in these columns
    tm_vals = [s for s in df.loc[~df['DeepTMHMM'].isnull()]['DeepTMHMM'].unique() if "TMhelix" in s]
    periplasm_vals = [s for s in df.loc[~df['DeepTMHMM'].isnull()]['DeepTMHMM'].unique() if "periplasm" in s]
    
    
    # categorize
    df.loc[df['DeepLoc2'].str.contains("Extracellular"), 'Location'] = "Extracellular"
    df.loc[ (df['SP'] == 'Y') & (~df['DeepTMHMM'].isin(tm_vals + periplasm_vals)), 'Location'] = 'Extracellular'
    
    df.loc[df['DeepTMHMM'].isin(tm_vals), 'Location'] = 'Transmembrane'

    df['Location'] = df['Location'].fillna('Intracellular')

    df.to_csv(output_file_path, index=False)
    
    return(df)