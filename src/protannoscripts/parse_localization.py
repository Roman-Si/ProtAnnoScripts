import os
import pandas as pd

from .signalp6 import signalp6_to_df
from .deeploc2 import deeploc2_to_df
from .tmhmm import tmhmm_to_df
from .netgpi import netgpi_to_df


def annotate_location(
        unifire_df = None, 
        functional_info_dir = "data/functional/location/",  
        signalp_output_file = "signalp6/prediction_results.txt.gz", tmhmm_sp_output_file = "signalp6/tmhmm2.out.gz", netgpi_output_file = "signalp6/netgpi.out.gz",
        tmhmm_output_file = "tmhmm2/tmhmm2.out.gz",
        deeploc2_output_file = "deeploc2/deeploc2.csv.gz", deeploc_extracellular_threshold = 0.52753906, deeploc_membrane_threshold = 0.64638672) -> pd.DataFrame:
    """
    Give the path to the localization annotations directory, the relative path of each output file to this directory and a dataframe with the parsed UniFire annotations with columns ['proteinId', 'unifire_loc'] if available. 
    This function will classify these proteins with the following order:
    1. SignalP6 && NetGPI ->  Extracellular GPI-anchored
    2. DeepLoc2 Extracellular column > 0.52753906 -> Extracellular
    3. unifire_loc -> Intracellular, Extracellular, Membrane
    4. SignalP6 && noTM -> Extracellular
    5. TMHMM | DeepLoc2 -> Membrane
    6. Deeploc(Nucles|Cytoplasm) -> Intracellular
    7. Rest remains blank but is probably Intracellular
    
    Returns the dataframe with assigned annotation in 'Location' column but contains all the info as well.
    """

    ### Load data
    signalp_df = signalp6_to_df(os.path.join(functional_info_dir + signalp_output_file))
    netgpi_df = netgpi_to_df(os.path.join(functional_info_dir + netgpi_output_file))
    deeploc_df = deeploc2_to_df(os.path.join(functional_info_dir + deeploc2_output_file))
    # For proteins with signal peptide (SP) replace the initial tmhmm annotation with that of the SP-stripped sequence 
    tmhmm_df = tmhmm_to_df(os.path.join(functional_info_dir + tmhmm_output_file)).set_index('proteinId')
    tmhmm_sp_df = tmhmm_to_df(os.path.join(functional_info_dir + tmhmm_sp_output_file)).set_index('proteinId')
    tmhmm_df = tmhmm_sp_df.combine_first(tmhmm_df)
    tmhmm_df = tmhmm_df.reset_index()

    ### Merge data
    # If dataframe with unifire_loc is not available will start with the deeploc2 output
    if unifire_df is not None:
        location_df = unifire_df[['proteinId', 'unifire_loc']]
    else:
        location_df = deeploc_df[['proteinId']].copy()
        location_df['unifire_loc'] = pd.NA
        
    for dataframe in [signalp_df, netgpi_df, tmhmm_df[['proteinId', 'ExpAA', 'PredHel', 'Topology']], deeploc_df]:
        location_df = location_df.merge(dataframe, on = "proteinId", how = "left")

    ### Annotate with the rules
    # 1. NetGPI after SignalP6 -> Extracellular;GPI
    location_df.loc[(location_df["NetGPI"] == "Y") , "Location"] = "Extracellular;GPI"
    # 2. Filter for deeploc_extracellular_threhold and select those that start with Extracellular
    location_df.loc[(location_df['Extracellular'] > deeploc_extracellular_threshold) & ( location_df['DeepLoc2'].str.startswith("Extracellular") ) & (location_df['Location'].isnull()), 'Location'] = location_df['DeepLoc2']
    # 3. UniFire
    if unifire_df is not None:
        location_df.loc[location_df['Location'].isnull(), 'Location']  = location_df["unifire_loc"]
    # 4. SignalP & no TMhelix -> Extracellular
    location_df.loc[(location_df['SignalP6'] == 'SP') & (location_df['PredHel'] == 0) & (location_df['Location'].isnull()), 'Location'] = 'Extracellular'
    # 5. TMHMM | DeepLoc2 for membrane proteins
    location_df.loc[(location_df['Location'].isnull()) & (location_df['Cell membrane'] > deeploc_membrane_threshold) & ( location_df['DeepLoc2'].str.startswith("Cell membrane") ), 'Location'] = 'Membrane'
    location_df.loc[(location_df['Location'].isnull()) & (location_df['PredHel'] > 0) , 'Location'] = 'Membrane'
    # 6. Intracellular: Proteins without any other annotation and DeepLoc2 starts with Cytoplasm | Nucleus
    location_df.loc[ (location_df['Location'].isnull()) & ( ( location_df['DeepLoc2'].str.startswith("Cytoplasm") ) | ( location_df['DeepLoc2'].str.startswith("Nucleus") ) ), "Location"] = "Intracellular"
    
    return location_df
