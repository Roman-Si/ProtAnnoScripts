import pandas as pd
import numpy as np


def interproscan_to_df(interproscan_output_file: str) -> pd.DataFrame:
    """
    Parse the tsv interproscan output and return pandas dataframe 
   
    Parameters:
    interproscan_output_file (str): Path to the interproscan tsv output.
    
    Returns:
    pandas.DataFrame: DataFrame with interproscan output except the columns MD5, Status and Date
    """ 
    ipr_header = ['proteinId', 'MD5', 'Pr_length', 'Analysis', 'Sig_acc', 'Sig_descr', 'Start_loc', 'Stop_loc', 'E Value', 'Status', 'Date', 
                  'IPR_acc', 'IPR_descr', 'GO', 'Pathway']
    df = pd.read_csv(interproscan_output_file, sep = '\t', header=None, names=ipr_header, na_values=['-'], dtype={'proteinId': object})
    df['Domain_length'] = abs(df['Stop_loc'] - df['Start_loc'])
    df['Domain_cov'] = (df['Domain_length'] / df['Pr_length']) * 100
    df['Domain_cov'] = df['Domain_cov'].round(2)
    df.sort_values(by=['proteinId', 'E Value', 'Domain_cov'], ascending=[True, True, False], ignore_index=False, inplace = True)
    return(df[['proteinId', 'Pr_length', 'Domain_length', 'Domain_cov',  'Analysis', 'Sig_acc','Sig_descr', 'Start_loc', 'Stop_loc', 'E Value', 
               'IPR_acc', 'IPR_descr', 'GO', 'Pathway']])

def interproscan_df_to_tidy(interproscan_df: pd.DataFrame, interpro_entry_list: str) -> pd.DataFrame:
    """
    Make the interproscan dataframe of interproscan tidy. 
   
    Parameters:
    interproscan_df (pd.DataFrame): The dataframe from the interproscan function
    interpro_entry_list (str): Path to the entry list file, for the release 96 was downloaded from https://ftp.ebi.ac.uk/pub/databases/interpro/releases/96.0/entry.list
    
    Returns:
    pandas.DataFrame: DataFrame with one row per protein. Different InterPro annotation types are in respective columns. All the annotations without respective InterPro annotation are in columns "Sig_acc" and "Sig_descr"
    """ 

    ### Parse InterPro annotations
    # here only the interpro annotations
    ipr_interproscan = interproscan_df[~interproscan_df['IPR_acc'].isnull()][['proteinId', 'IPR_acc', 'IPR_descr']]
    ipr_interproscan.drop_duplicates(keep='first', inplace=True, ignore_index=True)
    
    # read the InterPro entry list
    interpro_db = pd.read_csv(interpro_entry_list, sep = "\t")
    interpro_db.columns = ['IPR_acc', 'ENTRY_TYPE', 'IPR_descr']

    # Concatenate unique 'IPR_acc' values for each 'proteinId'
    ipr_acc_concatenated = ipr_interproscan.groupby('proteinId')['IPR_acc'].unique().str.join('; ').reset_index()

    # Pivot the dataframe using pivot_table and aggregate duplicate values with a semicolon
    ipr_interproscan = ipr_interproscan.merge(interpro_db[['IPR_acc', 'ENTRY_TYPE']], on='IPR_acc', how='left')
    ipr_interproscan = ipr_interproscan.pivot_table(index='proteinId', columns='ENTRY_TYPE', values='IPR_descr', 
                                                aggfunc=lambda x: '; '.join(x))
    ipr_interproscan.reset_index(inplace=True)

    # Merge the concatenated 'IPR_acc' with the pivoted dataframe
    ipr_interproscan = ipr_interproscan.merge(ipr_acc_concatenated, on='proteinId', how='left')

    ### Take Sig annotations without InterPro annotation excluding Coils and MobiDBLite
    interproscan_sig = interproscan_df.loc[(interproscan_df['IPR_acc'].isnull()) & 
                                       (~interproscan_df['Analysis'].isin(["Coils", "MobiDBLite"]))][['proteinId', 'Sig_acc', 'Sig_descr']]

    interproscan_sig = interproscan_sig.replace(np.nan, '', regex=True)
    interproscan_sig = interproscan_sig.groupby('proteinId').agg(lambda x: ';'.join(set(x))).reset_index()

    ### Merge the dataframes 
    ipr_interproscan = ipr_interproscan.merge(interproscan_sig, on = "proteinId", how = "outer")

    return ipr_interproscan[['proteinId','Family', 'Homologous_superfamily', 'Domain',  'Active_site', 
                                'Binding_site', 'Conserved_site',  'PTM', 'Repeat', 'IPR_acc', 'Sig_acc', 'Sig_descr']]



def interproscan_df_to_go_df(interproscan_df: pd.DataFrame) -> pd.DataFrame:
    """
    Create a tidy GO dataframe with the GO terms of each proteinId from the interproscan dataframe.
   
    Parameters:
    interproscan_df (pd.DataFrame): The dataframe from the interproscan function
    
    Returns:
    pandas.DataFrame: DataFrame with one row per protein with GO annotations. Suitable format for GOATOOLS with ";" delimiter between GO annotations.
    """ 

    ### Parse InterPro annotations
    go_df = interproscan_df[['proteinId', 'GO']]
    go_df = go_df.loc[~go_df['GO'].isnull()]
    go_df = go_df.replace('\|', ";", regex = True)
    go_df = go_df.groupby('proteinId').agg(lambda x: ';'.join(set(x))).reset_index()
    # revome dublicate GO terms for a protein
    go_df['GO'] = go_df['GO'].str.split(';')
    go_df['GO'] = go_df['GO'].apply(lambda x: list(set(x)))
    go_df['GO'] = [';'.join(map(str, l)) for l in go_df['GO']]

    return go_df