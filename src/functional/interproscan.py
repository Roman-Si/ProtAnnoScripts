import pandas as pd


def interproscan(interproscan_output_file):
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