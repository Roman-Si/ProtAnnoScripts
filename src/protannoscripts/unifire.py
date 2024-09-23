import pandas as pd
import numpy as np
pd.options.mode.copy_on_write = True


def process_unifire_dataframe(df: pd.DataFrame):
    """
    Takes as input an ARBA or UniRule dataframe. 
    
    Returns:
        Tuple of DataFrames: protein_names, protein_families, ec_numbers, localizations, enzyme_keywords, gene_names
    

    ARBA does not give gene names.
    Example usage:
    >>> arba_name, arba_fam, arba_ec, arba_loc, arba_keyws, arba_genes = process_unifire_dataframe(arba_df)
    >>> unirule_name, unirule_fam, unirule_ec, unirule_loc, unirule_keyws, unirule_gene = process_unifire_dataframe(unirule_df)

    """
    # get the first recommended protein name
    protein_names = df.loc[df['AnnotationType'] == 'protein.recommendedName.fullName'].drop_duplicates(subset='ProteinId', keep='first', ignore_index=False)
    # get gene names, only in unirule
    if "gene.name.primary" in df['AnnotationType'].unique():
        gene_names = df.loc[df['AnnotationType'] == "gene.name.primary"].drop_duplicates(subset='ProteinId', keep='first', ignore_index=False)
    else:
        gene_names = None
    
    # get the protein family
    protein_families = df.loc[df['AnnotationType'] == 'comment.similarity']
    protein_families.loc[:, 'unifire_family'] = protein_families['Value'].str.split('elongs to the ').str[1].str.split(r"\. ").str[0] # the "\." split is for superfamilies where you want the first sentence with the family only
    # get EC
    ec_numbers = df.loc[df['AnnotationType'] == 'protein.recommendedName.ecNumber'].drop_duplicates(subset=['ProteinId', 'Value'], keep='first', ignore_index=False)[['ProteinId', 'Value']]
    # get localization-related columns
    localizations = df.loc[df['AnnotationType'].isin(['comment.subcellular_location', 'feature.TOPO_DOM', 'feature.INTRAMEM'])]
    # get enzyme keywords
    enzyme_keywords = df.loc[df['AnnotationType'] == 'keyword']
    enzyme_keywords = enzyme_keywords.loc[enzyme_keywords['Value'].str.endswith('ase') & (enzyme_keywords['AnnotationType'] == 'keyword')][['ProteinId', 'Value']]
    
    return protein_names, protein_families[['ProteinId', 'unifire_family']], ec_numbers, localizations, enzyme_keywords, gene_names


def assign_location(df: pd.DataFrame,location_df: pd.DataFrame, words: list, location_type: str):
    """
    Give a dataframe with ProteinIds, list of words and associated location type (Intracellular, Membrane, Extracellular)
    In case of already existing unifire_loc annotation append the new one with ";"

    Args:
        df (dataframe): DataFrame with ProteinId column of all proteins present in unifire output
        location_df (dataframe): DataFrame with the location-related rows of unirule and ARBA
        words (list): List with keywords related with the given location_type
        location_type (str): Location type, one of ["Intracellular", "Membrane", "Extracellular"]
        
    
    Searches the dataframe with the unifire localization-related rows for these words and assign unifire_loc in the given dataframe
    """
    matching_proteins = location_df[location_df['Value'].str.lower().str.contains('|'.join(words))] # 'I' is an OR operator within str.contains
    
    # exclude the cytoplasmic vesicles from the intracellular annotations
    if location_type == "Intracellular":
        matching_proteins = matching_proteins[~matching_proteins['Value'].str.lower().str.contains("cytoplasmic vesicle")] 
    
    matching_protein_ids = matching_proteins['ProteinId']
    
    # Check if 'unifire_loc' column exists, and create it if not
    if 'unifire_loc' not in df:
        df['unifire_loc'] = ''
        
    # Iterate through matching proteins
    for protein_id in matching_protein_ids:
        if protein_id in df['ProteinId'].values:
            
            new_loc = location_type
            
            # Check if 'unifire_loc' is already assigned
            existing_loc = df.loc[df['ProteinId'] == protein_id, 'unifire_loc'].values[0]
            # if it exists combine the two locations
            if existing_loc and location_type not in existing_loc:
                new_loc = existing_loc + ';' + new_loc
            # assign the new_loc
            df.loc[df['ProteinId'] == protein_id, 'unifire_loc'] = new_loc
    
    

def merge_unifire_annotations(arba: str, unirule: str) -> pd.DataFrame:
    """
    Merges ARBA and UniRule dataframes and extracts UniFire annotations.

    Args:
        arba (str): Path to the ARBA output.
        unirule (str): Path to the UniRule output.
    
    Returns:
        pandas.DataFrame: DataFrame with UniFire annotations: protein name, gene name, protein family, EC, localization, enzyme keywords
    """
    
    # read the arba and unirule dataframes
    arba_df = pd.read_csv(arba, sep='\t', dtype={'ProteinId': object})
    unirule_df = pd.read_csv(unirule, sep='\t', dtype={'ProteinId': object})
    
    # create a dataframe with all detected ProteinId
    protein_ids = list(set(arba_df['ProteinId']).union(set(unirule_df['ProteinId'])))
    df = pd.DataFrame({'ProteinId': protein_ids})

    
    arba_name, arba_fam, arba_ec, arba_loc, arba_keyws, arba_genes = process_unifire_dataframe(arba_df)
    unirule_name, unirule_fam, unirule_ec, unirule_loc, unirule_keyws, unirule_genes = process_unifire_dataframe(unirule_df)

    ### Location
    location = pd.concat([arba_loc, unirule_loc], axis=0, ignore_index=True)
    
    assign_location(df, location, ['secret', 'extracellular', 'cell surface', 'cell wall'], "Extracellular")
    assign_location(df, location,['membrane'], "Membrane")
    assign_location(df,location, ['cytoplasm', 'nucleus'], "Intracellular")
    df.replace({'unifire_loc' : ""},  np.nan, inplace=True)
    
    ### Process names, family, EC, and keywords here
    
    # Names
    #give priority to unirule
    names = pd.concat([unirule_name, arba_name], axis=0).drop_duplicates(subset='ProteinId',keep='first', inplace=False, ignore_index=False)
    names.loc[:, 'proteinName'] = names['Value']

    # Genes
    # ARBA does not give gene names.
    genes = unirule_genes.drop_duplicates(subset='ProteinId',keep='first', inplace=False, ignore_index=False)
    genes.loc[:, 'geneName'] = genes['Value']

    # Family
    family = pd.concat([unirule_fam, arba_fam], axis=0)
    family = family.groupby('ProteinId').agg(lambda x: ';'.join(set(x))).reset_index()
    
    # EC
    ec = pd.concat([unirule_ec, arba_ec], axis=0).drop_duplicates(subset=['ProteinId', 'Value'],keep='first', inplace=False, ignore_index=False)
    ec = ec.groupby('ProteinId').agg(lambda x: ';'.join(set(x))).reset_index()
    ec.loc[:, 'unifire_EC'] = ec['Value']
    
    # Keywords
    keyws = pd.concat([unirule_keyws, arba_keyws], axis=0).drop_duplicates(subset=['ProteinId', 'Value'],keep='first', inplace=False, ignore_index=False)
    keyws = keyws.groupby('ProteinId').agg(lambda x: ';'.join(set(x))).reset_index()
    keyws.loc[:, 'enzyme_keywords'] = keyws['Value']
    
    # merge them to the df dataframe
    merge_list = [names[['ProteinId', 'proteinName']], genes[['ProteinId', 'geneName']], family, ec[['ProteinId', 'unifire_EC']], keyws[['ProteinId', 'enzyme_keywords']]]
    for data in merge_list:
        df = pd.merge(df, data, how='left', on='ProteinId')
    
    
    return df


def parse_unifire_for_gag(arba: str, unirule: str) -> pd.DataFrame:
    """
    Extracts gene and product names from unifire for GAG. 
    If for some reason you have only unirule output just give it twice.
    Steps:
    Merges ARBA and UniRule dataframes with merge_unifire_annotations.
    Prioritizing UniRule extracts gene name if relevant, converts it to lowercase and converts relevant proteinId to geneId
                         extracts product from proteinName
                         if no proteinName extracts product from family annotation

    Args:
        arba (str): Path to the ARBA output.
        unirule (str): Path to the UniRule output.
    
    Returns:
        pandas.DataFrame: DataFrame with UniFire annotations: protein name, gene name, protein family, EC, localization, enzyme keywords
    """
    
    df = merge_unifire_annotations(arba, unirule)
    
    long_df = pd.melt(df[['ProteinId', 'proteinName', 'geneName','unifire_family' ]], id_vars=['ProteinId'], 
                  value_vars=['proteinName', 'geneName', 'unifire_family'],
                  var_name='type', value_name='sacc')
    long_df = long_df[~long_df['sacc'].isnull()]
    
    # For family annotations ending with word "family" append "protein"
    long_df.loc[long_df['type'] == "unifire_family", 'sacc'] = long_df.loc[long_df['type'] == "unifire_family", 'sacc'].str.replace(r'family', 'family protein', regex=False)
    # conver proteinId to geneId for gene names
    long_df.loc[long_df['type'] == "geneName", 'ProteinId'] = long_df.loc[long_df['type'] == "geneName", 'ProteinId'].str.split("-").str[0]
    # for transporters remove the tcdb acession
    long_df['sacc'] = long_df['sacc'].str.replace(r'\(TC .*[0-9]\)', '', regex=True)
    # now change proteinName and family to product and keep first
    long_df['type'] = long_df['type'].replace({'proteinName': 'product', 'geneName': 'name', 'unifire_family' : 'product'})

    long_df.columns = ['qacc', 'type', 'sacc']
    # split multiple family annotations in separate rows
    long_df = long_df.assign(sacc=long_df['sacc'].str.split(';')).explode('sacc', ignore_index=True)

    ### Some product polishing
    # These are unifire specific
    # remove complex subunit annotation
    long_df = long_df[~long_df['sacc'].str.contains(r'complex .* kDa subunit family protein', na=False)]
    # remove the kDa from product name
    long_df['sacc'] = long_df['sacc'].str.replace(r' \d+ kDa ', ' ', regex=True)
    # Remove everything after " homolog"
    long_df['sacc'] = long_df['sacc'].str.replace(r' homolog.*', '', regex=True)
    # Remove "fungal " and everything before it
    long_df['sacc'] = long_df['sacc'].str.replace(r'.*(fungal|archaeal|eukaryotic) ', '', regex=True, case = False)

    long_df = long_df.drop_duplicates(subset=['qacc', 'type'], keep='first')
    return long_df

    
