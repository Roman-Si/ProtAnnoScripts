import pandas as pd


def kobas_annotate_to_df(kobas_annotate_file: str) -> pd.DataFrame:
    """
    Parse the kobas annotate output and return pandas dataframe with one row for each protein
   
    Parameters:
    kobas_annotate_file (str): Path to the kobas annotate output.
    
    Returns:
    pandas.DataFrame: DataFrame with kobas annotate output with columns "proteinId", "KEGG" that has the pathway names and "KEGG_orthologs"
    """ 
    
    proteinIDs = []
    pathways = []
    orthologs = []

    # Read the file and split it into sections using "////" delimiter
    with open(kobas_annotate_file, "r") as file:
        sections = file.read().split("////\n")[1:]

        # Iterate through sections to extract pathway information
        for section in sections:
            lines = section.strip().split("\n")
        
            # Extract gene IDs
            gene_id = [line.split("\t")[-1] for line in lines if line.startswith("Query:")]
        
            # Extract pathway information
            pathway_info = [line.split("\t") for line in lines if 'KEGG PATHWAY' in line]
        
            # If there is pathway information, process it
            if pathway_info:
                proteinIDs.append(gene_id[0])
                pathways.append([pathway[1] for pathway in pathway_info])
                orthologs.append([pathway[3] for pathway in pathway_info])

    df = pd.DataFrame({'proteinId': proteinIDs, "KEGG": pathways, "KEGG_orthologs": orthologs})
    # Join the list of pathways into a single string separated by ';'
    df['KEGG'] = df['KEGG'].apply(lambda pathways: ';'.join(pathways))
    # Join the list of KEGG orthologs into a single string separated by ';' after removal of duplicate ortholog IDs
    df['KEGG_orthologs'] = df['KEGG_orthologs'].apply(lambda orthologs: ';'.join(set(orthologs)))


    return(df)