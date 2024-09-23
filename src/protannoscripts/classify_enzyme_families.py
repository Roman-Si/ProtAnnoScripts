import json
import pandas as pd


def search_annocol_for_keywords(
    df: pd.DataFrame, anno_column: str, keywords: list, avoid_keywords: list
) -> pd.DataFrame:
    """
    Filters a DataFrame based on keywords and avoid keywords in a specified annotation column.

    Args:
        df (pd.DataFrame): The DataFrame to filter.
        anno_column (str): The column in which to search for keywords.
        keywords (list): List of keywords to include in results.
        avoid_keywords (list): List of keywords to exclude from results.

    Returns:
        pd.DataFrame: Filtered DataFrame.
    """
    unique_annotations = df[anno_column].dropna().unique()

    matching_names = [
        annotation
        for annotation in unique_annotations
        if any(keyword in str(annotation).lower() for keyword in keywords)
    ]

    to_remove = [
        annotation
        for annotation in matching_names
        if any(keyword in str(annotation).lower() for keyword in avoid_keywords)
    ]

    df = df.loc[
        df[anno_column].isin(
            [annotation for annotation in matching_names if annotation not in to_remove]
        )
    ]

    return df


def classify_proteins_to_families(
    df: pd.DataFrame,
    protein_id_column: str,
    anno_column: str,
    enzyme_families_keywords_json: str,
) -> dict:
    """
    Classify proteins in df Dataframe to families in the categories JSON file.

    The classification follows a hierarchical order specified in the JSON file. Once a protein is assigned to a family, it is not assigned to any subsequent family, even if relevant keywords are present.

    Args:
        df (pd.DataFrame): The DataFrame to filter.
        protein_id_column (str): The column with the protein IDs.
        anno_column (str): The column in which to search for keywords.
        enzyme_families_keywords_json (str): Path to json file with enzyme families and their keywords. Example file in ../data/enzyme_family_keywords.json

    Returns:
        dict: A dictionary where each key represents an enzyme family, and the associated value is a list of classified proteins.

    Example usage:
    >>> family_assignments = classify_proteins_to_families(interproscan_df, "proteinId", "IPR_descr", "../data/enzyme_family_keywords.json")
    >>> print(family_assignments)
    {'ureases': ['protein1', 'protein2'], 'amidases': ['protein3', 'protein4'], ...}
    """

    with open(enzyme_families_keywords_json, "r") as file:
        enzyme_categories = json.load(file)

    classified_proteins = set()  # To track already classified proteins
    family_assignments = {}

    for family, family_data in enzyme_categories.items():
        keywords = family_data["keywords"]
        avoid_keywords = family_data["avoid_keywords"]

        # Filter proteins for the current family
        family_proteins = (
            search_annocol_for_keywords(df, anno_column, keywords, avoid_keywords)[
                protein_id_column
            ]
            .unique()
            .tolist()
        )

        # Remove already classified proteins
        unclassified_proteins = [
            protein for protein in family_proteins if protein not in classified_proteins
        ]

        # Assign proteins to the family
        family_assignments[family] = unclassified_proteins

        # Update the set of classified proteins
        classified_proteins.update(unclassified_proteins)

    return family_assignments


def merge_protein_family_classifications(
    dataframes_with_columns, enzyme_families_keywords_json: str
) -> dict:
    """
    Classify proteins using the provided dataframes and columns.

    The classification follows a hierarchical order following the order of the dataframes.


    Args:
        dataframes_with_columns (list of tuples): List of tuples where each tuple contains:
            - dataframe (pandas.DataFrame): The dataframe to classify proteins.
            - protein_id_column (str): The name of the column containing protein IDs.
            - anno_column (str): The name of the column with annotations.
        enzyme_families_keywords_json (str): Path to the JSON file with enzyme families and their keywords.

    Returns:
        dict: A dictionary where each key represents an enzyme family, and the associated value is a list of classified proteins.


    Example usage:
    # Gather the data
    >>> sprot_df = parse_sprot_fasta(swissprot_fasta_file)
    >>> blast_sprot_df = blast_sprot(blast_sprot_output_file, sprot_df)
    >>> unifire_df = merge_unifire_annotations(arba_output_file, unirule_output_file)
    >>> interproscan_df = interproscan(interproscan_output_file)

    # Recommended hierarchical order for dataframes
    >>> dataframes_with_columns = [
                                    (unifire_df, 'ProteinId', 'proteinName'),
                                    (blast_swissprot, 'proteinId', 'Protein name'),
                                    (unifire_df, 'ProteinId', 'unifire_family'),
                                    (unifire_df, 'ProteinId', 'enzyme_keywords'),
                                    (interproscan_df, 'proteinId', 'IPR_descr'),
                                    (interproscan_df, 'proteinId', 'Sig_descr')]
    >>> enzyme_classifications_dict = merge_protein_family_classifications(dataframes_with_columns, "../data/enzyme_family_keywords.json")
    """
    # Initialize a dictionary to store classified proteins by family
    classified_proteins = {}

    # Set to store proteins classified at least once
    classified_proteins_set = set()

    # Iterate through the provided dataframes
    for df, protein_id_column, anno_column in dataframes_with_columns:
        family_assignments = classify_proteins_to_families(
            df, protein_id_column, anno_column, enzyme_families_keywords_json
        )

        # Update the classified_proteins dictionary with unclassified proteins
        for family, unclassified_proteins in family_assignments.items():
            # Filter out proteins that have been classified in any family
            unclassified_proteins = [
                protein
                for protein in unclassified_proteins
                if protein not in classified_proteins_set
            ]

            if family not in classified_proteins:
                classified_proteins[family] = unclassified_proteins
            else:
                # Append unclassified proteins
                classified_proteins[family].extend(unclassified_proteins)

            # Update the set of proteins classified at least once
            classified_proteins_set.update(unclassified_proteins)

    return classified_proteins
