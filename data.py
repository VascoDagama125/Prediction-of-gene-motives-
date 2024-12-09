import pandas as pd

FILE_PATH = "/home/paweldyngosz/Dokumenty/Praca_Dyplomowa/arabidopsis_thaliana_tair_gtf.xlsx"
OUTPUT_FILE_PATH = "/home/paweldyngosz/Dokumenty/Praca_Dyplomowa/merged_file.csv"
SEPARATOR = '\t'
INCLUDE_HEADER = False
INCLUDE_INDEX = False

try:
    # Read the Excel file
    excel_data = pd.ExcelFile(FILE_PATH)

    try:
        gene_data = excel_data.parse(0)
        gene_ontology = excel_data.parse(1)
    except ValueError as e:
        print(f"Error parsing sheet: {e}")
        raise  # or handle accordingly

    # Filter rows where 'feature' is 'gene' and extract necessary fields
    if 'feature' not in gene_data.columns or 'attribute' not in gene_data.columns:
        raise KeyError("Necessary columns ('feature', 'attribute') not found in gene data.")

    is_gene = gene_data['feature'] == 'gene'
    gene_id = gene_data.loc[is_gene, 'attribute'].str.extract(r'ID=([^;]+)', expand=False)
    tss = gene_data.loc[is_gene, 'start']

    # Create dataframes and merge
    gene_frame = pd.DataFrame({'gene_id': gene_id, 'TSS': tss})
    gene_ontology_df = gene_ontology[['GO_ID', 'gene_id']]  # Directly select relevant columns

    if gene_frame.empty or gene_ontology_df.empty:
        raise ValueError("One of the datasets is empty.")

    # Merge the two DataFrames on 'gene_id'
    merged_df = pd.merge(gene_frame, gene_ontology_df, on='gene_id', how='inner')

    if merged_df.empty:
        print("Warning: Merged DataFrame is empty.")

    # Save the merged DataFrame to a CSV file
    merged_df.to_csv(OUTPUT_FILE_PATH, sep=SEPARATOR, header=True, index=INCLUDE_INDEX)

    print("Merging and saving completed successfully.")

except FileNotFoundError as e:
    print(f"File not found: {FILE_PATH}. {e}")
except pd.errors.EmptyDataError:
    print("No data found in the Excel file.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")
