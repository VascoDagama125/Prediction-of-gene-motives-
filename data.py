import pandas as pd

FILE_PATH = "/home/paweldyngosz/Dokumenty/Praca_Dyplomowa/arabidopsis_thaliana_tair_gtf.xlsx"
OUTPUT_FILE_PATH = "/home/paweldyngosz/Dokumenty/Praca_Dyplomowa/merged_file.csv"
SEPARATOR = '\t'

INCLUDE_INDEX = False

try:
    # Odczyt pliku Excel
    excel_data = pd.ExcelFile(FILE_PATH)

    try:
        gene_data = excel_data.parse(0)
        gene_ontology = excel_data.parse(1)
    except ValueError as e:
        print(f"Error parsing sheet: {e}")
        raise  

    # Filtruj wiersze, w których „feature” to „gene” i wyodrębnia niezbędne pola.
    if 'feature' not in gene_data.columns or 'attribute' not in gene_data.columns:
        raise KeyError("Necessary columns ('feature', 'attribute', start, end, feature) not found in gene data.")
    is_gene = gene_data['feature'] == 'gene'
    #Filtruje attribute dla gene_id i gene feature
    gene_id = gene_data.loc[is_gene, 'attribute'].str.extract(r'ID=([^;]+)', expand=False)
    #Filtruje dane o chromosomach
    chromosome = gene_data.loc[is_gene, 'seq_name']
    #Filtruje dane o nici
    strand = gene_data.loc[is_gene, 'strand']
    #Filtruje dane o koordynatach start i end
    start = gene_data.loc[is_gene, 'start']
    end = gene_data.loc[is_gene, 'end']
    #Filtruje feature
    feature = gene_data.loc[is_gene, 'feature']
    # Tworzy dataframes i scala
    gene_frame = pd.DataFrame({'gene_id': gene_id, 'seq_name': chromosome, 'strand': strand, 'start': start, 'end': end, 'feature': feature })
    gene_ontology_df = gene_ontology[['GO_ID', 'gene_id']]  

    if gene_frame.empty or gene_ontology_df.empty:
        raise ValueError("One of the datasets is empty.")

    # Scal DataFrame na „gene_id
    merged_df = pd.merge(gene_frame, gene_ontology_df, on='gene_id', how='inner')
    # Usuwa duplikaty ze scalonej DataFrame
    merged_df = merged_df.drop_duplicates()
    
    if merged_df.empty:
        print("Warning: Merged DataFrame is empty.")

    # Zapisanie scalonej ramki danych do pliku CSV
    merged_df.to_csv(OUTPUT_FILE_PATH, sep=SEPARATOR, header=True, index=INCLUDE_INDEX)

    print("Merging and saving completed successfully.")

except FileNotFoundError as e:
    print(f"File not found: {FILE_PATH}. {e}")
except pd.errors.EmptyDataError:
    print("No data found in the Excel file.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")
