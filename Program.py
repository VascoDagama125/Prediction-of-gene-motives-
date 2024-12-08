import pandas as pd

# Load the data from an Excel file
file = "arabidopsis_thaliana_tair_gtf.xlsx"
tss_data = pd.read_excel(file, sheet_name=0)

# Extract GeneID from the 'attribute' column
tss_data = tss_data[tss_data['feature'] == 'gene']
tss_data['GeneID'] = tss_data['attribute'].str.extract(r'ID=([^;]+)')

# Adjust BED format fields:
# Assuming 'TSS' is the transcription start site (1-based), convert it to 0-based BED convention.

# Create a BED file with the desired columns
bed_data = tss_data[['seq_name', 'start', 'end', 'GeneID', 'strand']]

# Save to a BED file
bed_data.to_csv('newfile1.bed', sep='\t', header=False, index=False)

