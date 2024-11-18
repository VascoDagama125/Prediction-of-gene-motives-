import pandas as pd

file = "arabidopsis_thaliana_tair_gtf.xlsx"
df = pd.read_excel(file, sheet_name = 0)

#print(selected_columns_gtf. info())

#filter_df = df[df['feature'] == 'CDS']
#filtered_df filters df to include only rows where feature column's values are in the specified list
filtered_df = df[df['feature'].isin(['CDS', 'three_prime_UTR', 'five_prime_UTR'])]
#filters specified feature data and connects with start and end coordinates
selected_columns = filtered_df.loc[:,["feature", "start", "end"]]
print(selected_columns)
print(selected_columns.info())

#filter_df = df[df['feature'] == 'CDS']


