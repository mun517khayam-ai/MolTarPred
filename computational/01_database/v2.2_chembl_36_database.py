# Github repository code reference: chembl_34_database

import pandas as pd
import sqlalchemy
from sqlalchemy import text

# Read the file
file = '/Users/macbook/chembl/Code/v2/v2_chembl_36_filter.csv'
df = pd.read_csv(file)

# Check if all entries in the 'standard_units' column are 'nM'
if not df['standard_units'].eq('nM').all():
    raise ValueError("Not all entries in 'standard_units' are 'nM'")

# Filter rows where the value in 'standard_value' column is greater than 0
filtered_df = df[df['standard_value'] > 0][['molecule_chembl_id', 'target_chembl_id']]

# Connect to the PostgreSQL database
engine = sqlalchemy.create_engine('postgresql://macbook@localhost:5432/chembl_36')  # Use your own postgres password and localhost

# Execute SQL query to get molecule_chembl_id and SMILES
query = text("""
SELECT molecule_dictionary.chembl_id, compound_structures.canonical_smiles
FROM molecule_dictionary
JOIN compound_structures ON molecule_dictionary.molregno = compound_structures.molregno
""")
with engine.connect() as conn:
    result = conn.execute(query)
    chembl_smiles_df = pd.DataFrame(result.fetchall(), columns=['chembl_id', 'canonical_smiles'])

# Merge the filtered data with SMILES data
result_df = filtered_df.merge(chembl_smiles_df, left_on='molecule_chembl_id', right_on='chembl_id', how='left')

# Remove the chembl_id column and rename columns to Ligand-id, SMILES, and Target-id
result_df = result_df[['molecule_chembl_id', 'canonical_smiles', 'target_chembl_id']]
result_df.columns = ['Ligand-id', 'SMILES', 'Target-id']

# Remove 'CHEMBL' string
result_df['Ligand-id'] = result_df['Ligand-id'].str.replace('CHEMBL', '', regex=False)
result_df['Target-id'] = result_df['Target-id'].str.replace('CHEMBL', '', regex=False)

# Drop rows with missing values in the SMILES column
result_df = result_df.dropna(subset=['SMILES'])

# Merge Target-id by Ligand-id, separating with colons
result_df = result_df.groupby('Ligand-id').agg({
    'SMILES': 'first',  # Keep the first SMILES for each Ligand-id
    'Target-id': lambda x: ':'.join(x.unique())  # Concatenate Target-id
}).reset_index()

# Save the results to a CSV file
output_file_path = '/Users/macbook/chembl/Code/v2/v2_chembl_36_database.csv'
result_df.to_csv(output_file_path, index=False)

print("Standardized data exported to CSV file.")