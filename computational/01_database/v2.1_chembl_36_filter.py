# GitHub Repository code reference: Chembl_34_filter.py

import sqlalchemy
import pandas as pd

# Connect to the PostgreSQL database using SQLAlchemy
engine = sqlalchemy.create_engine('postgresql://macbook@localhost:5432/chembl_36') # Use your own postgres password and localhost

# Define the query to retrieve data
query = """
SELECT
    compound_records.compound_key AS compound_id,
    molecule_dictionary.chembl_id AS molecule_chembl_id,
    molecule_dictionary.molregno,
    activities.standard_value,
    activities.standard_units,
    activities.standard_type,
    target_dictionary.chembl_id AS target_chembl_id,
    target_dictionary.pref_name AS target_name,
    target_dictionary.organism AS target_organism
FROM
    activities
JOIN
    assays ON activities.assay_id = assays.assay_id
JOIN
    target_dictionary ON assays.tid = target_dictionary.tid
JOIN
    compound_records ON activities.record_id = compound_records.record_id
JOIN
    molecule_dictionary ON compound_records.molregno = molecule_dictionary.molregno
WHERE
    activities.standard_value <= 10000 AND
    activities.standard_units = 'nM' AND
    (activities.standard_type = 'IC50' OR
     activities.standard_type = 'Ki' OR
     activities.standard_type = 'EC50') AND
    -- 'Unchecked' (CHEMBL612545) is a curation placeholder for assays whose target was
    -- never assigned: no protein components, ~2.3M unrelated activities under one id.
    target_dictionary.target_type != 'UNCHECKED'
-- Postgres returns rows unordered, and the drop_duplicates below keeps the FIRST row of
-- each colliding group, so without a total order the exported set is not reproducible
-- (~23k molecule-target pairs churn between identical runs).
ORDER BY
    activities.activity_id
"""

# Execute the query and load the results into a Pandas DataFrame
df = pd.read_sql_query(query, engine)

# Data cleaning
# Exclude data on multi-protein complexes and non-specific experiments
df = df[~df['target_name'].str.contains('multiple|complex', case=False, na=False)]

# Remove duplicate entries.
# NB: compound_id is compound_records.compound_key -- the paper-local label ("1", "5a"),
# which 47,877 distinct molecules share. Deduplicating on it collapses unrelated compounds
# from different publications and destroyed ~497k real molecule-target pairs (37%).
# The molecule identity is molecule_chembl_id.
df = df.drop_duplicates(subset=['molecule_chembl_id', 'target_chembl_id'])

# Export the cleaned data to a CSV file
df.to_csv('/Users/macbook/chembl/Code/Github_target_prediction/Chembl_36_filter.csv', index=False)

# Print confirmation of file export
print("Filtered data exported to CSV file.")