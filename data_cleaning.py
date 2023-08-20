import pandas as pd

# Load raw bioactivity data from the CSV file
df = pd.read_csv('Data/bioactivity_data.csv')

# Filter out rows where standard_value or canonical_smiles is NaN
df2 = df.dropna(subset=['standard_value', 'canonical_smiles'])

# Remove duplicates based on canonical_smiles
df2_nr = df2.drop_duplicates(['canonical_smiles'])

# Select relevant columns
selection = ['molecule_chembl_id', 'canonical_smiles', 'standard_value', 'name']
df3 = df2_nr[selection]

# Save the preprocessed data to a new CSV file
df3.to_csv('Data/aromatase_02_bioactivity_data_preprocessed.csv', index=False)

# Reload the preprocessed data (Note: this step might be unnecessary if you continue using df3 directly)
df4 = pd.read_csv('Data/aromatase_02_bioactivity_data_preprocessed.csv')

# Determine the bioactivity class based on standard_value
def classify_bioactivity(value):
    if value >= 10000:
        return "inactive"
    elif value <= 1000:
        return "active"
    else:
        return "intermediate"

bioactivity_class = df4['standard_value'].apply(classify_bioactivity)

# Append the bioactivity class to the dataframe
df5 = df4.copy()
df5['class'] = bioactivity_class

# Save the curated data to a new CSV file
df5.to_csv('Data/aromatase_03_bioactivity_data_curated.csv', index=False)
