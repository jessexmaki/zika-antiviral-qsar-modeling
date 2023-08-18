# The following lines use wget to download required files from the given URLs.
# wget is a free utility that allows for non-interactive downloading of files from the web.
#! wget https://github.com/dataprofessor/bioinformatics/raw/master/padel.zip
#! wget https://github.com/dataprofessor/bioinformatics/raw/master/padel.sh

# Importing the pandas library, which is useful for data manipulation and analysis.
import pandas as pd

# Read the bioactivity data from the CSV file into a dataframe.
df3 = pd.read_csv('Data/aromatase_04_bioactivity_data_3class_pIC50.csv')

# Select the 'canonical_smiles' and 'molecule_chembl_id' columns from the dataframe.
selection = ['canonical_smiles','molecule_chembl_id']
df3_selection = df3[selection]

# Save the selected data to a new tab-separated file without headers.
df3_selection.to_csv('molecule.smi', sep='\t', index=False, header=False)

# Read descriptor data from the CSV file into another dataframe.
df3_X = pd.read_csv('descriptors_output.csv')

# Extract the 'pIC50' column from the original dataframe.
df3_Y = df3['pIC50']

# Combine the descriptor data and 'pIC50' column into a single dataframe.
dataset3 = pd.concat([df3_X,df3_Y], axis=1)

# Display the combined dataframe.
dataset3

# Save the combined data to a new CSV file.
dataset3.to_csv('Data/aromatase_06_bioactivity_data_3class_pIC50_pubchem_fp.csv', index=False)
