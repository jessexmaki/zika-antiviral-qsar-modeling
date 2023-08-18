# Install Miniconda3 and setup environment
# Download Miniconda installer script for Linux
#! wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh

# Make the installer script executable
#! chmod +x Miniconda3-py37_4.8.2-Linux-x86_64.sh

# Run the installer in batch mode, force overwriting if necessary, and install to /usr/local
#! bash ./Miniconda3-py37_4.8.2-Linux-x86_64.sh -b -f -p /usr/local

# Install rdkit package using conda
#! conda install -c rdkit rdkit -y

# Add the installed packages to system path
import sys
sys.path.append('/usr/local/lib/python3.7/site-packages/')

# Import required libraries
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

# Load the dataset
df = pd.read_csv('Data/aromatase_03_bioactivity_data_curated.csv')

# Remove the canonical_smiles column
df_no_smiles = df.drop(columns='canonical_smiles')

# Extract and process SMILES data
smiles = []
for i in df.canonical_smiles.tolist():
  cpd = str(i).split('.')  # Splitting at '.' to get individual compounds
  cpd_longest = max(cpd, key = len)  # Get the longest compound (in terms of characters)
  smiles.append(cpd_longest)
smiles = pd.Series(smiles, name = 'canonical_smiles')

# Combine the processed SMILES data with the original dataframe
df_clean_smiles = pd.concat([df_no_smiles, smiles], axis=1)

# Function to compute Lipinski descriptors
def lipinski(smiles, verbose=False):
    moldata= []
    for elem in smiles:
        mol=Chem.MolFromSmiles(elem) 
        moldata.append(mol)
       
    baseData= np.arange(1,1)
    i=0  
    for mol in moldata:        
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
        row = np.array([desc_MolWt, desc_MolLogP, desc_NumHDonors, desc_NumHAcceptors])   
    
        if(i==0):
            baseData=row
        else:
            baseData=np.vstack([baseData, row])
        i=i+1      
    
    columnNames=["MW","LogP","NumHDonors","NumHAcceptors"]
    descriptors = pd.DataFrame(data=baseData,columns=columnNames)
    
    return descriptors

# Apply the lipinski function
df_lipinski = lipinski(df_clean_smiles.canonical_smiles)

# Combine the original data with the Lipinski descriptors
df_combined = pd.concat([df, df_lipinski], axis=1)

# Get the descriptive statistics for the standard_value column
df_combined.standard_value.describe()

# Functions to normalize the standard values and calculate pIC50
def norm_value(input):
    norm = []
    for i in input['standard_value']:
        if i > 100000000:
            i = 100000000  # Cap the maximum value at 100,000,000
        norm.append(i)
    input['standard_value_norm'] = norm
    return input.drop('standard_value', axis=1)

def pIC50(input):
    pIC50 = []
    for i in input['standard_value_norm']:
        molar = i*(10**-9)  # Convert nM to M
        pIC50.append(-np.log10(molar))  # Compute pIC50 value
    input['pIC50'] = pIC50
    return input.drop('standard_value_norm', axis=1)

# Normalize the standard values and calculate pIC50 values
df_norm = norm_value(df_combined)
df_final = pIC50(df_norm)

# Save the cleaned data
df_final.to_csv('Data/aromatase_04_bioactivity_data_3class_pIC50.csv')

# Filter out intermediate class and save the 2-class dataset
df_2class = df_final[df_final['class'] != 'intermediate']
df_2class.to_csv('Data/aromatase_05_bioactivity_data_2class_pIC50.csv')
