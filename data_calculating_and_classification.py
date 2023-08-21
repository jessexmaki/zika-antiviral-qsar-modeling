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
df = pd.read_csv('Data/zika_03_bioactivity_data_curated.csv')

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
df_final.to_csv('Data/zika_04_bioactivity_data_3class_pIC50.csv')

# Filter out intermediate class and save the 2-class dataset
df_2class = df_final[df_final['class'] != 'intermediate']
df_2class.to_csv('Data/zika_05_bioactivity_data_2class_pIC50.csv')

# Importing necessary libraries
import seaborn as sns
sns.set(style='ticks')
import matplotlib.pyplot as plt

# Set the default size of all plots
plt.figure(figsize=(5.5, 5.5))

# Plotting a bar chart for frequency of each bioactivity class
sns.countplot(x='class', data=df_2class, edgecolor='black')

# Adding labels for the x and y axes
plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('Frequency', fontsize=14, fontweight='bold')

# Saving the plot as a PDF
plt.savefig('Data/plot_bioactivity_class.pdf')

# Set the default size for the next plot
plt.figure(figsize=(5.5, 5.5))

# Scatter plot of MW vs LogP, colored by class and sized by pIC50
sns.scatterplot(x='MW', y='LogP', data=df_2class, hue='class', size='pIC50', edgecolor='black', alpha=0.7)

# Adding labels for the x and y axes
plt.xlabel('MW', fontsize=14, fontweight='bold')
plt.ylabel('LogP', fontsize=14, fontweight='bold')

# Adding a legend to the plot
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)

# Set the default size for the next plot
plt.figure(figsize=(5.5, 5.5))

# Boxplot for pIC50 values across bioactivity classes
sns.boxplot(x = 'class', y = 'pIC50', data = df_2class)

# Adding labels for the x and y axes
plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('pIC50 value', fontsize=14, fontweight='bold')

# Saving the plot as a PDF
plt.savefig('Data/plot_ic50.pdf')


# Define a function for performing the Mann-Whitney U test
def mannwhitney(descriptor, verbose=False):
    # Import necessary libraries and modules
    from numpy.random import seed
    from numpy.random import randn
    from scipy.stats import mannwhitneyu

    # Seed the random number generator for reproducibility
    seed(1)

    # Extract active samples for the given descriptor
    selection = [descriptor, 'class']
    df = df_2class[selection]
    active = df[df['class'] == 'active']
    active = active[descriptor]

    # Extract inactive samples for the given descriptor
    selection = [descriptor, 'class']
    df = df_2class[selection]
    inactive = df[df['class'] == 'inactive']
    inactive = inactive[descriptor]

    # Perform the Mann-Whitney U test
    stat, p = mannwhitneyu(active, inactive)

    # Interpret the result
    alpha = 0.05
    if p > alpha:
        interpretation = 'Same distribution (fail to reject H0)'
    else:
        interpretation = 'Different distribution (reject H0)'

    # Compile the results into a DataFrame
    results = pd.DataFrame({'Descriptor':descriptor,
                            'Statistics':stat,
                            'p':p,
                            'alpha':alpha,
                            'Interpretation':interpretation}, index=[0])
    
    # Save the results as a CSV file
    filename = 'Data/mannwhitneyu_' + descriptor + '.csv'
    results.to_csv(filename)

    return results

# Run the Mann-Whitney U test for pIC50 descriptor
mannwhitney('pIC50')

# Plotting boxplots for other descriptors as well, following the same pattern
plt.figure(figsize=(5.5, 5.5))
sns.boxplot(x = 'class', y = 'MW', data = df_2class)
plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('MW', fontsize=14, fontweight='bold')
plt.savefig('Data/plot_MW.pdf')
mannwhitney('MW')

plt.figure(figsize=(5.5, 5.5))
sns.boxplot(x = 'class', y = 'MW', data = df_2class)
plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('MW', fontsize=14, fontweight='bold')
plt.savefig('Data/plot_MW.pdf')
mannwhitney('LogP')

plt.figure(figsize=(5.5, 5.5))
sns.boxplot(x = 'class', y = 'LogP', data = df_2class)
plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('LogP', fontsize=14, fontweight='bold')
plt.savefig('Data/plot_LogP.pdf')
mannwhitney('NumHDonors')

plt.figure(figsize=(5.5, 5.5))
sns.boxplot(x = 'class', y = 'NumHDonors', data = df_2class)
plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('NumHDonors', fontsize=14, fontweight='bold')
plt.savefig('Data/plot_NumHDonors.pdf')

plt.figure(figsize=(5.5, 5.5))
sns.boxplot(x = 'class', y = 'NumHAcceptors', data = df_2class)
plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('NumHAcceptors', fontsize=14, fontweight='bold')
plt.savefig('Data/plot_NumHAcceptors.pdf')
mannwhitney('NumHAcceptors')

# Zip the resulting CSV and PDF files
# ! zip -r results.zip . -i *.csv *.pdf
