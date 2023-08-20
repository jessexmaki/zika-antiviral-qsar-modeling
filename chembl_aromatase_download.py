import pandas as pd
from chembl_webresource_client.new_client import new_client
import os
import csv

# Initialize ChEMBL client
chembl_client = new_client 

# Search ChEMBL for aromatase targets
target = chembl_client.target
target_query = target.search('aromatase')
targets = pd.DataFrame.from_dict(target_query)

# Get first target ChEMBL ID  
selected_target = targets.target_chembl_id[0]

# Get bioactivity data for target filtered by IC50 assay type
activity = chembl_client.activity
res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
df = pd.DataFrame.from_dict(res)

# Get list of molecule Chembl IDs
chembl_ids = df['molecule_chembl_id'].tolist() 

# Lookup molecule names using ChEMBL API
molecule = chembl_client.molecule
names = []
for id in chembl_ids:
  result = molecule.get(id)
  names.append(result['molecule_properties']['full_molformula'])

# Add molecule names to dataframe
df['name'] = names

# Create output directory if it doesn't exist
if not os.path.exists('Data'):
  os.makedirs('Data')

# Write dataframe to CSV  
try:
  with open('Data/bioactivity_data.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    df.to_csv(f, index=False)
except Exception as e: 
  print("Error saving CSV:")
  print(e)