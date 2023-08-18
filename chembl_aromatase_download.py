import pandas as pd  
from chembl_webresource_client.new_client import new_client  

# Searching for 'aromatase' targets in ChEMBL
target = new_client.target
target_query = target.search('aromatase')
targets = pd.DataFrame.from_dict(target_query)

# Selecting the first ChEMBL ID
selected_target = targets.target_chembl_id[0]

# Fetching activity data for the selected target with type 'IC50'
activity = new_client.activity
res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
df = pd.DataFrame.from_dict(res)

# Saving data to CSV
df.to_csv('Data/bioactivity_data_raw.csv', index=False)

