import pandas as pd
import subprocess
import pickle

# Molecular descriptor calculator
def desc_calc():
    bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_output.csv"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

# Model building
def build_model(input_data):
    # Reads in saved regression model
    load_model = pickle.load(open('aromatase_model.pkl', 'rb'))
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    return prediction

# Load data from input_smiles.txt
load_data = pd.read_table('Testing/input_smiles.txt', sep=' ', header=None)
print("Loaded input data:")
print(load_data.head())

# Save smiles to a file for descriptor calculation
load_data.iloc[:, 0].to_csv('molecule.smi', sep='\t', header=False, index=False)

print("\nCalculating descriptors...")
desc_calc()

# Load the calculated descriptors
desc = pd.read_csv('descriptors_output.csv')

# Subset descriptors (you may need to adjust this part if your previous model used a different set of descriptors)
Xlist = list(pd.read_csv('descriptor_list.csv').columns)
desc_subset = desc[Xlist]

# Apply trained model to get predictions
predictions = build_model(desc_subset)

# Ensure we're only assigning the number of predictions that match the number of rows in load_data
load_data['pIC50'] = predictions[:len(load_data)]

# Save the data with the predictions
load_data.to_csv('output_with_predictions.csv', index=False)

print("Predictions saved to output_with_predictions.csv.")
