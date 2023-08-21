# Bioinformatics-Project

This project utilizes bioactivity data on Zika virus inhibitors from the ChEMBL database. The raw bioactivity data is preprocessed to remove duplicates, filter out missing values, and classify compounds into active/inactive/intermediate categories based on potency cutoffs. Key physicochemical properties are calculated including pIC50, molecular weight, LogP, hydrogen bond donors/acceptors. The curated dataset is analyzed using statistical tests and visualizations to compare property distributions between active and inactive compounds.

Molecular descriptors are computed from SMILES strings using the PaDEL software package. Descriptors with low variance are filtered out to derive a final set of informative descriptors for model building.

A Random Forest regression model is developed to predict the pIC50 antiviral activity from the descriptors. 

The trained regression model can be applied to new compounds by taking their SMILES string, computing descriptors using PaDEL, and making a prediction using the saved Random Forest model. This allows rapid prediction of new Zika inhibitors without needing to synthesize and test each compound experimentally.

Overall, this workflow demonstrates a typical QSAR modeling approach leveraging public bioactivity data, physicochemical properties, molecular descriptors, machine learning algorithms, and rigorous validation to derive predictive models for accelerated antiviral drug discovery. The code encapsulates data cleaning, analysis, model development, and application to new compounds.

-----
***LINKS***
PaDel Software: https://github.com/dataprofessor/bioinformatics/
