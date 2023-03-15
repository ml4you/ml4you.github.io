---
title: 'Using PaDELPy to generate molecular fingerprints for machine learning-based QSAR'
date: 2023-03-15
permalink: /blog/2023/03/using-padelpy-generate-molecular-fingerprints-machine-learning-based-qsar
excerpt_separator: <!--more-->
toc: true
tags:
  - machine learning
  - drug discovery
  - qsar
---
PaDELPy is a Python library that wraps the PaDEL-Descriptor molecular descriptor calculation software and can be used to build scientific machine learning models. Machine learning models are created by training an algorithm to recognize patterns in data and can be either supervised or unsupervised. There are many machine learning algorithms, such as classification and regression, and they can be implemented in languages such as Python or R. The efficiency and accuracy of both the algorithm and the model can be analyzed and calculated.
<!--more-->

![Figure1](https://user-images.githubusercontent.com/7014404/225259643-df0568cd-1cfe-4395-aa7e-980902108f25.png)

## What Is PaDELPy?
PaDELPy is a Python library that simplifies the use of the PaDEL-Descriptor molecular descriptor calculation software to calculate molecular fingerprints for scientific machine learning models. It is a wrapper for the Java-based PaDEL-Descriptor, eliminating the need for running a Java file and making implementation quicker. 
## Getting Started with the Code
In this article, we will try to create a scientific machine learning model to predict molecular activity in the HCV Drug dataset using Random Forest. The code is based on the PaDELPy library's creators and the dataset can be downloaded from a provided [link](https://github.com/chaninlab/hcvpred/blob/master/HCV_NS5B_Curated.csv).   
To install the PaDELPy library, execute the following line of code:
```
#installing the library
!pip install padelpy
```
Next, we will load and configure our calculator model using the required PaDELPy files, which are only available in XML format.
```
#Downloading the XML data files
!wget https://github.com/dataprofessor/padel/raw/main/fingerprints_xml.zip
!unzip fingerprints_xml.zip
#listing and sorting the downloaded files
import glob
xml_files = glob.glob("*.xml")
xml_files.sort()
xml_files
```
Output:
```
['AtomPairs2DFingerprintCount.xml',
 'AtomPairs2DFingerprinter.xml',
 'EStateFingerprinter.xml',
 'ExtendedFingerprinter.xml',
 'Fingerprinter.xml',
 'GraphOnlyFingerprinter.xml',
 'KlekotaRothFingerprintCount.xml',
 'KlekotaRothFingerprinter.xml',
 'MACCSFingerprinter.xml',
 'PubchemFingerprinter.xml',
 'SubstructureFingerprintCount.xml',
 'SubstructureFingerprinter.xml']
#Creating a list of present files
FP_list = ['AtomPairs2DCount',
 'AtomPairs2D',
 'EState',
 'CDKextended',
 'CDK',
 'CDKgraphonly',
 'KlekotaRothCount',
 'KlekotaRoth',
 'MACCS',
 'PubChem',
 'SubstructureCount',
 'Substructure']
```
We will now create a data dictionary that includes all the loaded and available data files to obtain a key-value pair.
```
#Creating Data Dictionary
fp = dict(zip(FP_list, xml_files))
fp
```
Output:
```
{'AtomPairs2D': 'AtomPairs2DFingerprinter.xml',
 'AtomPairs2DCount': 'AtomPairs2DFingerprintCount.xml',
 'CDK': 'Fingerprinter.xml',
 'CDKextended': 'ExtendedFingerprinter.xml',
 'CDKgraphonly': 'GraphOnlyFingerprinter.xml',
 'EState': 'EStateFingerprinter.xml',
 'KlekotaRoth': 'KlekotaRothFingerprinter.xml',
 'KlekotaRothCount': 'KlekotaRothFingerprintCount.xml',
 'MACCS': 'MACCSFingerprinter.xml',
 'PubChem': 'PubchemFingerprinter.xml',
 'Substructure': 'SubstructureFingerprinter.xml',
 'SubstructureCount': 'SubstructureFingerprintCount.xml'}
```
After setting up all the required PaDELPy files for calculation, our next step is to load the dataset that needs to be processed.
```
#Loading the dataset
import pandas as pd
df = pd.read_csv('https://raw.githubusercontent.com/dataprofessor/data/master/HCV_NS5B_Curated.csv')
```
```
#Loading data head
df.head()
```
![Figure2](https://user-images.githubusercontent.com/7014404/225259647-e6996ee3-f4f2-45a9-93ad-6d8299c304cd.png)
```
#Loading data tail
df.tail(2)
```
![Figure3](https://user-images.githubusercontent.com/7014404/225259650-4b67958c-cd8f-4d5a-b630-6cdc42e548a2.png)

In order to calculate the molecular descriptor using PaDEL, it is necessary to prepare the data by concatenating the two relevant columns from the dataset. This concatenated data will serve as input to our model.
```
#Concatenating necessary columns
df2 = pd.concat( [df['CANONICAL_SMILES'],df['CMPD_CHEMBLID']], axis=1 )
df2.to_csv('molecule.smi', sep='\t', index=False, header=False)
df2
```
![Figure4](https://user-images.githubusercontent.com/7014404/225259651-e4e3bae9-a8b9-4238-88e2-c02a16b7e1c6.png)
To compute the molecular fingerprint using PaDEL, we can choose from 12 available fingerprint types. To calculate all 12 types, we need to modify the input argument for descriptor types to any of the options in the fp dictionary variable.
```
#listing the dictionary pairs
fp
```
Output:
```
{'AtomPairs2D': 'AtomPairs2DFingerprinter.xml',
 'AtomPairs2DCount': 'AtomPairs2DFingerprintCount.xml',
 'CDK': 'Fingerprinter.xml',
 'CDKextended': 'ExtendedFingerprinter.xml',
 'CDKgraphonly': 'GraphOnlyFingerprinter.xml',
 'EState': 'EStateFingerprinter.xml',
 'KlekotaRoth': 'KlekotaRothFingerprinter.xml',
 'KlekotaRothCount': 'KlekotaRothFingerprintCount.xml',
 'MACCS': 'MACCSFingerprinter.xml',
 'PubChem': 'PubchemFingerprinter.xml',
 'Substructure': 'SubstructureFingerprinter.xml',
 'SubstructureCount': 'SubstructureFingerprintCount.xml'}
```
To calculate the molecular fingerprint, we need to load the required file. In this case, we will be using PubChem.
```
#Importing PubChem
fp['PubChem']
```
Setting up the module to calculate the molecular fingerprint.
```
#Setting the fingerprint module
from padelpy import padeldescriptor
fingerprint = 'Substructure'
fingerprint_output_file = ''.join([fingerprint,'.csv']) #Substructure.csv
fingerprint_descriptortypes = fp[fingerprint]
padeldescriptor(mol_dir='molecule.smi', 
                d_file=fingerprint_output_file, #'Substructure.csv'
                #descriptortypes='SubstructureFingerprint.xml', 
                descriptortypes= fingerprint_descriptortypes,
                detectaromaticity=True,
                standardizenitro=True,
                standardizetautomers=True,
                threads=2,
                removesalt=True,
                log=True,
                fingerprints=True)

```
Displaying the calculated fingerprints.
```
descriptors = pd.read_csv(fingerprint_output_file)
descriptors
```
![Figure5](https://user-images.githubusercontent.com/7014404/225259654-91453085-9eea-4129-a763-6fe458ce7634.png)

Next, we will utilize the processed data and Random Forest algorithm to create a classification model.
```
X = descriptors.drop('Name', axis=1)
y = df['Activity'] #feature being predicted

#removing the low variance features
from sklearn.feature_selection import VarianceThreshold
 
def remove_low_variance(input_data, threshold=0.1):
    selection = VarianceThreshold(threshold)
    selection.fit(input_data)
    return input_data[input_data.columns[selection.get_support(indices=True)]]
 
X = remove_low_variance(X, threshold=0.1)
X
```
![Figure6](https://user-images.githubusercontent.com/7014404/225259660-deff3ff9-cd5d-484e-b7e4-db47bc8c7dcf.png)

The processed data appears to be well-organized and informative, highlighting the most effective drugs. Therefore, we can proceed to generate predictions using our model.
```
#Splitting into Train And Test
from sklearn.model_selection import train_test_split
 
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

#Printing Shape
X_train.shape, X_test.shape

((462, 18), (116, 18))

#Implementing Random Forest
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import matthews_corrcoef
 
model = RandomForestClassifier(n_estimators=500, random_state=42)
model.fit(X_train, y_train)
```
Output:
```
RandomForestClassifier(bootstrap=True, ccp_alpha=0.0, class_weight=None,
                       criterion='gini', max_depth=None, max_features='auto',
                       max_leaf_nodes=None, max_samples=None,
                       min_impurity_decrease=0.0, min_impurity_split=None,
                       min_samples_leaf=1, min_samples_split=2,
                       min_weight_fraction_leaf=0.0, n_estimators=500,
                       n_jobs=None, oob_score=False, random_state=42, verbose=0,
                       warm_start=False)
```
We can make predictions on the drug molecular activity from the created model by passing our test data through the model's predict function. This will output the predicted molecular activity of the drugs in the test set. Here's an example code snippet:
```
y_train_pred = model.predict(X_train)
y_test_pred = model.predict(X_test)
```
Calculating performance metrics of train split  using matthews correlation coefficient,
```
mcc_train = matthews_corrcoef(y_train, y_train_pred)
mcc_train
```
Output:
0.833162800019916

Calculating performance metrics of test split using matthews correlation coefficient,
```
mcc_test = matthews_corrcoef(y_test, y_test_pred)
mcc_test
```
Output:
0.5580628933757674
```
#performing cross validation
from sklearn.model_selection import cross_val_score
rf = RandomForestClassifier(n_estimators=500, random_state=42)
cv_scores = cross_val_score(rf, X_train, y_train, cv=5)
cv_scores
```
Output:
array([0.83870968, 0.80645161, 0.86956522, 0.86956522, 0.81521739])
```
#calcutating mean from the five fold 
mcc_cv = cv_scores.mean()
mcc_cv
```
Output:
0.8399018232819074
```
#implementing metric test in a single dataframe
model_name = pd.Series(['Random forest'], name='Name')
mcc_train_series = pd.Series(mcc_train, name='MCC_train')
mcc_cv_series = pd.Series(mcc_cv, name='MCC_cv')
mcc_test_series = pd.Series(mcc_test, name='MCC_test')
performance_metrics = pd.concat([model_name, mcc_train_series, mcc_cv_series, mcc_test_series], axis=1)
performance_metrics
```
![Figure7](https://user-images.githubusercontent.com/7014404/225259664-4c0cabf5-a6e1-4d2f-9533-493881c0780b.png)
Based on the performance metrics of the Random Forest model used for predicting molecular drug activity, it appears to be effective on the current dataset. However, it is possible to test the performance of other algorithms as well.

## Conclusion

This article provided an explanation of the meaning and importance of machine learning models. Additionally, we utilized the PaDELPy library to create a scientific learning model which employs the PaDEL Descriptor to compute molecular fingerprints. Finally, we evaluated the model's performance metrics.




