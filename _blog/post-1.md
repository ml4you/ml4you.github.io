---
title: 'How to use PubChem database for machine learning in drug discovery'
date: 2022-09-07
permalink: /blog/post-1
excerpt_separator: <!--more-->
toc: false
tags:
  - machine learning
  - cheminformatics
  - career
---

**What is PubChem database**  
PubChem is an American database of chemical molecules managed by the National Center for Biotechnology Information (NCBI), branch of the National Library of Medicine of the United States under the authority of the National Institutes of Health (NIH).
PubChem lists several million compounds by putting a large amount of data of various kinds online for each substance free of charge: chemical, biochemical, pharmacological, production, toxicological, environmental, etc.  
**What is PubChemPy**  
PubChemPy offers a way to use and interact with PubChem database directly with Python. It allows chemical searches by name, substructure and similarity, chemical standardization, conversion between chemical file formats, depiction and retrieval of chemical properties.
**Installing PubChemPy**  
The easiest and recommended way to install is using pip by simply running the code below:

```
!pip install pubchempy
```


The first step is to search the ChEMBL database using keywords of a target protein of interest, it is possible to run a search using other keywords related to diseases, compounds or assays. In this tutorial, we are going to search for Acetylcholinesterase as illustrated in figure 2.

![Figure1](https://raw.githubusercontent.com/yboulaamane/yboulaamane.github.io/984c2b8dff1a546b4fd9ad6d2cdb6e57e572851f/_blog/1post-1.png)
Figure 1: ChEMBL search result example

Notice that our search resulted on 24 targets, it is important to choose the right protein for the right organism of the study of interest. For this example, we are interested in human Acetylcholinesterase corresponding to the ID: CHEMBL220. 
After clicking on the target, we will be sent to another page containing all the data concerning the selected target such as: name and classification, drugs and clinical candidates, activity charts, etc. 
Scroll down to activity charts and notice the pie chart on the left concerning all the associated bioactivity data compiled from the literature and their distribution according to the activity type.  
 
 
![Figure1](https://raw.githubusercontent.com/yboulaamane/yboulaamane.github.io/984c2b8dff1a546b4fd9ad6d2cdb6e57e572851f/_blog/1post-2.png)
Figure 2: Activity charts and distribution of activity types of the selected target, CHEMBL220

Upon observation of the activity chart, we can quickly determine which activity type is the most reported in the literature, in this case it refers to half-maximal inhibitory concentrations (IC50) which have been reported 8205 times. 
Once we click on the desired activity type, we can download the entire dataset in CSV or TSV Format containing various informations such as ChEMBL ID for each compound, SMILES, Standard Type and Standard Value referring to the activity type and value respectively.
Data curation
Note that it is necessary to remove any unwanted data before proceeding with data curation. In this case, we are only interested in the compound’s IDs, Smiles, Standard Type and Standard Value. It is possible to perform this task with any CSV reader such as Google Sheets or Microsoft Excel. 
Once we have performed the primary cleaning on our data, we can import it on Google Colab or Jupyter Notebook using the code below:


**Import necessary libraries**
```
import pandas as pd
```
**Read the dataset**
```
x=pd.read_csv('ache.csv')
```
**Display the dataset**
```
x
```

Output:

![Figure3](https://raw.githubusercontent.com/yboulaamane/yboulaamane.github.io/984c2b8dff1a546b4fd9ad6d2cdb6e57e572851f/_blog/1post-3.png)

 


**Remove duplicate compounds**

When dealing with a large dataset of compounds, it is very likely to find a great deal of duplicates with different reported activities due to different conditions of each laboratory. However, it is possible to deal with this issue by averaging all reported activity by calculating their mean values using the code below:

```
x["mean_value"]=x[['Molecule ChEMBL ID', 'Smiles','Standard Type','Standard Value']].groupby(['Molecule ChEMBL ID'])['Standard Value'].transform("mean")
```
The next step is to merge all the duplicate compounds into one, for this reason we can use the code below to remove all duplicates while keeping only the first one.
```
x=x.drop_duplicates("Molecule ChEMBL ID", keep="first")
```

Output:

![Figure4](https://raw.githubusercontent.com/yboulaamane/yboulaamane.github.io/984c2b8dff1a546b4fd9ad6d2cdb6e57e572851f/_blog/1post-4.png)
 
It is possible to find some compounds on the dataset with no available activity, notice on the sheet above that some activities are marked with “NaN” which stands for “Not a Number” in computer science, it is therefore necessary to remove them before proceeding. We can simply run the code below:

```
x=x.dropna()
x
```

Output:

![Figure5](https://raw.githubusercontent.com/yboulaamane/yboulaamane.github.io/984c2b8dff1a546b4fd9ad6d2cdb6e57e572851f/_blog/1post-5.png)
 
**Data classification**

Once we have curated our data, now it is possible to classify compounds in order to apply it for machine learning classification models. For this reason, we need to define an activity cutoff to define our active and inactive compounds. In the case of enzyme inhibition, the literature indicates that most potent enzyme inhibitors have activities in the nanomolar range. For this reason, we can proceed by setting a threshold of 1000 nM corresponding to 1 μM or lower for defining our active compounds. 

Run the code below to create a variable with all active compounds:

```
active=x.loc[x['mean_value']<=1000]
```

We do the same for inactive compounds by setting a cutoff of 10 000 nM (10 μM) or higher using the code below:

```
inactive=x.loc[x['mean_value']>10000]
```

**Data labelling**

Now that we have defined our active and inactive compounds, it is necessary to label the data in order to combine the entire dataset. We will simply refer to active compounds as “1” and inactive compounds as “0”.
Run the code below:

```
active["Class"]=1
inactive["Class"]=0
```


Now we can proceed to combining the entire dataset.

Run the code below:

```
combined=pd.concat([active,inactive],axis=0)
combined
```


Output:
 
![Figure6](https://raw.githubusercontent.com/yboulaamane/yboulaamane.github.io/984c2b8dff1a546b4fd9ad6d2cdb6e57e572851f/_blog/1post-6.png)

Finally, we can save our dataset for further use.
Run the code below:

```
combined.to_csv("ache_labelled.csv", index=None)
```


**Conclusion**

This article’s aim was to demonstrate an alternative way to retrieve bioactivity data from ChEMBL without using code as opposed to the previous article. Furthermore, data curation and data classification was covered in detail as it is a necessary step and can highly impact the performance of machine learning models. If you found this article useful, follow the blog for more tutorials in the future.

