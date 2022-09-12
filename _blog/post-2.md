---
title: 'Chemical databases for machine learning in drug discovery'
date: 2022-09-12
permalink: /blog/post-2
excerpt_separator: <!--more-->
toc: true
tags:
  - machine learning
  - cheminformatics
  - career
---

The ever-increasing bioactivity data that are produced nowadays allow exhaustive data mining and knowledge discovery approaches that change chemical biology research. 
A wealth of cheminformatics tools, web services, and applications therefore exists that supports a careful evaluation and analysis of experimental data to draw conclusions that can influence the further development of chemical probes and potential lead structures.
<!--more-->

## 1.	ChEMBL database
ChEMBL is a manually curated database of bioactive molecules with drug-like properties. It brings together chemical and bioactivity data to aid the translation of information into effective new drugs.
Bioactivity data is reported in Ki, Kd, IC50, % of inhibition and EC50. Data can be filtered and analyzed to develop compound screening libraries for lead identification during the drug discovery process. 
The availability of curated bioactivity data on small drug-like molecules opens new opportunities for data-driven drug discovery and makes it possible to apply machine learning methodologies to pharmaceutical research field.

 
Figure 1: Web interface of ChEMBL database (https://www.ebi.ac.uk/chembl/)

## 2. BindingDB 
BindingDB is a public database of experimental binding affinities of interactions between macromolecules and small molecules. It contains more than 1.2 million binding data points for more than 6,400 protein targets and around 550,000 small molecules. Furthermore, affinities for protein-protein, protein-peptide as well as host-guest interactions are provided. 
BindingDB provides a wide range of searches, possible queries, tools and datasets.
BindingDB also provides a virtual compound screening tool with which the user has the possibility to screen an external dataset of compounds for similar bioactivity. 
BindingDB also provides a virtual compound screening tool with which the user has the possibility to screen an external dataset of compounds for similar bioactivity. Access is given by a download option (SD-file, tab-separated value (TSV) or Oracle dump) or by programmatic access (RESTful API, structured URLs, or KNIME).

### 2.1.	Find My Compound’s Target tool
Find My Compound’s Target is a tool integrated in BindingDB that aims to predict the target of a small molecule of interest or possible off-targets. The query compound is first compared to other compounds in the database and targets of these compounds are selected if the compounds’ similarity is above the chosen cut-off and the affinity is respectively above a certain threshold.
### 2.2.	Find Compounds for My Target
Find Compounds for My Target is another tool that tries to find compounds for a specific target. It features advanced search and extraction of curated information from various data sources such as ChEMBL database, PubChem bioassays, D3R, etc. 

  
Figure 2: Web interface of BindingDB (https://www.bindingdb.org/)


## 3.	Pubchem Bioassays Database
PubChem is an American database of chemical molecules managed by the National Center for Biotechnology Information (NCBI), branch of the National Library of Medicine of the United States under the authority of the National Institutes of Health (NIH).
PubChem lists several million compounds by putting a large amount of data of various kinds online for each substance free of charge: chemical, biochemical, pharmacological, production, toxicological, environmental, etc.
PubChem BioAssay stores the activity data of small molecules or RNAi and contains curated parts of ChEMBL for which flags for active or inactive compounds are assigned depending on whether the IC50, EC50 or Ki is above 50 μM or not. The data can be accessed and analyzed via a broad range of provided web services and tools (Figure1). Besides using a name, a smiles code can be used or a structure can be drawn to search for an identical molecule, a similar molecule or a substructure. This
leads to information about bioassay results or substance descriptions. Variations of assay results can be analyzed using detailed description of the performed experiments. The data is additionally clustered, e.g., according to the protein or gene target, the type of assay (e.g., cell-based, protein-protein interaction), an assay project or more complex kinds of relationships like target similarity or common active compounds. 

### 3.1.	PubChemPy
PubChemPy offers a way to use and interact with PubChem database directly with Python. It allows chemical searches by name, substructure and similarity, chemical standardization, conversion between chemical file formats, depiction and retrieval of chemical properties. For more information on installing and using PubChemPy package visit the official website (https://pubchempy.readthedocs.io/en/latest/index.html).

 
Figure 3: PubChem web interface (https://pubchem.ncbi.nlm.nih.gov/)
## 4.	PDBbind Database
The PDBbind database is a comprehensive collection of experimentally measured binding affinity data (Kd, Ki, and IC50) for the protein-ligand complexes deposited in the Protein Data Bank (PDB). It thus provides a link between energetic and structural information of protein-ligand complexes, which is of great value to various studies on molecular recognition occurred in biological systems.
The basic information of each complex in PDBbind is totally open for browsing. Users are however required to register for access under a license agreement in order to utilize the full functions provided on this web site or to download the contents of PDBbind. The registration is free of charge to all academic and industrial users.
 
Figure 4: PDBbind web interface (http://www.pdbbind.org.cn/index.php) 




## 5.	BRENDA Enzyme Database
BRENDA is an enzyme database. It is maintained and developed by the Institute of Biochemistry of the University of Cologne. Data on enzyme functions are taken directly from the primary literature. The database covers 40 entries, with information about enzyme nomenclature, reactions, specificity, structure, method of isolation or preparation, references in scientific literature and cross-references for the sequence or 3D structure.
The database is accessible free of charge for academic and non-profit uses, commercial uses need to acquire a license. To use the database, it is necessary to register by email. The database can be searched by EC nomenclature, enzyme name, organism or an advanced search combining these entries.
 
Figure 5: Brenda Enzymes Database web interface (https://www.brenda-enzymes.org/)

## Conclusion
This article has sought to provide an overview of freely available chemical databases that can be readily used for machine learning approaches that support the major research trends in drug discovery. The ever-increasing amount of data and the improvement of analytical tools hold the potential to transform the drug development, leading to new treatments, improved patient outcomes, and lower costs.
