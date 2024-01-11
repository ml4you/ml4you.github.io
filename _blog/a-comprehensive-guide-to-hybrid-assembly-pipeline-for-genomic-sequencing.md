---
title: 'A Comprehensive Guide to Hybrid Assembly Pipeline for Genomic Sequencing'
date: 2024-01-11
permalink: /blog/2024/01/a-comprehensive-guide-to-hybrid-assembly-pipeline-for-genomic-sequencing
excerpt_separator: <!--more-->
toc: true
tags:
  - genomics
  - bioinformatics
---

Microbes, life's unseen workhorses, hold immense potential for bioremediation, medicine, and understanding our planet. Yet, their intricate workings remain largely a mystery. This is where microbial genomics steps in, offering a powerful tool to decode their genetic language.

<!--more-->

## 1. Introduction
Genomic sequencing has revolutionized our understanding of microbial diversity and function. In this guide, we'll walk through a hybrid assembly pipeline, combining long reads (ONT) and short reads (Illumina), using various tools for quality control, assembly, polishing, and assessment.


## 2. Requirements

**1. Using Anaconda for Package Management**

Anaconda is a powerful package manager and environment manager that simplifies the process of installing, managing, and updating software packages. It is particularly useful for managing bioinformatics tools and their dependencies.

Download and install Anaconda from the [official website](https://www.anaconda.com/download).  

Create a new conda environment for your project:
```
conda create --name myenv
conda activate myenv
```

Before diving into the pipeline, ensure you have the necessary tools installed:

**2. Quality Control:**
    - NanoPlot: [GitHub](https://github.com/wdecoster/NanoPlot)  
    - FastQC: [GitHub](https://github.com/s-andrews/FastQC)  
    - Filtlong: [GitHub](https://github.com/rrwick/Filtlong#installation)  

**3. Assembly:**
    - Flye: [GitHub](https://github.com/fenderglass/Flye/blob/flye/docs/INSTALL.md)  

**4. Polishing:**
    - Medaka: [GitHub](https://github.com/nanoporetech/medaka)  
    - BWA: [GitHub](https://github.com/lh3/bwa)  
    - PolyPolish: [Bioconda](https://bioconda.github.io/recipes/polypolish/README.html)  

**5. Assembly Assessment:**
    - BUSCO: [User Guide](https://busco.ezlab.org/busco_userguide.html#conda-package)  
    - QUAST: [Website](https://quast.sourceforge.net/install.html)  

## 3. Hybrid Assembly Pipeline

### Step 1: Quality Control

#### For Long Reads (LR - ONT):
```
NanoPlot --fastq LR_input.fastq --N50 --verbose --outdir 1-NanoPlot_LR_raw/ -t 8
```
#### For Short Reads (SR - Illumina):

### Step 2: Filter Long Reads
```
filtlong -1 SR_input_1.fastq -2 SR_input_2.fastq --min_length 1000 --keep_percent 90 LR_input.fastq > LR_filtered.fastq
NanoPlot --fastq LR_filtered.fastq --N50 --verbose --outdir 2-NanoPlot_LR_filtered/ -t 8
```

### Step 3: Long Reads Assembly using Flye
```
python3 Flye/bin/flye --nano-corr LR_filtered.fastq --out-dir Flye/ --threads 8 --scaffold -g 6m
```

### Step 4: First Polishing using Medaka
```
conda activate medaka
medaka_consensus -i LR_filtered.fastq -d Flye/assembly.fasta -o Polish1/ -m r941_min_fast_g303 -t 8
```

### Step 5: Second Polishing using Polipolish
```
mkdir Polish2
bwa index Polish1/consensus.fasta
bwa mem -t 8 -a Polish1/consensus.fasta SR_input_1.fastq > Polish2/alignments_1.sam
bwa mem -t 8 -a Polish1/consensus.fasta SR_input_2.fastq > Polish2/alignments_2.sam
polypolish_insert_filter.py --in1 Polish2/alignments_1.sam --in2 Polish2/alignments_2.sam --out1 Polish2/filtered_1.sam --out2 Polish2/filtered_2.sam
polypolish Polish1/consensus.fasta Polish2/filtered_1.sam Polish2/filtered_2.sam > final_assembly.fasta
```

### Step 6: Assembly Quality Assessment
```
conda activate busco
busco -i final_assembly.fasta -l bacteria_odb10 -m genome -o busco_final_assembly
quast.py -o quast_final_assembly -t 8 final_assembly.fasta
```

## 4. Annotation for 16S rRNA using Prokka and RAST

For further annotation, consider using tools like Prokka and RAST to identify and annotate features such as 16S rRNA genes. These tools can provide additional insights into the functional elements of your assembled genome.

### Prokka 

Prokka is a versatile tool for bacterial genome annotation. It predicts protein-coding genes, rRNAs, tRNAs, and other features.

```
conda install -c conda-forge -c bioconda prokka
prokka --outdir prokka_annotation --prefix final_assembly final_assembly.fasta
```

Prokka will generate annotation files in the specified directory, providing detailed information about the genomic features.

### RAST

Visit the RAST website (http://rast.nmpdr.org/) to submit your genome for annotation. RAST offers a user-friendly interface for functional annotation of bacterial and archaeal genomes.

## 5. Bottom line

This hybrid assembly pipeline, coupled with annotation tools like Prokka and RAST, empowers researchers to unravel the intricacies of microbial genomes. The combination of long and short reads, along with rigorous quality control and assessment, ensures a reliable and accurate representation of genomic information.

