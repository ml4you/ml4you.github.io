---
title: 'Supervised vs. unsupervised methods in machine learning'
date: 2022-10-16
permalink: /blog/2022/10/supervised-vs-unsupervised-methods-machine-learning
excerpt_separator: <!--more-->
toc: true
tags:
  - machine learning
  - cheminformatics
  - career
---

The increasing volume of biomedical data in chemistry and life sciences requires development of new methods and approaches for their analysis. New approaches have proved to show improvement and accelerate the joint drug discovery and development processes. 
<!--more-->
The accumulation of large datasets allows for better integration of machine learning and artificial intelligence to build and integrate more accurate models for predicting the bioactivity and the pharmacokinetics of new drugs in the pharmaceutical field. 
Machine learning approaches are divided into three broad categories, which correspond to learning patterns, depending on the nature of the "signal" or "feedback" available to the learning system.

## 1. Supervised learning

Machine learning models are supervised by loading them with knowledge so that we can have it predict future instances. Teaching the model requires training it with some data from a labeled dataset. 

| ![Figure1](https://user-images.githubusercontent.com/7014404/225256631-3a7927f8-0f0e-4c3a-9220-70bc0496e4db.png) |
|:--:|
| <b>Figure 1: Example of a chemical dataset viewed with Pandas.</b> |

The names up here which are called: molecule_chembl_id and smiles are called attributes. Other names such as standard_value represents the numerical values for each sample, whereas the class name represents a categorical value which can be either 1 (active) or 0 (inactive). The columns are called features which include the data. If we plot this data, and look at a single data point on a plot, it'll have all of these attributes that would make a row on this chart also referred to as an observation. Looking directly at the value of the data, you can have two kinds. The first is numerical, when dealing with machine learning, the most commonly used data is numeric. The second is categorical, that is its non-numeric because it contains characters rather than numbers. In this case, it's categorical because this dataset is made for classification.

| ![Figure2](https://user-images.githubusercontent.com/7014404/225256633-ae720fd3-3c21-4cc5-bda4-4c4e8bdea928.png) |
|:--:|
| <b>Figure 2: Supervised learning.</b> |

There are two types of supervised learning techniques. They are classification, and regression. Classification is the process of predicting a discrete class label, or category. Wheras, regression is the process of predicting a continuous value as opposed to predicting a categorical value in classification.


## 2. Unsupervised learning

Unsupervised models are exactly what they sound like, the models are left to work on their own to discover information that may not be visible to the human eye. It means, the unsupervised algorithm trains on the dataset, and draws conclusions on unlabeled data. Generally speaking, unsupervised learning has more difficult algorithms than supervised learning since we know little to no information about the data, or the outcomes that are to be expected. Dimension reduction, density estimation, market basket analysis, and clustering are the most widely used unsupervised machine learning techniques. Dimensionality reduction, and/or feature selection, play a large role in this by reducing redundant features to make the classification easier.

| ![Figure3](https://user-images.githubusercontent.com/7014404/225256636-fbbe6ac3-8726-46ff-aad4-9025c25bd75a.png) |
|:--:|
| <b>Figure 3: Unsupervised learning tasks. Image by Dmytro Nikolaiev (medium.com/@andimid).</b> |

Market basket analysis, on the other hand, is a modeling technique based upon the theory that if you buy a certain group of items, you're more likely to buy another group of items. 
Density estimation is a very simple concept that is mostly used to explore the data to find some structure within it. 
Finally, clustering is consired to be one of the most popular unsupervised machine learning techniques used for grouping data points, or objects that are somehow similar. 
Clustering analysis has many applications in different domains, whether it be a bank's desire to segment his customers based on certain characteristics, or helping an individual to organize in-group his, or her favorite types of music. Broadly though, clustering is used mostly for discovering structure, summarization, and anamoly detection. 

## Bottom line

The biggest difference between supervised and unsupervised learning is that supervised learning deals with labeled data while unsupervised learning deals with unlabeled data. Supervised models employ machine learning algorithms for classification and regression, whereas in unsupervised learning, we have methods such as clustering. In comparison to supervised learning, unsupervised learning has fewer models and fewer evaludation methods that can be used to ensure the outcome of the model is accurate. As such, unsupervised learning creates a less controllable environment as the machine is creating outcomes for us.
