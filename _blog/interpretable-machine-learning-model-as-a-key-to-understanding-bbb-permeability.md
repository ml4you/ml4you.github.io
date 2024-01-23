---
title: 'Interpretable Machine Learning model as a key to understanding BBB permeability'
date: 2024-01-23
permalink: /blog/2024/01/interpretable-machine-learning-model-as-a-key-to-understanding-bbb-permeability
excerpt_separator: <!--more-->
toc: true
tags:
  - machine learning
  - cheminformatics
  - career
---

The blood-brain barrier (BBB) is a vital selective barrier in the central nervous system. Assessing the permeability of compounds across the BBB is crucial for drug development targeting the brain. While clinical experiments are accurate, they are time-consuming and costly. Computational methods offer an alternative for predicting BBB permeability.
<!--more-->

## 1. Downloading the dataset

In the first step of our tutorial, we initiate the process by downloading the essential dataset. This dataset, curated by Meng et al. in 2021, is a valuable resource comprising over 7000 compounds and 1613 chemical descriptors, calculated using Mordred fingerprints. To ensure a seamless experience, execute the provided command in your Python environment to obtain the dataset from the specified URL. This dataset serves as the foundation for our exploration into machine learning applications for predicting chemical drug properties, particularly focusing on aqueous water solubility and BBB permeability.

```
!wget https://github.com/theochem/B3DB/raw/87240af2b4e585d56f9681a6426af6b7f2940e96/B3DB/B3DB_classification_extended.tsv.gz
```

In the subsequent code snippet, we unpack the compressed dataset file for ease of access and analysis. The Python script utilizes the gzip library to decompress the file named "B3DB_classification_extended.tsv.gz." The decompressed file is then saved with the name "B3DB_classification_extended.tsv." This extraction process ensures that the dataset is in a readable format for further exploration and manipulation. After running this code, a confirmation message will be displayed, indicating that the extraction was successful. This step is crucial in preparing the dataset for subsequent machine learning analyses, allowing us to delve into predicting chemical drug properties with enhanced clarity and convenience.

```
import gzip
import shutil

input_file_path = "B3DB_classification_extended.tsv.gz"
output_file_path = "B3DB_classification_extended.tsv"

with gzip.open(input_file_path, 'rb') as f_in:
    with open(output_file_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

print(f"File '{input_file_path}' has been successfully extracted to '{output_file_path}'.")
```
In the following code snippet, we leverage the power of the pandas library to seamlessly import the extracted dataset into a structured and manipulable DataFrame. The dataset, now stored as "B3DB_classification_extended.tsv," is read into the variable 'df' using the 'read_csv' function from pandas. The 'sep' parameter is set to '\t' to indicate that the data is tab-separated, ensuring proper parsing. By displaying the DataFrame 'df,' we gain a preliminary glimpse into the dataset's structure and content. This step marks a pivotal moment as we transition from data acquisition to data exploration, setting the stage for in-depth analyses and insights into the chemical properties encapsulated within the dataset. With the dataset loaded into memory, we are ready to unleash the capabilities of pandas for comprehensive data exploration and preprocessing, paving the way for subsequent machine learning applications.
```
df = pd.read_csv("B3DB_classification_extended.tsv", sep='\t')
df
```
![image](https://github.com/yboulaamane/yboulaamane.github.io/assets/7014404/f5bbc971-d805-405b-91de-48dfd032de01)


## 2. Curating the dataset

Next, we use the 'dropna' method from pandas to efficiently handle missing values within our DataFrame ('df'). The 'axis=1' parameter signifies that the operation is applied along columns, effectively removing any columns containing NaN (Not a Number) values. This step is crucial for ensuring the cleanliness and completeness of our dataset, setting the foundation for accurate and robust machine learning analyses. By executing this line, we enhance the dataset's quality and prepare it for subsequent feature selection, model training, and predictions.
```
df = df.dropna(axis=1)
```
![image](https://github.com/yboulaamane/yboulaamane.github.io/assets/7014404/90e8705d-76be-4ddb-942c-ba7dcd71e5ee)

## 3. Labelling the dataset
In this code snippet, a new column named 'labels' is added to the DataFrame 'df.' The values in this column are determined based on the 'BBB+/BBB-' column. The 'apply' function, combined with a lambda function, assigns a binary label: 0 if 'BBB-' and 1 if 'BBB+'. This line of code is instrumental in preparing the dataset for a supervised machine learning task, where we aim to predict the binary outcome of blood-brain barrier (BBB) permeability. The 'labels' column now serves as the target variable, facilitating the training and evaluation of machine learning models for predicting BBB permeability.
```
df['labels'] = df['BBB+/BBB-'].apply(lambda x: 0 if x == 'BBB-' else 1)
```

## 4. Selecting chemical descriptors
In this code snippet, a subset of the DataFrame 'df' is extracted to form the 'features' DataFrame. The 'iloc' method is employed to select columns within a specific index range, denoted by 'abc_column_index' and 'mZagreb2_column_index.' This operation isolates the columns containing the features used for machine learning analysis. By creating the 'features' DataFrame, we focus on the relevant input variables necessary for training our machine learning models to predict blood-brain barrier (BBB) permeability. This step is pivotal in delineating the predictor variables from the target variable, facilitating streamlined model development and analysis.
```
features = df.iloc[:, 6:738]
```
Next, we extract the binary classification labels for blood-brain barrier (BBB) permeability from the "labels" column of the DataFrame 'df' and assigns them to the variable 'labels.'

```
labels = df["labels"]
```

## 5. Building the model
To build the ML model, we import essential modules from scikit-learn, a popular machine learning library in Python. The modules include RandomForestClassifier, svm, SVC, train_test_split for data splitting, cross_val_score, cross_val_predict for cross-validation, confusion_matrix for confusion matrix computation, roc_auc_score for ROC AUC score calculation, and classification_report for generating a classification report. These modules collectively provide a robust foundation for building, training, and evaluating machine learning models.

```
from sklearn.ensemble import RandomForestClassifier
from sklearn import svm
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score, cross_val_predict
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_auc_score
from sklearn.metrics import classification_report
```

In the next line of code, the dataset is split into training and testing sets using the 'train_test_split' function from scikit-learn. The 'features' and 'labels' variables represent the input features and target labels, respectively. The 'test_size' parameter is set to 0.3, allocating 30% of the data for testing. The 'random_state' ensures reproducibility, and 'shuffle=True' randomizes the data before splitting. The 'stratify=labels' parameter ensures that the class distribution is maintained in both the training and testing sets, which is crucial for balanced model training and evaluation.

```
X_train, X_test, y_train, y_test=train_test_split(features, labels, test_size=0.3, random_state=42, shuffle=True, stratify=labels)
```

Then, a RandomForestClassifier model is instantiated using scikit-learn. The model is then trained on the training data ('X_train' for features and 'y_train' for labels) using the 'fit' method. This marks a crucial step in the machine learning workflow, where the algorithm learns patterns from the training data to make predictions on new, unseen data. The variable 'rf' now holds the trained Random Forest Classifier ready for evaluation and prediction tasks.

```
rf = RandomForestClassifier()
rf = rf.fit(X_train, y_train)
```

Next, we employ the trained Random Forest Classifier ('rf') to make predictions on the testing data ('X_test'). The predicted values are stored in the variable 'y_pred,' representing the model's anticipated outcomes for the corresponding features in the testing set. This step allows us to assess the model's performance by comparing its predictions against the actual labels in the testing data.

## 6. Evaluating the model's performance

The trained Random Forest Classifier ('rf') is then employed to make predictions on the testing data ('X_test'). The predicted values are stored in the variable 'y_pred,' representing the model's anticipated outcomes for the corresponding features in the testing set. This step allows us to assess the model's performance by comparing its predictions against the actual labels in the testing data.

```
y_pred = rf.predict(X_test)
```
The next line prints a concise classification report, summarizing the performance metrics of the Random Forest Classifier on the test set.

```
print(classification_report(y_pred, y_test))
```
![image](https://github.com/yboulaamane/yboulaamane.github.io/assets/7014404/0b95329d-f1f8-4d30-8ad4-252a6bbb9337)

Afterwards, we calculate the ROC AUC score, a key performance metric for binary classification models like the Random Forest Classifier. The 'roc_auc_score' function from scikit-learn computes the area under the Receiver Operating Characteristic (ROC) curve, providing a single value to gauge the model's ability to distinguish between the two classes. The resulting score is printed as "ROC AUC Score."

```
roc_auc = roc_auc_score(y_test, y_pred)
print("ROC AUC Score:", roc_auc)
```

## 7. Calculating the most important features

The line of code below retrieves the feature importances calculated by the trained Random Forest Classifier. The attribute 'feature_importances_' provides insights into the contribution of each feature in making predictions. The resulting array contains importance scores corresponding to the features used in the model. Analyzing these scores can help identify the most influential features in predicting blood-brain barrier (BBB) permeability.
```
rf.feature_importances_
```
In the next code snippet, a new DataFrame named 'xfeatures' is created, consolidating information about the features and their respective importance values obtained from the Random Forest Classifier. The DataFrame consists of two columns: "features," representing the feature names extracted from the original dataset, and "Imp_values," containing the corresponding feature importances calculated by the model. This DataFrame provides a clear and structured summary, making it easy to analyze and interpret the significance of each feature in predicting blood-brain barrier (BBB) permeability.
```
xfeatures=pd.DataFrame({"features":features.columns, "Imp_values":rf.feature_importances_})

```
We sort the 'xfeatures' DataFrame in descending order based on the "Imp_values" column, providing a ranked view of feature importances. The resulting DataFrame showcases the features in decreasing order of importance, enabling a quick identification of the most influential factors in predicting blood-brain barrier (BBB) permeability according to the Random Forest Classifier.

```
xfeatures=xfeatures.sort_values("Imp_values", ascending=False)
xfeatures
```
![image](https://github.com/yboulaamane/yboulaamane.github.io/assets/7014404/fb2c60dc-41dc-4048-8d34-82f2615c8be0)


## 8. Interpreting the model using SHAP explainer

We'll use the 'pip' package manager to install the 'shap' library. The 'shap' library is often used for explaining the output of machine learning models, providing insights into the contribution of each feature to individual predictions. Once installed, the 'shap' library can be imported and utilized in the analysis.

```
!pip install shap
```
The 'shap' library is employed to create an explainer ('explain') for the trained Random Forest Classifier ('rf'). The explainer is then used to compute Shapley values ('shapvalues') for the features in the testing set ('X_test'). Shapley values offer insights into the impact of each feature on individual predictions, providing a valuable tool for interpreting and understanding the model's decision-making process.

```
import shap
explain=shap.Explainer(rf)
shapvalues=explain.shap_values(X_test)
```
Then we generate a summary plot using the 'shap' library, visualizing the Shapley values for each feature across the entire testing set ('X_test'). The plot provides a concise overview of feature importance and their impact on model predictions, aiding in the interpretation of the Random Forest Classifier's decision-making process.

```
shap.summary_plot(shapvalues,X_test)
```
![image](https://github.com/yboulaamane/yboulaamane.github.io/assets/7014404/d584b45b-c037-4594-9b2d-c8fcc5dce254)

We can also create a summary plot using the 'shap' library, specifically focusing on the Shapley values for the first class in the binary classification. The plot visualizes the impact of each feature on model predictions for this specific class across the testing set ('X_test'). This targeted summary aids in understanding the contributions of individual features to predictions related to the first class.
```
shap.summary_plot(shapvalues[0], X_test)
```
![image](https://github.com/yboulaamane/yboulaamane.github.io/assets/7014404/cd1822f7-0384-476a-81c6-8081bb822991)

Finally, we generate a dependence plot using the 'shap' library. The plot focuses on the most important feature "TopoPSA" and illustrates how its values impact the Shapley values and, consequently, the model predictions for the first class in the binary classification. Dependence plots are valuable for visualizing the relationship between a specific feature and the model output, enhancing interpretability.

```
shap.dependence_plot("TopoPSA", shapvalues[0], X_test)
```
![image](https://github.com/yboulaamane/yboulaamane.github.io/assets/7014404/9ef5df90-cfd1-48a6-9387-fa8a8bd9f70e)

## 9. Bottom line

In conclusion, the significance of explainable machine learning (ML) becomes evident in unraveling complex models' decision-making processes. In our analysis of blood-brain barrier (BBB) permeability prediction, employing the 'shap' library allowed us to interpret the Random Forest Classifier's output. Notably, our findings highlight "TopoPSA" as the most crucial feature influencing the model's predictions for BBB permeable compounds. This discovery aligns with existing literature, where molecular polar surface area (TopoPSA) has been consistently linked to BBB permeability. Our transparent and interpretable ML approach not only provides valuable insights into predictive features but also reinforces the importance of understanding model decisions for informed decision-making in drug development.






