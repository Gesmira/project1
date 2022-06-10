---
date: "2016-04-27T00:00:00Z"
external_link: ""
image:
  caption: Photo by rawpixel on Unsplash
  focal_point: Smart
links:
- icon: twitter
  icon_pack: fab
  name: Follow
  url: 
slides: example
summary: Using machine learning 
tags:
- R
title: Predicting Breast Cancer Using Gene Expression Data in R
url_code: ""
url_pdf: ""
url_slides: ""
url_video: ""
---


# Overview 

In this project, I predicted whether or not a patient has breast cancer only based on their gene expression data. 

From a practical standpoint, this would be a very useful model for both hospitals or a business idea! Cancer is the second leading cause of death in the U.S. [(CDC, 2020)](https://www.cdc.gov/cancer/dcpc/research/update-on-cancer-deaths/index.htm#:~:text=Cancer%20was%20the%20second%20leading,females%20and%20317%2C731%20among%20males.) and as we know, our health care system is overburdened [(Forbes, 2021)](https://www.forbes.com/sites/williamhaseltine/2021/07/26/overwhelmed-us-hospital-systems-a-look-into-the-future/). Thus, having a machine learning model that can detect cancer simply based on gene expression, would expedite the testing process, require less doctor time, and overall improve patient outcomes. 

# About the Model

I used 3 breast cancer datasets to train a model on predicting breast cancer. The models used are: 

- KNN

  - Performs well for continuous, numeric data
  - Simple and easy to understand
  - Has performed well for gene expression classification in the past

- SVM
  
  - Performs well for binary classification
  - Performs well for high number of features
  - Uses distance metric 

- Random Forest

  - Performs well for many data points
  - Bagging
  - Has performed well for gene expression previously

Lastly, I will combine all three into an ensemble model. 

# About the Data

I chose 3 datasets with gene expression data measuring both cancerous and normal breast tissue. The data were from published articles, but were obtained from the Gene Expression Omnibus website hosted by NCBI.

- [Clarke C, Madden SF, Doolan P, Aherne ST et al. Correlating transcriptional networks to breast cancer survival: a large-scale coexpression analysis. Carcinogenesis 2013 Oct;34(10):2300-8.](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568)

- [Planche A, Bacac M, Provero P, Fusco C et al. Identification of prognostic molecular features in the reactive stroma of human breast and prostate cancer. PLoS One 2011;6(5):e18640. ](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26910)

- [Richardson AL, Wang ZC, De Nicolo A, Lu X et al. X chromosomal abnormalities in basal-like human breast cancer. Cancer Cell 2006 Feb;9(2):121-32.](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7904)


```{r, message = FALSE, warning = FALSE}
packages = c("dplyr", "devtools", "affy", "edgeR", "corrplot", "GEOquery", "missForest", 
             "impute", "randomForest", "DMwR", "caret", "gmodels", "pROC")
if_not_install_package <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

install_github("vqv/ggbiplot")
library(ggbiplot)
```


First, I will download the data files from the GEO website using the GEOquery package. The `getGEOSuppFiles()` function downloads the appropriate .CEL files that we will need. 

```{r, message = FALSE, warning = FALSE}
filePath_1 = getGEOSuppFiles("GSE42568")
filePath_2 = getGEOSuppFiles("GSE26910")
filePath_3 = getGEOSuppFiles("GSE7904")

```


Since the files are in a compressed tar format, they must be 'untarred' so they can be read in. We will get the current path so the path can be pasted with the file names when we load in the data. 

```{r, message = FALSE, warning = FALSE}
path <- getwd()
untar(tarfile = row.names(filePath_1), exdir = paste0(path, "/GSE42568_RAW"))
untar(tarfile = row.names(filePath_2), exdir = paste0(path, "/GSE26910_RAW"))
untar(tarfile = row.names(filePath_3), exdir = paste0(path, "/GSE7904_RAW"))

```



The Datasets will be read in using the `ReadAffy()` Function from the [Affy Package](https://www.bioconductor.org/packages/release/bioc/html/affy.html). 

```{r, message = FALSE, warning = FALSE}
cancer_1<- ReadAffy(celfile.path = paste0(path, "/GSE42568_RAW/"), compress = TRUE)
cancer_2<- ReadAffy(celfile.path = paste0(path, "/GSE26910_RAW/"), compress = TRUE)
cancer_3<- ReadAffy(celfile.path = paste0(path, "/GSE7904_RAW/"), compress = TRUE)

```

# Inspect Objects
```{r, message = FALSE, warning = FALSE}
cancer_1
```
