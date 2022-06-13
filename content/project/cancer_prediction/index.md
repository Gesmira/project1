---
title: "Predicting Breast Cancer Using Gene Expression Data in R"
external_link: ''
date: '2022-05-03T00:00:00Z'
links:
- icon: twitter
  icon_pack: fab
  name: Follow
  url: null
slides: example
summary: Using machine learning
tags: R
image:
  caption: Photo by rawpixel on Unsplash
  focal_point: Smart
url_code: ''
url_pdf: ''
url_slides: ''
url_video: ''
---


# Overview 

In this project, I predicted whether or not a patient has breast cancer based on their gene expression data. 

Cancer is the second leading cause of death in the U.S. [(CDC, 2020)](https://www.cdc.gov/cancer/dcpc/research/update-on-cancer-deaths/index.htm#:~:text=Cancer%20was%20the%20second%20leading,females%20and%20317%2C731%20among%20males.) and our health care system is overburdened [(Forbes, 2021)](https://www.forbes.com/sites/williamhaseltine/2021/07/26/overwhelmed-us-hospital-systems-a-look-into-the-future/). Thus, having a machine learning model that can detect cancer simply based on gene expression, would expedite the testing process, require less doctor time, and overall improve patient outcomes. 

# About the Models

I used the following three models:

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

I also combined all three into an ensemble model. 

# About the Data

To train the models, I chose 3 datasets with gene expression data measuring both cancerous and normal breast tissue. The data were from published articles, but were obtained from the Gene Expression Omnibus website hosted by NCBI.

- [Clarke C, Madden SF, Doolan P, Aherne ST et al. Correlating transcriptional networks to breast cancer survival: a large-scale coexpression analysis. Carcinogenesis 2013 Oct;34(10):2300-8.](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568)

- [Planche A, Bacac M, Provero P, Fusco C et al. Identification of prognostic molecular features in the reactive stroma of human breast and prostate cancer. PLoS One 2011;6(5):e18640. ](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26910)

- [Richardson AL, Wang ZC, De Nicolo A, Lu X et al. X chromosomal abnormalities in basal-like human breast cancer. Cancer Cell 2006 Feb;9(2):121-32.](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7904)

# Steps 

### Install and load packages 

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

### Get Data 

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

### Inspect Objects
```{r, message = FALSE, warning = FALSE}
cancer_1
```

We see that this is an AffyBatch object. Thus, to get the count data so we can run exploratory analysis on the object, we will use the `exprs()` function.

```{r, message = FALSE, warning = FALSE}
# Dataset 1 
cancer_1_exprs <- exprs(cancer_1)

# Dataset 2
cancer_2_exprs <- exprs(cancer_2)
# Remove first 12 samples as they are from prostate, not from breast
cancer_2_exprs <- cancer_2_exprs[,-c(1:12)]

# Dataset 3
cancer_3_exprs <- exprs(cancer_3) 
# Remove organelle samples, as we only want breast tissue samples
cancer_3_exprs <- cancer_3_exprs[,-c(44:55)]

```

### Rename Samples

Taking a peak at the column names, we see that they are long with repetitive information. 

```{r, message = FALSE, warning = FALSE}
head(colnames(cancer_1_exprs))
```

Thus, we will rename the column names, representing each sample, so they only contain the GSM ID number, which is in the first 10 digits. 

```{r, message = FALSE, warning = FALSE}
colnames(cancer_1_exprs) <- substr(colnames(cancer_1_exprs),start=1,stop=10)
colnames(cancer_2_exprs) <- substr(colnames(cancer_2_exprs),start=1,stop=9)
colnames(cancer_3_exprs) <- substr(colnames(cancer_3_exprs),start=1,stop=9)
```

### Library Sizes

We next can see how many samples we have in each dataset. 

```{r, message = FALSE, warning = FALSE}
cat('Samples in 1st Dataset:', ncol(cancer_1_exprs), '\n')
cat('Samples in 2nd Dataset:', ncol(cancer_2_exprs), '\n')
cat('Samples in 3rd Dataset:', ncol(cancer_3_exprs))
```

Next, we look at the library sizes. The size of the library represents how many transcripts there were from the experiment for each sample within the dataset. Different library sizes means there is different sequencing depth per sample [(Penn State, 2018)](https://online.stat.psu.edu/stat555/node/13/#:~:text=Library%20Size,total%20number%20of%20mapped%20reads.).

```{r, fig.width=10, fig.height = 4}
barplot(colSums(cancer_1_exprs) /1e6, xlab = '', las = 2, 
        main = "Dataset 1: Library Size Per Sample")

barplot(colSums(cancer_2_exprs) /1e6, xlab = '', las = 2, 
        main = "Dataset 2: Library Size Per Sample")

barplot(colSums(cancer_3_exprs) /1e6, xlab = '', las = 2, 
        main = "Dataset 3: Library Size Per Sample")
```

We see that the library size varies greatly per sample within each dataset. This shows us that we will need to perform normalization. If we didn't, a gene may look like its more highly expressed just because it is part of a sample with a larger library size, and not because it's actually more highly expressed. 


### Sample Distribution

Next, we evaluate the distribution of the samples. We look at a density plot of the probe intensity for each sample.

We use the cpm() function to get the counts per million per sample. This puts our data in a more understandable format. We will use density plots with the log transformed counts of all samples. [(Gene Expression Analysis)](https://mkempenaar.github.io/gene_expression_analysis/chapter-3.html).

```{r, message = FALSE, warning = FALSE}
counts <- cpm(cancer_1_exprs) # get counts per million 
plotDensity(log2(counts), 
            lty=c(1:ncol(counts)), xlab='Log2(count)',
            main='Dataset 1: Expression Distribution')

counts_2 <- cpm(cancer_2_exprs) # get counts per million 
plotDensity(log2(counts_2), 
            lty=c(1:ncol(counts_2)), xlab='Log2(count)',
            main='Dataset 2: Expression Distribution')


counts_3 <- cpm(cancer_3_exprs) # get counts per million 
plotDensity(log2(counts_3), 
            lty=c(1:ncol(counts_3)), xlab='Log2(count)',
            main='Dataset 3: Expression Distribution')

```

We see that the samples have a similar distribution, but some normalization will still be necessary. 

### Data Cleaning and Shaping

Next, we'll use the rma() function to normalize the data, as well as a couple other steps. This function runs the RMA pipeline, also known as the Robust Multi-Array Average, which is the standard for gene expression data [(Introduction to PreProcessing: RMA)](https://www.usu.edu/math/jrstevens/stat5570/1.4.preprocess_4up.pdf). It performs three steps: 

1.  __Background Correction__ 
  - This removes any potential noise in the measurement, where measurements are impacted by their location in the array, which confounds their actual expression. 

2. __Normalization__
  - This normalizes across samples, correcting for the library size differences mentioned previously. The method is entitled quantile normalization, which puts every sample in the same distribution. 
  
3. __Summarization__
  - This combines the probe intensities of every probe set. 
  - Initially we have 1,354,896 columns w all the probe sets: 
  
```{r, message = FALSE, warning = FALSE}
nrow(cancer_1_exprs)
```

  - but multiple probes map to one gene. So we must rename the features to the gene names and combine the counts of the probes that map to the same gene. 
  
```{r, message = FALSE, warning = FALSE}
# Run the RMA pipeline 
normalized.data.1 <- rma(cancer_1)
normalized.data.2 <- rma(cancer_2)
normalized.data.3 <- rma(cancer_3)
# Extract the Counts
exprs_data_1<- exprs(normalized.data.1)
exprs_data_2<- exprs(normalized.data.2)
exprs_data_3<- exprs(normalized.data.3)
# Create a data frame with rows as samples and columns as genes 
exprs_df_1 <- as.data.frame(t(exprs_data_1))
exprs_df_2 <- as.data.frame(t(exprs_data_2))
exprs_df_3 <- as.data.frame(t(exprs_data_3))

# Since we are restarting from the count data, I will repeat the subsetting and column renaming done above 
# Remove first 12 samples as they are from prostate, not from breast
exprs_df_2 <- exprs_df_2[-c(1:12),]

# Remove organelle samples, as we only want breast tissue samples
exprs_df_3 <- exprs_df_3[-c(44:55),]

row.names(exprs_df_1) <- substr(row.names(exprs_df_1),start=1,stop=10)
row.names(exprs_df_2) <- substr(row.names(exprs_df_2),start=1,stop=9)
row.names(exprs_df_3) <- substr(row.names(exprs_df_3),start=1,stop=9)
```

### Identify and Replace Missing Values

As we can confirm, there are no missing values in our data. We check this by applying the `is.na()` function to every column in each dataset. 
```{r, message = FALSE, warning = FALSE}
apply(is.na(exprs_df_1), 2, which)
apply(is.na(exprs_df_2), 2, which)
apply(is.na(exprs_df_3), 2, which)
```

We will randomly create some missing values, and then impute them. 

We will use the [missForest package](https://cran.r-project.org/web/packages/missForest/missForest.pdf) to do so. Since we have over 6.6 Million data points, just in the first dataset, iI will only make 0.1% of each dataset as missing values. 
```{r, message = FALSE, warning = FALSE}
# Produce NAs at a .1% of the dataset 
exprs_df_1_na <- missForest::prodNA(exprs_df_1,noNA=0.001)
exprs_df_2_na <- missForest::prodNA(exprs_df_2,noNA=0.001)
exprs_df_3_na <- missForest::prodNA(exprs_df_3,noNA=0.001)
```

By looking at the sum of all missing values, we see that we have .1% of 6.6 Million data points that are now missing values in our first dataset. 
```{r, message = FALSE, warning = FALSE}
sum(rowSums(is.na(exprs_df_1_na)))
sum(rowSums(is.na(exprs_df_2_na)))
sum(rowSums(is.na(exprs_df_3_na)))
```

Next, we will use the [impute.knn](https://rdrr.io/bioc/impute/man/impute.knn.html) function from the impute package to impute our missing values using 10 Nearest Neighbors. The function finds the 10 nearest neighbors using euclidean distance by looking at columns in which the gene measurement is not missing. Then, the function averages the nearest neighbor values for that gene. 
```{r, message = FALSE, warning = FALSE}
# We take the transpose of the df as the impute.knn function takes row names as genes 
exprs_df_1_na_t <- as.data.frame(t(exprs_df_1_na))
exprs_df_2_na_t <- as.data.frame(t(exprs_df_2_na))
exprs_df_3_na_t <- as.data.frame(t(exprs_df_3_na))
```

Note that I did not include the actual imputation function output as the code lets out a warning message that is impossible to suppress, but the code is included in the .Rmd file. The command is like so: `exprs_df_1_na_t <- impute.knn(as.matrix(exprs_df_1_na_t) , k = 10)` for each dataframe. 

```{r, message = FALSE, warning = FALSE, include = FALSE}
# Run impute.knn w 10 nearest neighbors 
# Data must be in matrix form
exprs_df_1_na_t <- impute.knn(as.matrix(exprs_df_1_na_t) , k = 10)

# Dataset 2
exprs_df_2_na_t <- impute.knn(as.matrix(exprs_df_2_na_t) , k = 10)

# Dataset 3
exprs_df_3_na_t <- impute.knn(as.matrix(exprs_df_3_na_t) , k = 10)
```

```{r, message = FALSE, warning = FALSE}
# Recreate dataframe in original setup 
exprs_1_pp <- as.data.frame(t(exprs_df_1_na_t$data))
exprs_2_pp <- as.data.frame(t(exprs_df_2_na_t$data))
exprs_3_pp <- as.data.frame(t(exprs_df_3_na_t$data))
```

We then validate that the NA values were imputed.
```{r, message = FALSE, warning = FALSE}
sum(rowSums(is.na(exprs_1_pp)))
sum(rowSums(is.na(exprs_2_pp)))
sum(rowSums(is.na(exprs_3_pp)))
```

Next, we must add a column with the annotation of Cancer vs. Normal. This is not included previously as the downloaded data is only the counts data. Thus, the annotation is found from the original website where the data was downloaded from. 

```{r, message = FALSE, warning = FALSE}
# Add column of annotation for Dataset 1 
exprs_1_pp["Type"] <- 0
# First 17 samples are normal 
exprs_1_pp$Type[1:17] <- 'Normal'
# All the rest are cancer 
exprs_1_pp$Type[18:121] <- 'Cancer'
# Move Type to beginning of dataframe
exprs_1_pp <- exprs_1_pp %>% relocate(Type, .before = `1007_s_at`)


# Dataset 2 
exprs_2_pp["Type"] <- 0
# In this sample, every other sample is normal tissue, and every other is cancer 
exprs_2_pp$Type[c(1, 3, 5, 7, 9, 11)] <- 'Normal'
exprs_2_pp$Type[c(2, 4, 6, 8,  10, 12)] <- 'Cancer'
exprs_2_pp <- exprs_2_pp %>% relocate(Type, .before = `1007_s_at`)

# Dataset 3 
exprs_3_pp["Type"] <- 0
# First 42 samples are Cancerous
exprs_3_pp$Type[c(1:42)] <- 'Cancer'
# All the rest are normal
exprs_3_pp$Type[c(43:50)] <- 'Normal'
# Move Type to beginning of dataframe
exprs_3_pp <- exprs_3_pp %>% relocate(Type, .before = `1007_s_at`)
```

### Corellation 

Since there are 50,000 genes, it is computationally expensive to figure out the correlation between all of them and difficult to visualize. Furthermore,many genes are housekeeping genes and are irrelevant in differentiating between normal and cancerous samples, and because there are so many genes, it is likely that many will be highly correlated if they perform similar functions. 

I will just be selecting a subset of the most variable genes to conduct some preliminary correlation analysis. To do this, I will take the logcounts of all the genes, calculate the coefficient of variation per gene, look at some summary statistics of the distribution of variation for the genes, then select a cut off near the max distribution to get a subset of 20-30 genes. We will then calculate the correlation between these genes and visualize the correlation on a plot. 

We begin with the first dataset. 
```{r, message = FALSE, warning = FALSE}
# Get log counts per column 
logcounts <- cpm(exprs_1_pp[-1], log = TRUE)
# Calculate the coefficient of variation for each column
varcoff <- apply(logcounts, 2, function(x) sd(x) / mean(x))
# Look at summary statistics of all variation 
summary(varcoff)
```

We see that our max value for variation is 0.039. Thus, I will be using a cutoff of 0.025. 
```{r, message = FALSE, warning = FALSE}
# Get 10% of the top variable genes and their coefficient of variation 
variable <- varcoff[rank(varcoff) / length(varcoff) > 1 - .1]

# Get most variable genes (higher than .03)
variable_only <- variable[variable > .029]

# Get names of most variable genes 
variable_genes <- names(variable_only)

# Subset the initial dataset by the most variable genes 
exprs_1_pp_var <- exprs_1_pp[, variable_genes]

# Calculate correlations between these genes
cor_1 <- cor(exprs_1_pp_var)

# Visualize correlation
corrplot(cor_1, method="circle", tl.pos = "n")

```

Out of the variable genes, we see that many are highly positively correlated, as represented by the darker blue circles, and many are negatively correlated, as represented by the red circles. I have omitted labels as their overlap makes them unlegible and their names are unimportant. 

Although multicolinearity can present an issue for machine learning models, we will be conducting PCA to select our features. All PCs are uncorrelated, thus this will solve that problem. More details on this below. 

We repeat the same process for the other datsaets.

```{r, message = FALSE, warning = FALSE}
logcounts <- cpm(exprs_2_pp[-1], log = TRUE)
varcoff <- apply(logcounts, 2, function(x) sd(x) / mean(x))

variable <- varcoff[rank(varcoff) / length(varcoff) > 1 - .1]

variable_only <- variable[variable > .025]

variable_genes <- names(variable_only)

exprs_2_pp_var <- exprs_2_pp[, variable_genes]

cor_2 <- cor(exprs_2_pp_var)

corrplot(cor_2, method="circle", tl.pos = "n")
```

We see more genes are highly correlated in this dataset, but again that will be resolved with PCA. 

```{r, message = FALSE, warning = FALSE}
logcounts <- cpm(exprs_3_pp[-1], log = TRUE)
varcoff <- apply(logcounts, 2, function(x) sd(x) / mean(x))
variable <- varcoff[rank(varcoff) / length(varcoff) > 1 - .1]
variable_only <- variable[variable > .023]

variable_genes <- names(variable_only)

exprs_3_pp_var <- exprs_3_pp[, variable_genes]

cor_3 <- cor(exprs_3_pp_var)
corrplot(cor_3, method="circle", tl.pos = "n")
```

### Merge Data

Now, I'll merge all the data. First we check that all of the column names (or genes) are the same in all three datasets. 

```{r, message = FALSE, warning = FALSE}
setdiff(colnames(exprs_1_pp), colnames(exprs_2_pp))
setdiff(colnames(exprs_3_pp), colnames(exprs_2_pp))
```
Since they are, we can use rbind() to merge the dataframes

```{r, message = FALSE, warning = FALSE}
exprs_combined <- rbind(exprs_1_pp, exprs_2_pp, exprs_3_pp)
cat("Total Samples:", nrow(exprs_combined), "\n")
cat("Total Genes:", ncol(exprs_combined), "\n")
table(exprs_combined$Type)
```

There are significantly more cancerous samples than normal. Although this is expected as the experiments all measured mostly cancerous tissue, we will adjust for it below after splitting our datasets into training and testing. 


# Training and Testing Split 
We will split our data so that 70% is used for training and 30% is used for testing. We will create the test sets using random sampling as the data is ordered. 

```{r, message = FALSE, warning = FALSE}
# Setting a fixed seed
set.seed(187)
# 70% for training, 30% for testing
test_size = floor(0.7*nrow(exprs_combined)) 
inds <- sample(seq_len(nrow(exprs_combined)), size = test_size, replace=F) 
exprs.training <- exprs_combined[inds, ]
exprs.testing <- exprs_combined[-inds,]
```

We'll now take a peak at the composition of our training and testing sets. 

```{r, message = FALSE, warning = FALSE}
table(exprs.training$Type)
```

Although the sets are split 70/30, we see again that our 'Normal' variable is underrepresented. Thus, we will use the SMOTE algorithm from the `DMwR` package to perform over-sampling of our minority class. We set perc.over as 200, which means for every case in the original set 2 new cases will be created for the minority class (22*3 = 66). The function also underrepresents the majority variable.  

```{r, message = FALSE, warning = FALSE}
# Convert type to factor for use with SMOTE
exprs.training$Type <- as.factor(exprs.training$Type)
train.smote <- SMOTE(Type ~., exprs.training, perc.over = 200)
```

We look at the composition of our training set now and see it is much more balanced. 
```{r, message = FALSE, warning = FALSE}
table(train.smote$Type)
```


### PCA 

Next, we run principal component analysis only on the training set. PCA is a dimensional reduction technique that gives a set of "principal components" which aim to capture the variance in the data as a set of _uncorrelated_ variables which are linear combinations of the original variables. Thus, it's a great technique for this type of analysis where we have more than 50,000 genes, many of which are correlated. 


```{r, message = FALSE, warning = FALSE}
# Remove type column from dataset
exprs.pca.training.smote <- prcomp(train.smote[,2:54676], center = TRUE,scale. = TRUE)

```

For visualization purposes, I now plot the 1st to 2nd PCs, the 2nd to 3rd, and the 3rd to 4th. We see a bit of separation between the Cancer vs. Normal datasets start to form, especially in the 3rd to 4th plot. Using multiple principal components should capture the variance between the datasets well.

```{r, message = FALSE, warning = FALSE}

ggbiplot(exprs.pca.training.smote, choices= c(1,2),obs.scale = 1, var.scale = 1, 
         var.axes = FALSE,
         labels=train.smote$Type, ellipse = TRUE, groups = train.smote$Type)
```

```{r, message = FALSE, warning = FALSE}
ggbiplot(exprs.pca.training.smote, choices= c(2,3),obs.scale = 1, var.scale = 1, 
         var.axes = FALSE,
         labels=train.smote$Type, ellipse = TRUE, groups = train.smote$Type)
```

```{r, message = FALSE, warning = FALSE}
ggbiplot(exprs.pca.training.smote, choices= c(3,4),obs.scale = 1, var.scale = 1, 
         var.axes = FALSE,
         labels=train.smote$Type, ellipse = TRUE, groups = train.smote$Type)
```

Next, we plot the cumulative sum of the captured variance by the prinicipal components. We also plot an 'elbow plot' of the variance per principal component too. These plots help us understand how many principal components we need. Since it looks like 99% of the variance is covered by about 75 principal components, that is how many we will move forward with. 

```{r, message = FALSE, warning = FALSE}
plot(cumsum(exprs.pca.training.smote$sdev^2/sum(exprs.pca.training.smote$sdev^2)))
abline(h =0.99, col = "darkgreen")
```
```{r, message = FALSE, warning = FALSE}
screeplot(exprs.pca.training.smote, type = "l", npcs = 15, 
          main = "Screeplot of the first 10 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
```


We then take the 75 PCs and create a dataframe with them. Thus, this is a huge reduction from our original number of features (~50,000).

```{r, message = FALSE, warning = FALSE}
exprs.pca.training.df <- as.data.frame(exprs.pca.training.smote$x[,1:70])

```


### Get target variable labels for the training and testing data sets

Next, we do a bit rearranging of our target variables and dataframes for training and testing. 

```{r, message = FALSE, warning = FALSE}
# Change Target to Binary 0,1
train.smote$Type <- ifelse(train.smote$Type == "Cancer", 1, 0)

# Make vector of just training labels
exprs.training.labels <- train.smote$Type

# Make vector of just testing labels
exprs.testing.labels <- exprs.testing$Type

# Change Target to Binary 0,1
exprs.testing.labels <- ifelse(exprs.testing.labels == "Cancer", 1, 0)

# Remove "Type" column from testing dataframe
exprs.testing <- exprs.testing[,-1]

# Make training dataframe target a factor 
exprs.pca.training.df$Type <- as.factor(train.smote$Type)
```

Next, we project the PCA structure from the training set onto the testing set. Note that we did not run PCA on the entire dataframe before splitting, as this would have led to an information leak and would not have let the models be blind to the testing data. However, we will project the training PCA structure onto the testing set so that the features are in the same space. 

```{r, message = FALSE, warning = FALSE}
exprs.pca.testing <- predict(exprs.pca.training.smote, newdata = exprs.testing)
exprs.pca.testing.df <- as.data.frame(exprs.pca.testing)
```

# Model Training

## Cross Validation 

For all of our models, we will train them using the caret package. We will use cross validation, repeated 5 times to fine tune the parameters for our models. 
```{r, message = FALSE, warning = FALSE}
ctrl_caret <- trainControl(method = "cv", number = 5)
```

# KNN
The first model we will train is K Nearest Neighbors. This was chosen as it performs very well for continous, numeric data. The curse of dimensionality is an issue w KNN so it is good that we have already reduced our feature space. The KNN model is also simple and easy to understand, which is a good starting point as the later models will be more complicated. Lastly, it tends to perform well for gene expression classification in the past. 

As mentioned, we will train the algroithm using our cross validation control, then visualize the output. 

```{r, message = FALSE, warning = FALSE}
knnFit <- train(Type ~ ., data = exprs.pca.training.df, method = "knn", 
                trControl = ctrl_caret)
knnFit
```

We see that k=5 is our highest accuracy, so it is now used in the model. 

Next, we get the predictions using our kNN model and the testing data. 
```{r, message = FALSE, warning = FALSE}
cancer_knn_pred <- predict(knnFit,newdata = exprs.pca.testing.df )
```

# SVM

Our next model is a support vector machine. I chose this model as it performs very well for binary classifications which is what we are attempting to do. It also performs well with a high number of features, however, this is of less concern as we did do the dimensionality reduction. It also uses a distance metric as part of its algorithm, it will also be good for this numeric type of data. 

```{r, message = FALSE, warning = FALSE}
svmFit <- train(Type ~ ., data = exprs.pca.training.df, method = "svmLinear", 
                trControl = ctrl_caret)
svmFit
```


```{r, message = FALSE, warning = FALSE}
cancer_svm_pred <- predict(svmFit, exprs.pca.testing.df)
```


# Random Forest

Lastly, we will use a random forest model. I chose this model as it performs well when there are many data points and it has been previously used for gene expression analysis. Lastly, it has been known to be one of the most accurate machine learning models, so I wanted to see how it performs here.

Using the random forest model will also allow us to test out a bagging method. 
```{r, message = FALSE, warning = FALSE}
rfFit <- train(Type ~ ., data = exprs.pca.training.df, method = "rf", 
               trControl = ctrl_caret)
rfFit

```


```{r, message = FALSE, warning = FALSE}
cancer_rf_pred <- predict(rfFit, exprs.pca.testing.df)
```


# Evaluation of Models 

Simply, we first visualize our results using a CrossTable to determine the false positives and negatives. Accuracy alone is not a great metric as we have unbiased classes, so the minority class has little impact on the final accuracy, even if all of the minority class are being misclassified. 

## KNN 
```{r, message=FALSE}
# Compare the model predicted labels vs. the actual labels 
comparison_knn <- CrossTable(x = exprs.testing.labels, y = cancer_knn_pred, 
                             prop.chisq = FALSE)
```

## SVM 
```{r, message = FALSE, warning = FALSE}
comparison_svm <- CrossTable(x = exprs.testing.labels, y = cancer_svm_pred, 
                             prop.chisq = FALSE)
```

## RF
```{r, message = FALSE, warning = FALSE}
comparison_rf <- CrossTable(x = exprs.testing.labels, y = cancer_rf_pred, 
                            prop.chisq = FALSE)

```

Thus, we will begin by looking at precision, recall, and the F score for each model. Precision shows the number of true positive results over all positives (True or False), while recall is the number of true positive results over true positives and false negative results. The F1 score uses both of these metrics thus we should aim to have a high F1 score.

These can all be found using the _`confusionMatrix()`_ function in the _`caret`_ package. We also create an ROC curve to determine the area under the curve (AUC). ROC curves are great for evaluating binary classifiers. The axises represent the true positive rate and the true negative rate. Thus, a high area under the curve is indicative of a good model. An AUC of 0.5 is as good as random guessing. 

Below, we calculate these stats for each model, and put them all into a dataframe for easier visualization. 

```{r, message = FALSE, warning = FALSE}
stats_compare <- data.frame(row.names = c("Accuracy", "Precision", "Recall", 
                                          "F1 Score", "AUC"))


cm_knn <- confusionMatrix(as.factor(exprs.testing.labels), as.factor(cancer_knn_pred))
cm_svm <- confusionMatrix(as.factor(exprs.testing.labels), as.factor(cancer_svm_pred))
cm_rf <- confusionMatrix(as.factor(exprs.testing.labels), as.factor(cancer_rf_pred))

stats_compare$KNN <- c(cm_knn$overall[["Accuracy"]], cm_knn$byClass[["Pos Pred Value"]], 
                       cm_knn$byClass[["Sensitivity"]], cm_knn$byClass[["F1"]], 
                       roc(exprs.testing.labels, as.numeric(cancer_knn_pred))$auc)

stats_compare$SVM <- c(cm_svm$overall[["Accuracy"]], cm_svm$byClass[["Pos Pred Value"]], 
                       cm_svm$byClass[["Sensitivity"]], cm_svm$byClass[["F1"]], 
                       roc(exprs.testing.labels, as.numeric(cancer_svm_pred))$auc)
stats_compare$RF <- c(cm_rf$overall[["Accuracy"]], cm_rf$byClass[["Pos Pred Value"]],
                      cm_rf$byClass[["Sensitivity"]], cm_rf$byClass[["F1"]], 
                      roc(exprs.testing.labels, as.numeric(cancer_rf_pred))$auc)
stats_compare
```

For KNN for example, we see we have a higher precision, so normal samples are not often being mislabeled as cancerous, however, we have a low recall, which is indicated by the fact that some of our cancerous data sets are being classified as normal. This does not have too much of an impact on our our AUC as it is the majority class that is being misrepresented slightly. We see a similar issue, although slightly better, with RF. SVM performs the best out of the 3, as all of its metrics are high! 

I believe this may be the case as there might be a small amount of differentiation between some cancerous samples and normal samples, if it is an early cancerous case. Thus, the algorithm wouldn't be able to tell the difference between early cancerous samples and normal samples if there is a much larger difference between normal samples and much more cancerous samples. 

I did attempt to change the composition of the training data set using SMOTE to see if this would make a difference. Although there were some improvements on the accuracy on the training sets, the accuracy on the testing sets decreased, thus the model was overfitting to the training set. 

# Ensemble Model

We will now build an ensemble model using the above three models: KNN, SVM, and Random Forest. We start by building a function to get the most common class from the 3 predictions. We then build a function called `predictCancer` which will act as a weighted average ensemble model. It takes in a test case, runs the three models on it, then returns a weighted avergae of the predictions. 

```{r, message = FALSE, warning = FALSE}
test_df <- exprs.pca.testing.df
predictCancer <- function(test_df) {
  # Get test case into appropriate format for models
  test_df <- data.frame(t(test_df))
  # KNN
  test_df$knnpred <- predict(knnFit, as.data.frame(test_df))
  # SVM
  test_df$svmpred <-  predict(svmFit, as.data.frame(test_df))
  # RF
  test_df$rfpred <- predict(rfFit, as.data.frame(test_df))
  # Weighted Average
      # 2 times weight to knn bc it is effective at minimizing false positives
      # 3 times weight to svm because it was the most accurate
      # 1 times weight to rf because it was also decently accurate
      # Divide by 6 to get average 
  test_df$weighted_avg <- ((test_df$knnpred == 1)*2 +  (test_df$svmpred == 1)*3 
                           + (test_df$rfpred == 1))/6
  # Turn the number into a binary class prediction of 0 or 1 
  test_df$binary_pred <- factor(ifelse(test_df$weighted_avg > .5, 1, 0), levels = c(0,1))
  
}

# Apply the function to our test dataframe
test_df$prediction <- apply(test_df, 1, predictCancer)
test_df$prediction
```

```{r, message = FALSE, warning = FALSE}
comparison_ensemble <- CrossTable(x = exprs.testing.labels, y = as.factor(test_df$prediction), 
                                  prop.chisq = FALSE)
```

```{r, message = FALSE, warning = FALSE}
cm_ens <- confusionMatrix(as.factor(exprs.testing.labels), as.factor(test_df$prediction))

stats_compare$Ensemble <- c(cm_ens$overall[["Accuracy"]], 
                            cm_ens$byClass[["Pos Pred Value"]], 
                            cm_ens$byClass[["Sensitivity"]], cm_ens$byClass[["F1"]], 
                            roc(exprs.testing.labels, as.numeric(test_df$prediction))$auc)
stats_compare
```

We see that our ensemble performs just as well as the SVM model. In the future, testing this model on even more data if it becomes available would prove if the ensemble is actually any better, or if SVM just performs the best. If SVM just performs the best in all cases, then it would be better to only use this model as the other ones would lower its accuracy and it would be less computationally intensive to use only one model. 

# Conclusions 

Overall, we have 100% recall, but are lacking on precision. This means that we are mislabeling some normal samples as cancerous, however this is much preferred to labeling some cancerous samples as normal. If a cancerous sample is labeled as normal, that person will be unlikely to get follow up tests, and the cancer could become worse or fatal. If a normal sample is labeled as cancerous, that person will go and seek follow up testing, at which point it will be seen as normal. Thus, although there is a bit extra inconvenience for that person, it will not happen often, and it is much preferred to the alternative. 

Initially, I was concerned the batch effects of combining data from 3 different experiments would be large, but with normalization, dimensional reduction, and the models chosen, they were still able to perform well. Thus this is a robust model that could be applied to many experiments. 

