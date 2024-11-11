# CCA_musical_physiological_change_points

This repository contains the method for calculating Canonical Correlation Analysis (CCA) using R script. 
The CCA method implemented in this script was used to analyse the music-physiology association patterns in the dataset from the HeartFM study that is a part of the COSMOS project (https://cosmos.isd.kcl.ac.uk/?page_id=2).

CCA is a robust statistical method for determining the relationships between two data sets. It is a multidimensional analysis that determines the relationships between two sets of data, providing more complete and nuanced outputs than linear regression or multiple comparisons. 

INPUT DATA REQUIREMENTS

The input data are two matrices, X and Y, with the sizes m1 x n1 and m2 x n2, whereas m1 is equal to m2 (m1==m2), and it is the number of samples in the database. The number of columns (number of features in each matrix) can vary between the matrices.  

OUTPUT 

This script calculates the following output results:

- coefficients of the canonical functions (variates)
- canonical loadings
- Wilk's lambda values
- relative covariance value for each canonical function
- correlation coefficient for each canonical function
- the result of the CCA on permuted physiological data (matrix Y)


The script also visualises the results by plotting:
- Covariance and canonical correlation values at each canonical function,
- Scatter plots of physiological and musical variables were obtained for the first (U1, V1) and the second (U2, V2) canonical varieties.
- Mean values and 95%CI of physiological variable values (V1, V2) of data points from the correlation plot (scatter plot) for each piece
separately, showing the different contributions of change points to these results.
- Distributions of the correlation values of selected canonical variates obtained from the 1000 permutations. 
