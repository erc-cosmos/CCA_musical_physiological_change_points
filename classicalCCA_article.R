# Load necessary libraries
library(PMA)
library(rio)
library(CCP)
library(ggplot2)
library(boot)
library(CCA)
library(gridExtra)
library(smplot2)
library(mgcv)
library(rpart)
library(rpart.plot)
library(reshape2)
library(dunn.test)
library("Hmisc")
library(corrplot)
library("PerformanceAnalytics")
library(extrafont)
loadfonts(device = "all")
library(showtext)
library(vegan)
font_add("Abadi MT", "/Users/mateuszsolinski/Desktop/repositories/Cosmos/data_analysis/hypertensive_study/abadi-mt.ttf")
showtext_auto()

#calculate 95% CI
calculate_ci <- function(x) {
  n <- length(x)
  mean_x <- mean(x)
  sd_x <- sd(x)
  error <- qt(0.975, df = n - 1) * (sd_x / sqrt(n))
  lower_bound <- mean_x - error
  upper_bound <- mean_x + error
  return(c(mean = mean_x, lower = lower_bound, upper = upper_bound))
}

# Function to compute covariance between two matrices column-wise
compute_covariance <- function(scores_X, scores_Y) {
  n <- ncol(scores_X)
  covariance <- numeric(n)
  correlation <- numeric(n)
  
  for (i in 1:n) {
    covariance[i] <- cov(scores_X[, i], scores_Y[, i])
  }
  
  return(covariance)
}
# Function to compute covariance between two matrices column-wise
compute_correlation <- function(scores_X, scores_Y) {
  n <- ncol(scores_X)
  
  correlation <- numeric(n)
  
  for (i in 1:n) {
    correlation[i] <- cor(scores_X[, i], scores_Y[, i])
  }
  
  return(correlation)
}

# load matrices X (musical data) and Y (physiological data) containing values collected from Gaussian Kernel Functions created at the change points

# you can use this example
X <- import('./musical_featuresX_32_all.csv')
Y <- import('./physio_featuresY_32_all.csv')

#additional matrix for the analysis of the CCA for separate tracks 
clinical_date <-import('/Users/mateuszsolinski/Desktop/Cosmos_data_results/Hypertensive_Study/Results/hrv_analysis/Fast_response_RR_annotations/clinical_features_32_all.csv')



#calculate CCA 
ccClassical <- cc(X,Y)
rho <- ccClassical$cor #canonical correlations coefficents
 
# # Perform Sparse CCA
# perm.out <- CCA.permute(X,Y,typex="standard",typez="standard",nperms=100)
# print(perm.out)
#  
# result <- CCA(X,Y,typex="standard",typez="standard", K = dim(Y)[2], penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz, v=perm.out$v.init) 
# print(result, verbose = TRUE)

## Define number of observations, ## and number of dependent and independent variables:
N = dim(X)[1]
p = dim(X)[2]
q = dim(Y)[2]

# Wilks' Lambda = product(1 - canonical_correlations^2)
wilks_lambda <- prod(1 - rho^2)
print(paste("Wilks' Lambda:", wilks_lambda))
print(paste("R2 for the entire model:", 1-wilks_lambda))

## Calculate p-values using F-approximations of some test statistics:
p.asym(rho, N, p, q, tstat = "Wilks")
p.asym(rho, N, p, q, tstat = "Hotelling")
p.asym(rho, N, p, q, tstat = "Pillai")
p.asym(rho, N, p, q, tstat = "Roy")
# ## Plot the F-approximation for Wilks' Lambda, 
# ## considering 3, 2, or 1 canonical correlation(s): 
res1 <- p.asym(rho, N, p, q)

#permuation tests using Wilks lambda
perm.out_test = p.perm(X, Y, nboot = 999, rhostart = 1, type = "Wilks")
perm.out_test
 

#  
# Ensure the canonical vectors are numeric matrices
canonical_vectors_X_orig <- as.matrix(ccClassical$xcoef)
canonical_vectors_Y_orig <- as.matrix(ccClassical$ycoef)

# Project original data onto the canonical vectors 
# (convolution between input matrices and coefficients of the canonical functions)
X_canonical <- as.matrix(X) %*% canonical_vectors_X_orig
Y_canonical <- as.matrix(Y) %*% canonical_vectors_Y_orig

#calculate mean and standard deviation values for the canonical vectors
meanX=mean(X_canonical[,1])
meanY = mean(Y_canonical[,1])
sdX = sd(X_canonical[,1])
sdY = sd(Y_canonical[,1])


# Compute the covariance between the canonical variate scores
covariance_raw <- compute_covariance(X_canonical, Y_canonical)
# compute relative covariance of each canonical function (sum to 100%)
total_covariance <- sum(abs(covariance_raw))
proportion_covariance_explained_orig <- abs(covariance_raw) / total_covariance *100

# Compute the correlation between the canonical variate scores
correlation_orig <- compute_correlation(X_canonical, Y_canonical)


# visualise relative covariance and correlation
df <- data.frame(Component=1:length(proportion_covariance_explained_orig), VarianceExplained=proportion_covariance_explained_orig, CorrelationExplained = correlation_orig)

#covariance
ggplot(df, aes(x=Component, y=VarianceExplained)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue") +
  labs(title="Scree Plot of Covariance explained", x="Canonical Component", y="Variance Explained [%]") +
  theme_minimal()
mean(proportion_covariance_explained_orig)

#correlation
ggplot(df, aes(x=Component, y=CorrelationExplained)) +
  geom_point(color="green", size=3) +
  geom_line(color="green") +
  labs(title="Scree Plot of Canonical Correlation", x="Canonical Component", y="Canonical correlation") +
  theme_minimal()
mean(correlation_orig)

corr_cov_results <- rbind(proportion_covariance_explained_orig,correlation_orig)

# Calculate loadings
loadings_X <- cor(X, X_canonical)
loadings_Y <- cor(Y, Y_canonical)
loadings_X_Y = rbind(loadings_X,loadings_Y)

#export loadings and results (covariance and correlation) (as .csv files)
export(loadings_X_Y,"/Users/mateuszsolinski/Desktop/Cosmos_data_results/Hypertensive_Study/Results/hrv_analysis/Fast_response_RR_annotations/loadings_X_Y_classicalCCA.xlsx")
export(corr_cov_results,"/Users/mateuszsolinski/Desktop/Cosmos_data_results/Hypertensive_Study/Results/hrv_analysis/Fast_response_RR_annotations/corr_cov_results_classicalCCA.xlsx")

 
# Create a data frame
data <- data.frame('Y_canonical_1'=Y_canonical[,1], 
                   'X_canonical_1'=X_canonical[,1], 
                   'Y_canonical_2'=Y_canonical[,2], 
                   'X_canonical_2'=X_canonical[,2],
                   'track'=as.factor(clinical_date$track), 
                   'version'=as.factor(clinical_date$version), 
                   'trackWithoutVersion'=as.factor(clinical_date$trackWithoutVersion)
                   )

# Fit a GAM model
gam_model <- gam(Y_canonical_1 ~  track, data = data)
summary(gam_model)
gam_model <- gam(Y_canonical_1 ~  version + trackWithoutVersion, data = data)

# Display the model summary
summary(gam_model)


mean_canY_track_1 <- aggregate(Y_canonical_1 ~ track, data, calculate_ci)
mean_canY_track_2 <- aggregate(Y_canonical_2 ~ track, data, calculate_ci)
mean_canY_track_1$customisedColors = rep('white',nrow(mean_canY_track_1))
mean_canY_track_1$customisedColors[mean_canY_track_1$track=='Be_1']='magenta1'
mean_canY_track_1$customisedColors[mean_canY_track_1$track=='BG_2']='limegreen'
mean_canY_track_1$customisedColors[mean_canY_track_1$track=='BG_3']='green'
mean_canY_track_1$customisedColors[mean_canY_track_1$track=='BG_4']='darkgreen'
mean_canY_track_1$customisedColors[mean_canY_track_1$track=='Ch1_1']='khaki1'
mean_canY_track_1$customisedColors[mean_canY_track_1$track=='Ch1_2']='khaki3'
mean_canY_track_1$customisedColors[mean_canY_track_1$track=='Ch2_1']='cyan'
mean_canY_track_1$customisedColors[mean_canY_track_1$track=='De_3']='pink'
mean_canY_track_1$customisedColors[mean_canY_track_1$track=='Pr_3']='royalblue'
mean_canY_track_1$customisedColors[mean_canY_track_1$track=='Pr_2']='skyblue'
mean_canY_track_1$customisedColors[mean_canY_track_1$track=='Pr_1']='lightblue'
mean_canY_track_1$customisedColors[mean_canY_track_1$track=='SL_2']='navajowhite3'

mean_canY_track_2$customisedColors = rep('white',nrow(mean_canY_track_2))
mean_canY_track_2$customisedColors[mean_canY_track_2$track=='Be_1']='magenta1'
mean_canY_track_2$customisedColors[mean_canY_track_2$track=='Pr_4']='blue'
mean_canY_track_2$customisedColors[mean_canY_track_2$track=='Pr_3']='royalblue'
mean_canY_track_2$customisedColors[mean_canY_track_2$track=='BG_4']='tomato'
mean_canY_track_2$customisedColors[mean_canY_track_2$track=='BG_3']='gold'
mean_canY_track_2$customisedColors[mean_canY_track_2$track=='Ch1_4']='khaki4'
mean_canY_track_2$customisedColors[mean_canY_track_2$track=='SL_1']='moccasin'
mean_canY_track_2$customisedColors[mean_canY_track_2$track=='De_3']='pink'

# 
result_aov = aov(formula = as.formula('Y_canonical_1 ~ track'),data=data )
summary(result_aov)
tukey_result <- TukeyHSD(result_aov)
tukey_result


# Combine the levels of both categorical variables into a single factor
# combined_factor <- interaction(data$version, data$trackWithoutVersion)
# data['combined_factor'] = combined_factor
# Perform Kruskal-Wallis test on the combined factor
kruskal_test_combined <- kruskal.test(Y_canonical_1 ~ track, data = data)
# Print the result
print(kruskal_test_combined)
dunn_result <- dunn.test(data$Y_canonical_1, g = data$track,  kw = TRUE)
p_values <- dunn_result$P.adjusted

levels_cat1 <- levels(data$track)
p_value_matrix <- matrix(NA, nrow = length(levels_cat1), ncol = length(levels_cat1))
rownames(p_value_matrix) <- levels_cat1
colnames(p_value_matrix) <- levels_cat1

# Fill the matrix with p-values
for (i in 1:(length(levels_cat1) - 1)) {
  for (j in (i + 1):length(levels_cat1)) {
    p_value_matrix[i, j] <- p_values[which(dunn_result$comparisons == paste(levels_cat1[i], levels_cat1[j], sep = " - "))]
    p_value_matrix[j, i] <- p_value_matrix[i, j]
  }
}

# Convert the matrix to a data frame suitable for ggplot
p_value_melted <- melt(p_value_matrix, na.rm = TRUE)
p_value_melted$color <- ifelse(p_value_melted$value < 0.05, 1, 0)

# Plot the heatmap
sepP1 <- ggplot(p_value_melted, aes(x = Var1, y = Var2, fill = color)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red", na.value = "grey") + 
  labs(title = "Heatmap of Dunn's Test P-Values", x = "", y = "") + 
  theme_minimal() +   theme(axis.text.x = element_text(angle = 90, hjust = 1))  + theme(legend.position = "none") # Rotate x-axis labels

kruskal_test_combined <- kruskal.test(Y_canonical_2 ~ track, data = data)
# Print the result
print(kruskal_test_combined)
dunn_result <- dunn.test(data$Y_canonical_2, g = data$track,  kw = TRUE)
p_values <- dunn_result$P.adjusted

levels_cat1 <- levels(data$track)
p_value_matrix <- matrix(NA, nrow = length(levels_cat1), ncol = length(levels_cat1))
rownames(p_value_matrix) <- levels_cat1
colnames(p_value_matrix) <- levels_cat1

# Fill the matrix with p-values
for (i in 1:(length(levels_cat1) - 1)) {
  for (j in (i + 1):length(levels_cat1)) {
    p_value_matrix[i, j] <- p_values[which(dunn_result$comparisons == paste(levels_cat1[i], levels_cat1[j], sep = " - "))]
    p_value_matrix[j, i] <- p_value_matrix[i, j]
  }
}

# Convert the matrix to a data frame suitable for ggplot
p_value_melted <- melt(p_value_matrix, na.rm = TRUE)
p_value_melted$color <- ifelse(p_value_melted$value < 0.05, 1, 0)

# Plot the heatmap
sepP2 <- ggplot(p_value_melted, aes(x = Var1, y = Var2, fill = color)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red", na.value = "grey") + 
  labs(title = "Heatmap of Dunn's Test P-Values", x = "", y = "") + 
  theme_minimal() +   theme(axis.text.x = element_text(angle = 90, hjust = 1))  + theme(legend.position = "none") # Rotate x-axis labels

plot_corr_sep_pieces_pvalues <- grid.arrange(sepP1, sepP2, ncol = 2)

ggsave("/Users/mateuszsolinski/Desktop/Cosmos_data_results/Hypertensive_Study/Results/hrv_analysis/Fast_response_RR_annotations/plot_corr_sep_pieces_pvalues.pdf", plot = plot_corr_sep_pieces_pvalues, width = 10, height = 5, units = "in")


#prepare list of colors for each track to visualise the CCA outcomes
data$selectedTrack = rep('Other',nrow(data))
data$selectedTrack[data$track=='Pr_4']='Pr_4'
data$selectedTrack[data$track=='Pr_3']='Pr_3'
data$selectedTrack[data$track=='BG_4']='BG_4'
data$selectedTrack[data$track=='BG_3']='BG_3'

data$customisedColors_2 = rep('white',nrow(data))
data$customisedColors_1 = rep('white',nrow(data))

data$customisedColors_1[data$track=='Be_1']='magenta1'
data$customisedColors_1[data$track=='BG_2']='limegreen'
data$customisedColors_1[data$track=='BG_3']='green'
data$customisedColors_1[data$track=='BG_4']='darkgreen'
data$customisedColors_1[data$track=='Ch1_1']='khaki1'
data$customisedColors_1[data$track=='Ch1_2']='khaki3'
data$customisedColors_1[data$track=='Ch2_1']='cyan'
data$customisedColors_1[data$track=='De_3']='pink'
data$customisedColors_1[data$track=='Pr_3']='royalblue'
data$customisedColors_1[data$track=='Pr_2']='skyblue'
data$customisedColors_1[data$track=='Pr_1']='lightblue'
data$customisedColors_1[data$track=='SL_2']='navajowhite3'

data$customisedColors_2[data$track=='Be_1']='magenta1'
data$customisedColors_2[data$track=='Pr_4']='blue'
data$customisedColors_2[data$track=='Pr_3']='royalblue'
data$customisedColors_2[data$track=='BG_4']='tomato'
data$customisedColors_2[data$track=='BG_3']='gold'
data$customisedColors_2[data$track=='Ch1_4']='khaki4'
data$customisedColors_2[data$track=='SL_1']='moccasin'
data$customisedColors_2[data$track=='De_3']='pink'

#create the figure presented:
# first column: the scatter plots with the canonical functions values for musical variables (x-axis) and physiological variables (y-axis) for the first two canonical functions
# second column: the mean values of the physiological variables V (relatedo to physiological data) extracted from the first two canonical functions
 
pX1 <- ggplot(data, aes(x=X_canonical_1, y=Y_canonical_1)) +
  geom_point(shape = 21,          # Use shape 21 to allow fill and color
             fill = data$customisedColors_1,       # Set the fill color 
             color = "black",     # Set the outline color (black)
             size = 1.5,            # Set the size of points
             stroke = 0.5,        # Set the thickness of the outline
             alpha = 0.7) +
  xlim(-5, 5) +           # Set x-axis limits 
  ylim(-5, 5) +         # Set y-axis limits 
  # geom_line(color="blue") +
  stat_smooth(method=lm, size=0.5, color='black') +
  # sm_statCorr( color = "black") +
  labs(title="First canonical variate", x=bquote('Musical variable'-U[1]), y= bquote('Physiologial variable'-V[1])) +
  theme_classic() + theme(text = element_text(size = 12, family = "Abadi MT"), plot.title = element_text(size = 12))

# Create a boxplot with points in the background
ggplot(data, aes(x = as.factor(track), y = Y_canonical_1)) +
  geom_jitter(aes(color = as.factor(customisedColors_1)), width = 0.2, alpha = 0.6) +  # Add jittered points with transparency
  geom_boxplot(aes(fill = as.factor(customisedColors_1)), alpha = 0.6) +  # Add boxplot with transparency
  theme_minimal() + 
  labs(x = "Number of Cylinders", y = "Miles Per Gallon (mpg)", title = "Boxplot with Points in Background")


pX2 <- ggplot(mean_canY_track_1, aes(x = track, y =  Y_canonical_1[,1], fill = customisedColors)) +
  geom_point(shape = 21,size = 3,fill = mean_canY_track_1$customisedColors) + # Mean points
  geom_errorbar(aes(ymin = Y_canonical_1[,2], ymax = Y_canonical_1[,3]), width = 0.2) + # Confidence intervals
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.25) + 
  scale_y_continuous(breaks = seq(-2, 2, by = 0.5)) +
  labs(x = "", y = bquote('Physiologial variable'-V[1]), title = "Separate pieces") +
  theme_classic(base_size=12)  +   theme(text = element_text(size = 12, family = "Abadi MT"), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 0, size=10), plot.title = element_text(size = 12)) 

pY1 <- ggplot(data, aes(x=X_canonical_2, y=Y_canonical_2)) +
  geom_point(shape = 21,          # Use shape 21 to allow fill and color
             fill = data$customisedColors_2,       # Set the fill color 
             color = "black",     # Set the outline color (black)
             size = 1.5,            # Set the size of points
             stroke = 0.5,        # Set the thickness of the outline
             alpha = 0.5) +
  xlim(-5, 5) +           # Set x-axis limits 
  ylim(-5, 5) +         # Set y-axis limits 
  # geom_line(color="blue", size=1) +
  stat_smooth(method=lm, size=0.5, color = 'black') +
  # sm_statCorr( color = "black") +
  labs(title="Second canonical variate", x=bquote('Musical variable'-U[2]), y= bquote('Physiologial variable'-V[2])) +
  theme_classic() + theme(text = element_text(size = 12, family = "Abadi MT"), plot.title = element_text(size = 12))

pY2 <- ggplot(mean_canY_track_2, aes(x = track, y =  Y_canonical_2[,1], fill = customisedColors)) +
  geom_point(shape = 21,size = 3,fill = mean_canY_track_2$customisedColors) + # Mean points
  geom_errorbar(aes(ymin = Y_canonical_2[,2], ymax = Y_canonical_2[,3]), width = 0.2) + # Confidence intervals\
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.25) + 
  scale_y_continuous(breaks = seq(-2, 2, by = 0.5)) +
  labs(x = "", y =bquote('Physiologial variable'-V[2]), title = "Separate pieces") +
  theme_classic()  +   theme(text = element_text(size = 12, family = "Abadi MT"), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 0, size=10), plot.title = element_text(size = 12)) 

plot_corr_sep_pieces <- grid.arrange(pX1, pX2, pY1, pY2, ncol = 2, widths = c(1, 2))


##### PERMUATION TEST ANALYSIS ######

#lists initialisation
match_indicesALL <- list()

combined_data <- cbind(X, Y)
proportion_covariance_explainedALL <- list()
correlation_ALL <- list()

#1000 interations of physiological data permutation (matrix Y)
for (iterBoost in 1:1000) {

  print(iterBoost)

  X_perm <- X
  Y_perm <- Y[sample(nrow(Y)),] #permutation of the Y matrix

  #### CCA analysis in one iteration (on permuted data)
  ccClassical_perm <- cc(X_perm,Y_perm) 

  canonical_vectors_X <- as.matrix(ccClassical_perm$xcoef)
  canonical_vectors_Y <- as.matrix(ccClassical_perm$ycoef)
  
  X_canonical <- as.matrix(X_perm) %*% canonical_vectors_X
  Y_canonical <- as.matrix(Y_perm) %*% canonical_vectors_Y

  #calculate covariance and canonical correlations coefficients in one iteration and append them to the joint array
  correlation_explained <- compute_correlation(X_canonical, Y_canonical)
  covariance_explained <- compute_covariance(X_canonical, Y_canonical)
  total_covariance <- sum(abs(covariance_explained))
  proportion_covariance_explained <- abs(covariance_explained) / total_covariance *100
  proportion_covariance_explainedALL <- rbind(proportion_covariance_explainedALL,proportion_covariance_explained)

  correlation_ALL <- rbind(correlation_ALL,correlation_explained)
}

 

#reshape the matrix with covariance values into an array with sizes nxm, where n is the number of iterations (1000) and m is the number of canonical functions calculated
covMat = matrix(as.numeric(proportion_covariance_explainedALL), nrow = 1000, ncol =9 )
mean_matrix_cov <- apply(covMat,2, FUN = median)
# Calculate standard deviation across matrices
sd_matrix_cov <- apply(covMat, 2, FUN = sd)
# Calculate standard error of the mean (SEM)
n <- nrow(covMat)
sem_matrix <- sd_matrix_cov / sqrt(n)
# Calculate t-value for 95% confidence level (two-tailed)
t_value <- qt(0.975, df = n - 1)
# Calculate 95% confidence interval (CI)
lower_ci_cov <- mean_matrix_cov - t_value * sem_matrix
upper_ci_cov <- mean_matrix_cov + t_value * sem_matrix

#reshape the matrix with correlation values into an array with sizes nxm, where n is the number of iterations (1000) and m is the number of canonical functions calculated
corrMat = matrix(as.numeric(correlation_ALL), nrow = 1000, ncol =9 )
mean_matrix_corr <- apply(corrMat,2, FUN = mean)
# Calculate standard deviation across matrices
sd_matrix_corr <- apply(corrMat, 2, FUN = sd)
# Calculate standard error of the mean (SEM)
n <- nrow(corrMat)
sem_matrix <- sd_matrix_corr / sqrt(n)
# Calculate t-value for 95% confidence level (two-tailed)
t_value <- qt(0.975, df = n - 1)
# Calculate 95% confidence interval (CI)
lower_ci_corr <- mean_matrix_corr - t_value * sem_matrix
upper_ci_corr <- mean_matrix_corr + t_value * sem_matrix

#collect the parameters for covariance and correlation into one dataframe
dfToPlot_temp <- data.frame(
  lower_ci_temp_cov = lower_ci_cov,
  upper_ci_temp_cov = upper_ci_cov,
  means_temp_cov = mean_matrix_cov,
  means_temp_cov_cumsum = cumsum(mean_matrix_cov),
  groups_cov = c(1:ncol(covMat)),
  
  lower_ci_temp_corr = lower_ci_corr,
  upper_ci_temp_corr = upper_ci_corr,
  means_temp_corr = mean_matrix_corr,
  sd_matrix_corr = sd_matrix_corr,
  groups_corr = c(1:ncol(corrMat))
  
)

#visualise mean covariance values (over all iterations) calculated on the permuted data
p1 <- ggplot(dfToPlot_temp, aes(x=groups_cov, y=means_temp_cov, color=groups_cov)) +
  geom_pointrange(aes(ymin=lower_ci_temp_cov, ymax=upper_ci_temp_cov)) + theme_classic() +
  geom_hline(yintercept=0, linetype="dashed", color = "black")
p1
mean(dfToPlot_temp$means_temp_cov)

#visualise mean canonical correlation values (over all iterations) calculated on the permuted data
p1 <- ggplot(dfToPlot_temp, aes(x=groups_cov, y=means_temp_corr, color=groups_cov)) +
  geom_pointrange(aes(ymin=lower_ci_temp_corr, ymax=upper_ci_temp_corr)) + theme_classic() +
  geom_hline(yintercept=0, linetype="dashed", color = "black")
p1
mean(dfToPlot_temp$means_temp_corr)


## visualise distributions (as histograms) of the canonical correlation coefficients (from the first two canonical functions) obtained from 1000 iterations of permutation the original data set.

#export path
dirPath = "./permutation_tests_ALL_resampling.pdf"

# calculate maximum limit for y-axis
ymax = max(c(max(hist(corrMat[,1],20)$counts),
             max(hist(corrMat[,2],20)$counts)))

#create a .pdf file for exporting the figure
pdf(dirPath, width = 12, height = 4)  # Specify in inches
par(mfrow = c(1, 2),
    cex.main = 1,   # Title font size
    cex.lab = 1,    # Axis label font size
    cex.axis = 1,mar = c(5, 5, 4, 2) + 0.5, oma = c(0, 0, 2, 0))

#histogram for the canonical correlations from the first canonical function calculated from all permutations 
variateNb = 1 #first canonical function
h1 <- hist(corrMat[,variateNb],20,ylim=c(0,ymax),xlab = "Correlation", ylab = "Counts", col = "skyblue", main = expression(1^{st}~ 'variate'))

# original correlation
orig_corr_coef <- correlation_orig[variateNb]

#calculate the p-value as the number of the iteration for whom the correlation is larger than the original one
Pperm = sum(corrMat[,variateNb]>orig_corr_coef)/1000
orig_corr_coef_rounded = round(orig_corr_coef,3)
Pperm_rounded = round(Pperm,3)
# Add a vertical dashed line at the specified position
# abline(v = orig_corr_coef, col = "red", lty = 2, lwd = 2)

 
# add text label with the original value of the canonical correlation
label_text1 <- bquote(Corr["orig"] == .(orig_corr_coef_rounded))

# add text label with the original value of the p-value
label_text2 <- bquote(P["perm"] == .(Pperm_rounded))

# Add text above the line at the appropriate height
text(x = max(h1$breaks*0.90), 
     y = ymax * 0.8,  # Place the text near the top of the histogram
     labels = label_text1, 
     col = "black", cex = 1,
     pos = 3)  # `pos = 3` positions text above

text(x = max(h1$breaks*0.90), 
     y = ymax * 0.65,  # Place the text near the top of the histogram
     labels = label_text2, 
     col = "black", cex = 1,
     pos = 3)  # `pos = 3` positions text above
h1

# create the histogram for the second canonical function
variateNb = 2
h2 <- hist(corrMat[,variateNb],20,ylim=c(0,ymax),xlab = "Correlation", ylab = "Counts", col = "plum", main = expression(2^{nd}~ 'variate'))

orig_corr_coef <- correlation_orig[variateNb]

Pperm = sum(corrMat[,variateNb]>orig_corr_coef)/1000
orig_corr_coef_rounded = round(orig_corr_coef,3)
Pperm_rounded = round(Pperm,3)
# Add a vertical dashed line at the specified position
# abline(v = orig_corr_coef, col = "red", lty = 2, lwd = 2)

label_text1 <- bquote(Corr["orig"] == .(orig_corr_coef_rounded))
label_text2 <- bquote(P["perm"] == .(Pperm_rounded))


# Add text above the line at the appropriate height
text(x = max(h2$breaks*0.90), 
     y = ymax * 0.8,   
     labels = label_text1, 
     col = "black", cex = 1,
     pos = 3)  

text(x = max(h2$breaks*0.90), 
     y = ymax * 0.65,   
     labels = label_text2, 
     col = "black", cex = 1,
     pos = 3)   
h2

par(mfrow = c(1, 1))
dev.off()





