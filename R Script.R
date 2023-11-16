# Visualizing Gene Expression Data with GEO samples

# Load necessary libraries
library(GEOquery) # For querying GEO databases
library(limma)    # For linear models for microarray data
library(umap)     # For generating UMAP plots


# Fetching the dataset from GEO
gset <- getGEO("GSE5819", GSEMatrix = TRUE, getGPL = FALSE)
# Selecting the appropriate platform ID if multiple are present
if (length(gset) > 1) {
  idx <- grep("GPL4315", attr(gset, "names"))
} else {
  idx <- 1
}
gset <- gset[[idx]]

# Extracting expression data
ex <- exprs(gset)

# Normalizing data with log2 transformation if certain conditions are met
qx <- as.numeric(quantile(ex, probs = c(0.0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) {
  ex[ex <= 0] <- NaN
  ex <- log2(ex)
}

# Generating a box-and-whisker plot for expression data
dev.new(width = 3 + ncol(gset) / 6, height = 5)
par(mar = c(7, 4, 2, 1))
title <- paste("GSE5819", "/", annotation(gset), sep = "")
boxplot(ex, boxwex = 0.7, notch = TRUE, main = title, outline = FALSE, las = 2)
dev.off()

# Creating an expression value distribution plot
par(mar = c(4, 4, 2, 1))
title <- paste("GSE5819", "/", annotation(gset), " value distribution", sep = "")
plotDensities(ex, main = title, legend = FALSE)

# Plotting mean-variance trend of the data
ex <- na.omit(ex) # Removing missing values
plotSA(lmFit(ex), main = "Mean variance trend, GSE5819")

# Performing UMAP for dimensionality reduction and visualization
ex <- ex[!duplicated(ex), ] # Removing duplicate rows
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
plot(ump$layout, main = "UMAP plot, nbrs=15", xlab = "", ylab = "", pch = 20, cex = 1.5)
library("maptools") # For better label placement
pointLabel(ump$layout, labels = rownames(ump$layout), method = "SANN", cex = 0.6)
