library(tidyverse)

# create a matrix with 10 rows, 10 columns and random values
mat <- matrix(rnorm(100), nrow = 10, ncol = 10)

# insert random NAs in 30% of the cells
mat[sample(1:100, 30)] <- NA
colnames(mat) <- paste0("sample", 1:10)
# comput column-wise correlation of features using complete pairwise observations and plot it
cor_mat <- cor(mat, use = "pairwise.complete.obs")
# set the diagonal and upper triangle to NA
cor_mat[upper.tri(cor_mat)] <- NA
# plot as a heatmap without clustering
heatmap(cor_mat, Rowv = NA, Colv = NA, symm = TRUE, scale = "none", margins = c(5, 5))

# binarize matrix: 1s for values and 0s for NAs
mat_bin <- ifelse(is.na(mat), 0, 1)
# compute the crossproduct of the binarized matrix
mat_crossprod <- crossprod(mat_bin)
# remove diagonal and upper triangle
mat_crossprod[upper.tri(mat_crossprod)] <- NA
# plot as a heatmap without clustering and indicating number of observations within cells
heatmap(mat_crossprod, Rowv = NA, Colv = NA, symm = TRUE, scale = "none", margins = c(5, 5), col = rev(heat.colors(10)))
