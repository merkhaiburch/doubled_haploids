library(ggfortify)

# Plot results from pca
autoplot(prcomp(df.pca))

# Color dots based on env.lin from PCA plot
autoplot(prcomp(df.pca), data = df.dh, colour = 'env.lin', loadings = TRUE)

# Plot PCA w/colored env.lin and place eigenvectors
# Look how tassel branches stands out
autoplot(prcomp(df.pca), data = df.dh, colour = 'env.lin',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)

# lineage.generation
autoplot(prcomp(df.pca), data = df.dh, colour = 'lineage.generation',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)

# lineage
autoplot(prcomp(df.pca), data = df.dh, colour = 'lineage',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)

# block
autoplot(prcomp(df.pca), data = df.dh, colour = 'Block',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)



# factor analysis
d.factanal <- factanal(df.pca, factors = 3, scores = 'regression')

# Plot factor analysis
autoplot(d.factanal, data = df.dh, colour = 'env.lin')

# plot with eigenvectors
autoplot(d.factanal, label = TRUE, label.size = 3,
         loadings = TRUE, loadings.label = TRUE, loadings.label.size  = 3)

# Plotting k-means
set.seed(1)
autoplot(kmeans(df.pca, 15), data = df.pca, label = TRUE)

# Plotting cluster package
library(cluster)
autoplot(clara(df.pca, 5))
autoplot(fanny(df.pca, 5), frame = TRUE)

# Plotting Local Fisher Discriminant Analysis with {lfda} package
library(lfda)

model <- lfda(df.pca[-5], df.dh[,7], r = 5, metric = "plain")
autoplot(model, data = df.dh, frame = TRUE, frame.colour = 'lineage.generation')










