

# Pick the coordinates obtained for UT campus
fine.grid <- 40
x <- seq(-97.74280, -97.72277, length.out = fine.grid)
y <- seq(30.27794, 30.29534, length.out = fine.grid)
coord <- as.matrix(expand.grid(x, y))

rangex <- range(x)[2] - range(x)[1]
rangey <- range(y)[2] - range(y)[1]

# Create the centers of the three scale levels
level1 <- as.matrix(expand.grid(seq(min(x) + rangex/5, max(x) - rangex/5, length.out = 3), seq(min(y) + rangey/5, max(y) - rangey/5, length.out = 3)))
level2 <- as.matrix(expand.grid(seq(min(x) + rangex/10, max(x) - rangex/10, length.out = 6), seq(min(y) + rangey/10, max(y) - rangey/10, length.out = 6)))
level3 <- as.matrix(expand.grid(seq(min(x) + rangex/20, max(x) - rangex/20, length.out = 9), seq(min(y) + rangey/20, max(y) - rangey/20, length.out = 9)))

r1 <- dim(level1)[1]
r2 <- dim(level2)[1]
r3 <- dim(level3)[1]
# Plot the centers on the grid of data points
plot(coord, type = 'p', pch = 16, cex = 0.6, asp = 1)
points(level1, type = 'p', pch = 16, cex = 0.6, asp = 1, col = 'red')
points(level2, type = 'p', pch = 16, cex = 0.6, asp = 1, col = 'green')
points(level3, type = 'p', pch = 16, cex = 0.6, asp = 1, col = 'blue')

  
bisquare.basis <- function(coord, level1, level2, level3){
  
  n <- dim(coord)[1]
  # level_i is the matrix with r_i*2 centers
  r <- dim(level1)[1] + dim(level2)[1] + dim(level3)[1]
  S <- array(NA, dim = c(n, r))
  
  levels <- list(level1, level2, level3)

  scale <- c(5E-3, 2E-3, 1E-3)
    
  count <- 0
  for (i in 1:3){
    for (j in 1:dim(levels[[i]])[1]){
    count <- count+1
    h <- sqrt(colSums((levels[[i]][j,] - t(coord))^2))
    s <- (1-(h/scale[i])^2)^2
    s[h > scale[i]] <- 0
    S[,count] <- s
    }
  }
  return (S)
}

S <- bisquare.basis(coord, level1, level2, level3)


z <- rowSums(S[,1:r1])
z <- S[,2]
z <- matrix(z, fine.grid, fine.grid, byrow = T)

# Define colors
nrz <- nrow(z)
ncz <- ncol(z)
jet.colors <- colorRampPalette(c("blue", "red"))
nbcol <- 100
color <- jet.colors(nbcol)
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
facetcol <- cut(zfacet , nbcol)

# # Plot
# persp(x, y, z, col = color[facetcol], phi = 45, theta = 45)
# Contour plot
filled.contour(x, y, z, color.palette = heat.colors, asp = 1, plot.axes={points(level1)})


z <- rowSums(S[,(r1+1):(r1+r2)])
z <- matrix(z, fine.grid, fine.grid, byrow = T)

# Define colors
nrz <- nrow(z)
ncz <- ncol(z)
jet.colors <- colorRampPalette(c("blue", "red"))
nbcol <- 100
color <- jet.colors(nbcol)
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
facetcol <- cut(zfacet , nbcol)

# # Plot
# persp(x, y, z, col = color[facetcol], phi = 45, theta = 45)
# Contour plot
filled.contour(x, y, z, color.palette = heat.colors, asp = 1, plot.axes={points(level2)})


z <- rowSums(S[,(r1+r2+1):(r1+r2+r3)])
z <- matrix(z, fine.grid, fine.grid, byrow = T)

# Define colors
nrz <- nrow(z)
ncz <- ncol(z)
jet.colors <- colorRampPalette(c("blue", "red"))
nbcol <- 100
color <- jet.colors(nbcol)
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
facetcol <- cut(zfacet , nbcol)

# # Plot
# persp(x, y, z, col = color[facetcol], phi = 45, theta = 45)
# Contour plot
filled.contour(x, y, z, color.palette = heat.colors, asp = 1, plot.axes={points(level3)})







# #3dplot
# open3d()
# persp3d(x, y, z, aspect=c(1,1,1), col = "lightgreen", alpha=0.4, xlab = "X", ylab = "Y", zlab = "")

