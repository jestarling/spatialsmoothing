bisquare.basis <- function(coord, scale, level1, level2, level3){
  
  n <- dim(coord)[1]
  # level_i is the matrix with r_i*2 centers
  r <- dim(level1)[1] + dim(level2)[1] + dim(level3)[1]
  S <- array(NA, dim = c(n, r))
  
  levels <- list(level1, level2, level3)
  
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