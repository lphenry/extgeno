kinship <- function(N, L = 5)
{
  if(!require(Matrix))
    stop('Could not load package Matrix')
  
  stopifnot(!(N %% L)) # family size = L
  Nl <- N / L
  
  # kinship of a single family
  K0 <- matrix(0.5, L, L)
  diag(K0) <- 1
  K0[1, 2] <- 0
  K0[2, 1] <- 0
  
  K <- kronecker(Diagonal(Nl), K0)
  
  ids <- seq(N)
  rownames(K) <- ids
  colnames(K) <- ids
  
  return(K)
}

