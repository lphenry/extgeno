#This script simulates a quantitative phenotype that is affected by both host genome and microbiome, 
#under different heritabilities and degrees of microbiome transmission. 
#It then calculates the selection response R for different selection intensities S
#For simplicity, we simulate nuclear families with 3 F1s each. This simulation ran in ~1.7 days on a normal desktop computer
library(mvnfast)
library(MASS)

n <- 10000 #Population size, including both parents and F1s
h2s.host <- c(1e-4, .2, .4) 
h2s.micro <- c(1e-4, .2, .4)
t.df <- c(0,1,3,5,7) #df for the distribution of t (transmission parameter)

S <- seq(.1, 1.3, .2) #Selection intensities
n_sim <- 20 #Number of simulations per parameter combination
grid <- expand.grid(sim = 1:n_sim, h2s.host = h2s.host, h2s.micro = h2s.micro, t.df = t.df, S = S)
simRes <- data.frame(grid, R.obs = NA, S.obs = NA)

for (i in 1:nrow(grid)) {
  h2.host <- grid$h2s.host[i]
  h2.micro <- grid$h2s.micro[i]
  df <- grid$t.df[i]
  S <- grid$S[i]
  
  G <- kinship(N = n) #Generates kinship matrix
  
  #Host phenotypic effect
  g.host <- rmvn(1, rep(0, n), h2.host * G) #genetic effects
  
  #Microbiome phenotypic effect
  t <- 1/(rchisq(n = n/5, df = df) + 1) #Transmission parameters, one per family. 0 < t =< 1
  t.mat <- kronecker(diag(t), matrix(1, nrow = 5, ncol = 5))
  diag(t.mat) <- 1
  g.micro <- rmvn(1, rep(0, n), h2.micro * t.mat * G)
  
  #Phenotype
  y <- g.host[1,] + g.micro[1,] + rnorm(n = n, mean = 0, sd = sqrt(1 - h2.host - h2.micro)) 
  names(y) <- 1:length(y)
  
  
  #Select parents such that mean(selected parents) =~ S
  dads <- seq(from = 1, to = n - 4, by = 5) # dads
  moms <- seq(from = 2, to = n - 3, by = 5) # moms
  y.parents <- y[sort(c(dads, moms))]
  t <- -5
  while (mean(y.parents[y.parents > t]) < S) { 
    t <- t+.01
  }
  parents <- which(y.parents > t) #names(parents) are the idx in y. parents are the idx in y.parents
  #Remove "single parents"
  tmp <- as.numeric(names(parents))
  parents <- parents[tmp %in% (tmp + 1) | tmp %in% (tmp - 1)]
  
  
  #Calculate mean of F1s from selected parents = selection response R
  if(length(parents) > 6){
    #Pick F1s
    parents.yidx <- as.numeric(names(parents))
    moms <- parents.yidx[seq(from = 2, to = length(parents), by = 2)]
    f1.1 <- moms+1
    f1.2 <- moms+2
    f1.3 <- moms+3
    y.f1s <- y[c(f1.1, f1.2, f1.3)]
    
    #Save R and S
    simRes$R.obs[i] <- mean(y.f1s) #Observed R
    simRes$S.obs[i] <- mean(y.parents[parents]) #Observed S
  }
  
  message(i, '/', nrow(grid))
}
save(simRes, file = 'simulation.RData')

