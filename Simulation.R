
# 1. options / packages ---------------------------------------------------

options(scipen = 100)
if (!require(plyr)) install.packages('plyr'); require(plyr)
if (!require(dplyr)) install.packages('dplyr'); require(dplyr)
if (!require(Matrix)) install.packages('matrix'); require(Matrix)
if (!require(car)) install.packages('car'); require(car)
if (!require(MASS)) install.packages('MASS'); require(MASS)
if (!require(tseries)) install.packages('tseries'); require(tseries)


# 2. Utils ----------------------------------------------------------------

### 블록의 개수에 맞게 W행렬 생성하는 함수

getBlockW <- function(N, Nblock) {
  
  if (N %% Nblock == 0) {
    
    isDiagList = rep(list(matrix(1, 
                                 nrow = N/Nblock, 
                                 ncol = N/Nblock)), Nblock)  
    mList = rep(list(matrix(rbinom((N / Nblock) ^ 2, size = 1, prob = 0.3 * N^{-0.3}), 
                            nrow = N / Nblock, 
                            ncol = N / Nblock)), 
                Nblock)
  } 
  else {
    
    isDiagList = rep(list(matrix(1, 
                                 nrow = floor(N/Nblock), 
                                 ncol = floor(N/Nblock))),
                     Nblock - 1)
    isDiagList[[length(Nblock)]] = matrix(1, 
                                          nrow = N%%Nblock, 
                                          ncol = N%%Nblock)
    
    mList = rep(list(matrix(rbinom(floor(N / Nblock) ^ 2, size = 1, prob = 0.3 * N ^ {-0.3}), 
                            nrow = floor(N/Nblock), 
                            ncol = floor(N/Nblock))), 
                Nblock - 1)
    
    mList[[Nblock]] = matrix(rbinom(floor(N / Nblock) ^ 2, size = 1, prob = 0.3 * N ^ {-0.3}), 
                             nrow = floor(N / Nblock), 
                             ncol = floor(N / Nblock))
  }
  
  isDiag = bdiag(isDiagList)  
  offDiag = which(isDiag == 0, arr.ind = T)  
  
  mList = lapply(mList, function(M) {
    
    ind = which(rowSums(M) == 0)
    if (length(ind) > 0) 
      
      M[cbind(ind, sample(1:nrow(M), length(ind)))] = 1
    
    return(M)
  })
  
  bA = bdiag(mList)
  bA[offDiag] = rbinom(nrow(offDiag), size = 1, prob = 0.3/N)
  bA = as.matrix(bA)
  upperInd = which(upper.tri(bA), arr.ind = T)
  
  bA[upperInd[, 2:1]] = bA[upper.tri(bA)]
  diag(bA) = 0
  
  ind = which(rowSums(bA) == 0)  
  for (i in ind) {
    bA[i, sample(setdiff(1:N, i), 3)] = 1
  }
  
  bA = as(bA, "dgCMatrix")
  W = bA/rowSums(bA)
  
  return(as.matrix(W))
  
}

### OLS 추정함수

betaOLS <- function(Ymat, W, Z) {
  
  Ymat1 <- W %*% Ymat  
  
  Time <- ncol(Ymat) - 1  
  
  X <- cbind(rep(1, nrow(Ymat) * Time), 
             as.vector(Ymat1[, -ncol(Ymat)]),
             as.vector(Ymat[, -ncol(Ymat)]), 
             do.call("rbind", rep(list(Z), Time)))
  
  invXX <- solve(crossprod(X))  
  
  Yvec <- as.vector(Ymat[, -1])  
  
  thetaEst <- invXX %*% colSums(X * Yvec)  
  
  sigmaHat2 <- mean((Yvec - X %*% thetaEst)^2) 
  
  covHat <- invXX * sigmaHat2
  
  return(list(theta = thetaEst, 
              sigmaHat = sqrt(sigmaHat2),
              covHat = covHat,
              Ymat = Ymat,
              X = X)) 
}



### SB 추정함수

betaSB <- function(result, W, Z){
  
  Ymat <- result$Ymat
  theta <- result$theta
  
  error <- result$X %*% result$theta - result$Ymat[, -1] %>% as.vector()
  
  Ymat1 <- W %*% Ymat
  Time <- ncol(Ymat) - 1
  
  index <- tsbootstrap(1:Time, b = 1 / (0.05 * (Time / 100) ^ {-1/3}))
  Ymat_star <- Ymat[, index]
  Ymat1_star <- Ymat1[, index]
  
  X_star <- cbind(rep(1, nrow(Ymat) * Time),
                  as.vector(Ymat1_star[, -ncol(Ymat)]),
                  as.vector(Ymat_star[, -ncol(Ymat)]), 
                  do.call("rbind", rep(list(Z), Time)))
  
  Ytilda <- matrix(X_star %*% theta + sample(error , 1),
                   ncol = Time)
  
  SB_result <- betaOLS(Ytilda, W, Z)
  
  return(SB_result)
  
}

### 다변량 정규분포에 따라서 Y_1 ~ Y_T 생성하는 함수
### 모수인 평균벡터와 공분산행렬은 GetMu, GetGamma0Approx 함수를 통해 추정

getY <- function(Y0, Beta0, Beta, W, Time, sig = 1) {
  
  Ymat = matrix(0, nrow = length(Beta0), ncol = Time + 1) 
  Ymat[, 1] = as.vector(Y0)  
  
  for (i in 1:Time) {
    Ymat[, i + 1] = as.vector(Beta0 + Beta[1] * W %*% Ymat[, i] + Beta[2] * 
                                Ymat[, i] + rnorm(length(Beta0), sd = sig))  
  }
  
  return(Ymat)
}

getMu <- function(G, Beta0) {
  mu = solve(Diagonal(length(Beta0)) - G) %*% Beta0 
  return(mu)
}

getGamma0Approx <- function(G, sigma) {
  Gamma0 = (Diagonal(nrow(G)) + tcrossprod(G) + tcrossprod(G %*% G) + tcrossprod(G %*% G %*% G)) * sigma^2
  return(Gamma0)
}

simu.data <- function(W, beta0, Beta, Time, G, Z, sig = 1) {
  
  mu <- getMu(G, beta0 + Z %*% gamma0)
  Gamma0 <- getGamma0Approx(G, sigma = sig)
  Y0 <- mvrnorm(n = 1, mu = mu, Sigma = Gamma0)
  
  Ymat = getY(Y0 = Y0, beta0 + Z %*% gamma0, Beta, W, Time, sig = 1)
  
  return(Ymat)
}

### 실험 반복시키는 함수

fin_function <- function(beta, gamma0, W, G, Z, Ymat, Nrep){
  
  estimation <- list()
  
  for(i in 1:Nrep){
    
    W <- getBlockW(N = Nsize, Nblock = Nblock)
    G <- beta[2] * W + beta[3] * diag(1, nrow(W))
    Z <- mvrnorm(n = Nsize, mu = rep(0, nrow(ZSigma)), Sigma = ZSigma)
    Ymat = simu.data(W = W, 
                     beta0 = beta[1], 
                     Beta = beta[2:3], 
                     Time = Time, 
                     G = G, 
                     Z = Z, 
                     sig = 1)
    
    Thetaest <- betaOLS(Ymat, W, Z)
    Thetaest_SB <- betaSB(result = Thetaest, W = W, Z = Z)
    
    estimation[[i]] <- list(OLS = Thetaest$theta,
                            SB = Thetaest_SB$theta)
    
    if (i %% 100 == 1 | i == 1){
      print(paste0(i, '번 째 실험 중입니다.'))
    } 
    else if (i == Nrep){
      print(paste0(Nrep, '번의 실험이 모두 종료되었습니다.'))
    }
  }
  
  return(estimation)
  
}

getresult <- function(fin){
  
  result2 <- list()
  
  for(i in 1:8){
    
    index <- seq(from = i, to = nrow(fin), by = 8)
    result2[[i]] <- fin[index, ]
  }
  
  result2 <- result2 %>% 
    cbind.data.frame() %>% 
    `names<-`(nn)
}


### Average Length 구하는 함수 

AV <- function(result){
  
  low <- (2 * rep(params, each = 2)) - apply(result, 2, function(x) quantile(x, 0.95, type = 1))
  high <- (2 * rep(params, each = 2)) - apply(result, 2, function(x) quantile(x, 0.05, type = 1))
  (high - low) %>% round(4)
  
}

# 3. Simulation -----------------------------------------------------------

ZSigma <- 0.5^abs(outer(1:5, 1:5, "-")) 
beta <- c(0, 0.1, -0.2)
gamma0 <- c(-0.5, 0.3, 0.8, 0, 0)
params <- c(beta, gamma0)
Nrep <- 1000

nn <- c(paste0(rep(c('NA_', 'SB_'), 2), rep(paste0('beta', 1:3), each = 2)),
        paste0(rep(c('NA_', 'SB_'), 2), rep(paste0('gamma', 1:5), each = 2)))

Nsize <- 200; Nblock <- 20; Time <- 200
result <- fin_function(beta, gamma0, W, G, Z, Ymat, Nrep) %>% 
  bind_rows() %>% 
  getresult()
