# X:list
# K:the number of tissues
# G:vector,the number of cell types in each tissue

SRMDG = function(X, K, G, lambda, alpha, beta, model = "Gaussian", weights = "equal", penalize.diagonal = FALSE){
  
  p = dim(X[[1]])[2]
  
  cellType_num = sum(G)   #The number of cell type in all tissues
  n = rep(0, cellType_num)
  for (g in 1:cellType_num){
    n[g] =  dim(X[[g]])[1] #the number of samples in each cell type
  }
  
  S = matrix(list(), 1, cellType_num)
  
  if (model == "Gaussian"){
    for (g in 1:cellType_num){
      S[[g]] = Corr(X[[g]], method = "pearson")
    }
  }
  
  if (model == "nonparanormal_miss"){
    for (g in 1:cellType_num){
      tmp = X[[g]]
      tmp[tmp == 0] = NA
      S[[g]] = Corr(tmp, method = "kendall")
    }
  }
  
  if (weights=="equal"){
    n = rep(1,cellType_num)
  }else if (weights=="sample.size"){
    n = n
  }else{
    n = weights
  }
  
  
  result = SRMDG.admm(S, K, G, lambda, alpha, beta, n = n, penalize.diagonal = penalize.diagonal)
  
  result$Omega.bar = matrix(list(), 1, G[1])
  result$D.bar = matrix(list(), 1, G[1])
  result$Z.bar = result$Z.hat[[1]]!=0
  diag(result$Z.bar) <- 0
  for (g in 1:G[1]){
    result$Omega.bar[[g]] =  result$Omega.hat[[g]]!=0
    diag(result$Omega.bar[[g]]) <- 0
    result$D.bar[[g]] =  (result$Omega.hat[[g]]!=0)&(result$Z.bar==0)
    diag(result$D.bar[[g]]) <- 0
    
  }
  
  result
}


SRMDG.admm = function(S, K, G, lambda, alpha, beta, n = NULL, penalize.diagonal = FALSE, epsilon = 1e-5,
                     maxiter = 500, rho = 0.1,  rho.incr = 1.2, rho.max = 1e10){

  p = dim(S[[1]])[1]
  cellType_num = sum(G)
  
  if(is.null(n)){
    n = rep(1, cellType_num)
  }
  
  # initialize:
  Omega = matrix(list(), 1, cellType_num)
  Z = matrix(list(), 1, K)
  D = matrix(list(), 1, cellType_num)
  Q = matrix(list(), 1, cellType_num)
  
  A = diag(p)
  
  for (k in 1:K){
    Z[[k]] = diag(p)
  }
  
  for (g in 1:cellType_num){
    D[[g]] = diag(p)
    Omega[[g]] = diag(p)
    Q[[g]] = matrix(0, p, p)
  }
  
  
  for (i in 1:maxiter){
    
    Omega_prev = Omega
    sum_G = cumsum(G)
    k = 1
    for (g in 1:cellType_num){
      if (g > sum_G[1]){
        k = k + 1
        sum_G = sum_G[-1]
      }
      Omega[[g]] = expand(A + Z[[k]] + D[[g]] - (1/rho)*(n[g]*S[[g]] + Q[[g]]), rho, n[g])
    }
    
    
    sum_G = cumsum(G)
    k = 1
    B = matrix(list(), cellType_num, 1)
    for (g in 1:cellType_num){
      if (g > sum_G[1]){
        k = k + 1
        sum_G = sum_G[-1]
      }
      B[[g]] = Q[[g]]/rho + Omega[[g]] - Z[[k]] - D[[g]]
    }
    sum_G = cumsum(G)
    B = Reduce("+",B)
   # print(B)
    B = B/sum_G[length(sum_G)]
    A = group_lasso(B, lambda*alpha/rho, penalize.diagonal)
    
    
    M = matrix(list(), K, 1)
    l = 0
    for (k in 1:K){
      M[[k]] = matrix(0, p, p)
      for (g in 1:G[k]){
        if (k == 1){
          l = g
        }
        else{
          l = sum_G[k-1] + g
        }
        M[[k]] = M[[k]] + (Q[[l]]/rho + Omega[[l]] - A - D[[l]])
      }
      M[[k]] = M[[k]]/G[k]
      Z[[k]] = group_lasso(M[[k]], lambda*beta/rho, penalize.diagonal)
    }
    
    
    sum_G = cumsum(G)
    k = 1
    for (g in 1:cellType_num){
      if (g > sum_G[1]){
        k = k + 1
        sum_G = sum_G[-1]
      }
      D[[g]] = group_lasso((Q[[g]]/rho + Omega[[g]] - A - Z[[k]]), (lambda*(1 - alpha - beta))/rho, penalize.diagonal)
    }
    
    
    # updata dual variables Q
    sum_G = cumsum(G)
    k = 1
    for (g in 1:cellType_num){
      if (g > sum_G[1]){
        k = k + 1
        sum_G = sum_G[-1]
      }
      Q[[g]] =  Q[[g]] + rho*(Omega[[g]] - (A + Z[[k]] + D[[g]]))
    }
    
    # Check the convergence condition;
    diff_value_1 = 0
    diff_value_2 = 0
    norm_value = 0
    sum_G = cumsum(G)
    k = 1
    for (g in 1:cellType_num){
      if (g > sum_G[1]){
        k = k + 1
        sum_G = sum_G[-1]
      }
      diff_value_1 = diff_value_1 + sum(abs(Omega[[g]] - Omega_prev[[g]]))
      diff_value_2 = diff_value_2 + sum(abs(Omega[[g]] - (A + Z[[k]] + D[[g]])))
      norm_value  = norm_value + sum(abs(Omega[[g]]))
    }
    if (max(diff_value_1, diff_value_2) <= norm_value*epsilon){
      break
    }
    
    rho = min(rho*rho.incr,rho.max)
    
    
  }
  
  sum_G = cumsum(G)
  k = 1
  for (g in 1:cellType_num){
    if (g > sum_G[1]){
      k = k + 1
      sum_G = sum_G[-1]
    }
    Omega[[g]] = A + Z[[k]] + D[[g]]
  }
  
  out = list(Omega.hat = Omega, A.hat = A, Z.hat = Z, D.hat = D)
}






# Compute sample covariance matrices
Corr <- function(X, method = "pearson"){
  n = dim(X)[1]
  if(method=="pearson"){
    S = (n-1)*stats::cov(X, method = "pearson")/n
  }
  
  
  if (method == "kendall"){
    S = stats::cor(X, use = "pairwise.complete.obs", method = "kendall")
    S = sin(S*pi/2)
    S[is.na(S)] = 0
    diag(S) <- 1
    S = as.matrix(Matrix::nearPD(S, corr = TRUE)$mat)
  }
  S
}


# Expand operator
expand = function(A, rho, n){
  edecomp <- eigen(A)
  D <- edecomp$values
  U <- edecomp$vectors
  D2 <- 0.5*(D + sqrt(D^2 + 4*n/rho ))
  Omega <- U %*% diag(D2) %*% t(U)
  Omega
}


group_lasso =  function(A, lambda, penalize.diagonal=FALSE){
  m = dim(A)[1]
  n = dim(A)[2]
  X = matrix(0, m, n)
  for (i in 1:m){
    for (j in 1:n){
      X[i,j] = sign(A[i,j]) * max(abs(A[i,j]) - lambda, 0)
    }
  }
  if(!penalize.diagonal){
        diag(X) <- diag(A)
      }
  X
}

