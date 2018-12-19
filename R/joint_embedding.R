
# Author: Jesus Arroyo & Shangsi Wang
#' Function to perform joint embedding of graphs into d dimensions.
#' @param A a list of adjacency matrices with the same size n x n
#' @param d number of embedding dimensions
#' @param maxiter the maximum number of iterations on each step
#' @param Innitialize value for initialization. 1 = svd on the average. 2 = random, or an array of size n x d to set up to a specific value
#' @param verbose print output on each iteration. 0 (no output), 1 (output on outer loop) or 2 (output on inner loop)
#' @param large.and.sparse boolean option to indicate whether function for large and sparse graphs should be used (for best results, only set up to TRUE when elements in A are sparse and large matrices)
#' @return A list containing the result of the embeding
#' \item{objective}{A vector with the value of the objective function after each component has been found}
#' \item{lambda}{The loadings of the embedding}
#' \item{h}{The vector of latent positions of the embedding}
#' \item{iter}{Number of iterations for finding each component}
multidembed <- function(A,d,maxiter=20,Innitialize=1, verbose = 0, large.and.sparse = FALSE) {
  m <- length(A)
  n <- dim(A[[1]])[1]
  result <- list("objective" =rep(0,d),
                 "lambda" =matrix(0,m,d),
                 "h" =matrix(0,n,d),"iter"=rep(0,d))
  Gamma <- matrix(0, d, d)
  diag(Gamma) <- 1
  Psi <- matrix(0, d, m)
  H_c <- list()
  Resid<-A
  hprev <- NULL
  lambdaprev <- NULL
  for(k in 1:d){
    if(verbose > 0) cat("Iteration ", k,"\n")
    ##### Calculate one component
    if(!large.and.sparse) {
      ### Code for dense graphs ####################################################
      if(length(Innitialize)==1 && (Innitialize==1 | Innitialize ==2)){
        resultd <- onedembed(Resid,maxiter,Innitialize=Innitialize, verbose = verbose)
      } else {
        resultd <- onedembed(Resid,maxiter,Innitialize[,k], verbose = verbose)
      }
    }else {
      ### Code for dense graphs ####################################################
      if(length(Innitialize)==1 ){
        resultd <- onedembed_l(Resid, A, hprev, lambdaprev, 
                               maxiter,Innitialize=Innitialize, verbose = verbose)
      } else {
        resultd<-onedembed_l(Resid, A, hprev, lambdaprev, 
                             maxiter,Innitialize=Innitialize[,k], verbose = verbose)
      }
    }
    
    result$h[,k] <- as.vector(resultd$h)
    result$iter[k] <- resultd$iter
    
    ### Update lambda weights ###########################################
    if(k >1) {
      Gamma[1:(k-1), k] <- crossprod(result$h[,1:(k-1)], result$h[,k])^2
      Gamma[k, 1:(k-1)] <- Gamma[1:(k-1), k]
    }
    Psi[k, 1:m] <- sapply(A, function(a) 
      sum(t(result$h[,k]) %*% crossprod(a, result$h[,k])))
    result$lambda[,1:k] <- t(solve(Gamma[1:k, 1:k, drop=F], Psi[1:k,, drop = F]))
    #### Update residuals ###################################################
    hprev <- result$h[,1:k, drop = F]
    lambdaprev <- result$lambda[,1:k, drop = F]
    H_c[[k]] <- tcrossprod(resultd$h)
    # check if it is possible to get rid of this part for sparse
    for(i in 1:m){
      Resid[[i]]<-A[[i]] - Reduce(function(a,k) a + H_c[[k]]*result$lambda[i,k], 
                                  1:k ,init = matrix(0,n,n))
    }
    result$objective[k] <- sum(sapply(Resid, norm, type = "F")^2) 
  }
  return(result)
}


# Author: Jesus
# Wrapper for embedding method creating multiple random initialization values
# and running in parallel.
multidembed_random_parallel <- function(A, d, maxiter=20, rand_inits = 10, seed = 777, cl_numClusters = 12) {
  require(parallel)
  cl <- makeCluster(cl_numClusters)
  clusterExport(cl = cl, varlist = list("multidembed", "onedembed", "onedembed_l", "cptobj"))
  clusterExport(cl = cl, varlist = list("seed", "A", "d", "maxiter"), envir = environment())
  multiembed_random <- parLapply(cl = cl, 1:rand_inits, function(i) {
    set.seed(seed + i)
    require(Matrix)
    res <- multidembed(A, d,  Innitialize = 2)
    res
  })
  stopCluster(cl)
  objectives <- sapply(1:rand_inits, function(i) min(multiembed_random[[i]]$objective))
  
  return(list(best_result = multiembed_random[[which.min(objectives)]],
              results = multiembed_random,
              objectives = objectives))
}


# Author: Jesús Arroyo
#' Function to calculate the loadings for a given list of graphs A and 
#' @param A list of adjacency matrices of size n x n
#' @param h latent positions of JEG, of size n x d
#' @return Matrix of loadings for graphs in A
multiembed_test <- function(A, h) {
  Gamma <- crossprod(h)^2
  Psi <- sapply(A, function(a) 
    diag(t(h) %*% crossprod(a, h)))
  t(solve(Gamma, Psi))
}



#' Function to perform embedding into 1 dimension
#' @param R a list of graphs
#' @param maxiter the maximum number of iterations
#' @param Innitialize value for initialization. 1 = svd on the average. 2 = random
#' @return The result of one dimensional embedding of graphs
#' @export
onedembed <- function(R,maxiter=20,Innitialize=1, verbose = 0) {
  #browser()
  # initialize
  m<-length(R)
  n<-dim(R[[1]])[1]
  result<-list("objective" = 0,"lambda" =0,"h" =0,"iter"=0)
  Rbar<-matrix(0,n,n)
  Rbar <- Reduce("+", R)
  
  # default svd initialization
  #require(rARPACK)
  #s<-svds(Rbar,1)
  if(sum(abs(Rbar))/n^2<1e-10) {Innitialize = 2}else{
    require(irlba)
    s <- irlba(Rbar, nv = 1, nu = 1)
    h<-sqrt(s$d[1])*s$u    
    if(sum(abs(h))==0) {Innitialize =2}
  }
  
  # random initialize
  if (length(Innitialize)==1 && Innitialize==2){
    h <- as.matrix(rnorm(n))
  }
  
  # user input iniliatialization
  if(length(Innitialize)>1){
    h <- as.matrix(Innitialize)
  }
  
  h<-h/norm(h,type="F")
  
  lambda<-rep(0,m)
  lambda = sapply(R, function(a) sum(t(h) %*% crossprod(a, h)))
  
  obj_R_norm <- sum(sapply(R, norm, type = "F")^2)
  obj<-cptobj(R,lambda,h)
  objstart<-obj
  objpre<-Inf
  iter<-0
  while(iter<maxiter & ((objpre-obj)>(10^-7)*(1+abs(objstart)))){
    if(verbose > 1) cat("--iter ", iter, "--objpre-obj-", objpre-obj, "\n")
    iter<-iter+1
    # Gradient descent
    # Compute gradient of h
    Gradient<-matrix(0,n,1)
    for(i in 1:m){
      Gradient <- Gradient - (R[[i]]%*%h)*(lambda[i]*4) + h*(lambda[i]^2*4)
    }
    Gradient <- Gradient/m
    Gnorm <- norm(Gradient,type="F")
    stepsize<-1
    obj<-cptobj(R,lambda,h)
    objtmp<-Inf
    while(objtmp > obj - stepsize*Gnorm^2*0.01 & stepsize>10^-7 ){
      # backtracking line search 
      ht <- h - stepsize*Gradient
      htnorm <- norm(ht,type="F")
      ht <- ht/htnorm
      objtmp <- cptobj(R, lambda, ht)
      stepsize <- stepsize*0.5
    }
    
    if(stepsize<=10^-7 ) {
      objpre<-obj
    } else {
      objpre<-obj
      obj<-objtmp
      h<-ht
      # update lambda
      lambda <- sapply(R, function(a) sum(t(h) %*% crossprod(a, h)))
    }
  }
  
  obj <- cptobj(R,lambda,h)
  result$objective <- obj + obj_R_norm
  result$lambda <- lambda
  result$h <- h
  result$iter <- iter
  return(result)
}


#' Function to perform embedding into 1 dimension for sparse and large graphs
#' @param R a list of graphs
#' @param maxiter the maximum number of iterations
#' @param Innitialize value for initialization. 1 = svd on the average. 2 = random
#' @return The result of one dimensional embedding of graphs
#' @export
onedembed_l <- function(R, A, hprev=NULL, lambdaprev, maxiter=20,Innitialize=1, verbose = 0) {
  browser()
  # initialize
  m<-length(A)
  n<-dim(A[[1]])[1]
  result<-list("objective" = 0,"lambda" =0,"h" =0,"iter"=0)
  Rbar <- Reduce("+", R)
  
  # default svd initialization
  if(sum(abs(Rbar))/n^2<1e-10) {Innitialize = 2}else{
    require(irlba)
    s <- irlba(Rbar, nv = 1, nu = 1)
    h<-sqrt(s$d[1])*s$u    
    if(sum(abs(h))==0) {Innitialize =2}
  }
  
  # random initialize
  if(sum(abs(h))==0) {Innitialize =2}
  if (length(Innitialize)==1 && Innitialize==2){
    h <- as.matrix(rnorm(n))
  }
  
  # user input iniliatialization
  if(length(Innitialize)>1){
    h <- as.matrix(Innitialize)
  }
  
  h<-h/norm(h,type="F")
  
  lambda<-rep(0,m)
  if(is.null(hprev)) {
    lambda = sapply(A, function(a) sum(t(h) %*% crossprod(a, h))) 
  } else{
    hh <- crossprod(hprev, h)^2
    lambda = sapply(A, function(a) sum(t(h) %*% crossprod(a, h))) -
      as.vector(lambdaprev %*% hh) 
  }
  
  obj_A_norm <- sum(sapply(A, norm, type = "F")^2)
  obj<-cptobj(A,lambda,h)
  objstart<-obj + obj_A_norm
  objpre<-Inf
  iter<-0
  while(iter<maxiter & ((objpre-obj)>(10^-7)*(1+abs(objstart)))){
    if(verbose > 1) cat("--iter ", iter, "--objpre-obj-", objpre-obj, "\n")
    iter<-iter+1
    # Gradient descent
    # Compute gradient of h
    Gradient<-matrix(0,n,1)
    if(is.null(hprev)) {
      for(i in 1:m){
        Gradient <- Gradient - (A[[i]]%*%h)*(lambda[i]*4) + h*(lambda[i]^2*4)
      }
    }else {
      HHh <- hprev %*% diag(as.vector(crossprod(hprev, h)), nrow = ncol(hprev))
      for(i in 1:m){
        Gradient <- Gradient- (4*lambda[i])*(A[[i]] %*% h) +
          (4*lambda[i])* (HHh %*% lambdaprev[i,]) +
          h*(lambda[i]^2*4)
        
      }
    }
    
    Gradient <- Gradient/m
    Gnorm <- norm(Gradient,type="F")
    stepsize<-1
    obj<-cptobj(A,lambda,h) + obj_A_norm
    objtmp<-Inf
    while(objtmp > obj - stepsize*Gnorm^2*0.01 & stepsize>10^-7 ){
      # backtracking line search 
      ht <- h - stepsize*Gradient
      htnorm <- norm(ht,type="F")
      ht <- ht/htnorm
      objtmp <- cptobj(A, lambda, ht) + obj_A_norm
      stepsize <- stepsize*0.5
    }
    
    if(stepsize<=10^-7 ) {
      objpre<-obj
    } else {
      objpre<-obj
      obj<-objtmp
      h<-ht
      # update lambda
      if(is.null(hprev)) {
        lambda = sapply(A, function(a) sum(t(h) %*% crossprod(a, h))) 
      } else{
        hh <- crossprod(hprev, h)^2
        lambda = sapply(A, function(a) sum(t(h) %*% crossprod(a, h))) -
          as.vector(lambdaprev %*% hh) 
      }
    }
  }
  
  obj <- cptobj(A,lambda,h) + obj_A_norm
  result$objective <- obj
  result$lambda <- lambda
  result$h <- h
  result$iter <- iter
  return(result)
}



## calculate objective function
cptobj <- function(A,lambda,h) {
  m<-length(A)
  u <- sum(-2*(sapply(A, function(a) as.double(t(h) %*% crossprod(a, h))) * lambda) + (lambda^2))
  return(as.double(u))
}

