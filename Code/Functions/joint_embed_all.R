n<-50
m<-20
A<-list()
for (i in 1:m){
Ai<-matrix(runif(n^2,0,1),n,n)
A[[i]]<-(Ai+t(Ai))/2
A[[i]][1:10,1:10]<-A[[i]][1:10,1:10]+1
}
#result1<-multidembed(A,3,50)

cptobj_all <- function(A,lambda,h) {
  m<-length(A)
  obj<-0
  for(i in 1:m){
    obj<-obj+norm(A[[i]]-h%*%diag(lambda[i,])%*%t(h),type="F")^2
  }
  return(obj)
}



multidembed_all <- function(A,d,maxiter=20,h0=NULL) {
	m<-length(A)
	n<-dim(A[[1]])[1]
	result<-list("objective" =rep(0,d),"lambda" =matrix(0,m,d),"h" =matrix(0,n,d),"iter"=rep(0,d))
	
	
	###Innitialize obj, h, lambda, gradient
	Abar<-matrix(0,n,n)
	for(i in 1:m){
	  Abar<-Abar+A[[i]]
	}
	
	s<-svd(Abar,d,d)
	h<-s$u
	if (! is.null(h0)) {h <- h0
	hnorm<-apply(h, 2, function(x){sqrt(sum(x^2))})
	h<-h%*% diag(1/hnorm)}
	
	
	XX <- matrix(0,d,d)
	hh <- list()
	for(i in 1:d){
	  hh[[i]] <- h[,i]%*% t(h[,i])
	}
	
	for(i in 1:d){
	  for(j in 1:d){
	    XX[i,j] <-  sum(hh[[i]] * hh[[j]])
	  }
	}
	
	XY <- matrix(0,d,1)
	lambda <- matrix(0,m,d)
	for(j in 1:m){
	for(i in 1:d){
	  XY[i,1] <- sum(hh[[i]] * A[[j]])
	}
	 lambda[j,] <- solve(XX) %*% XY
	}
	
	obj<-cptobj_all(A,lambda,h)
	objstart<-obj
	objpre<-Inf
	
	iter<-0
	objv<-c(obj)
	
	while(iter<maxiter & ((objpre-obj)>(10^-7)*(1+objstart))){
	Gradiant<-matrix(0,n,d)
	for(i in 1:m){
	  Gradiant<-Gradiant-A[[i]]%*%h%*%diag(lambda[i,])*4+h%*%diag(lambda[i,]^2)*4
	}
	Gradiant<-Gradiant/m
	#Gnorm <- apply(Gradiant, 2, function(x){sqrt(sum(x^2))})
	#Gradiant<-Gradiant %*% diag(1 / Gnorm )
	#angle <- diag((t(Gradiant) %*% h)/Gnorm)
	#Gradiant<-Gradiant %*% diag(1 - angle^2 )
	Gnorm<-norm(Gradiant,type='F')
	print(Gnorm)
	
	
	objtmp<-Inf
	stepsize<-1
	while(objtmp>obj-stepsize*2*Gnorm^2*0.01 & stepsize>10^-7 ){
	  ht<-h-Gradiant*stepsize
	  htnorm<-apply(ht, 2, function(x){sqrt(sum(x^2))})
	  ht<-ht%*% diag(1/htnorm)
	  objtmp<-cptobj_all(A,lambda,ht)
	  stepsize<-stepsize*0.5
	}

	h<-ht
	for(i in 1:d){
	  hh[[i]] <- h[,i]%*% t(h[,i])
	}
	
	for(i in 1:d){
	  for(j in 1:d){
	    XX[i,j] <-  sum(hh[[i]] * hh[[j]])
	  }
	}
	
	lambda <- matrix(0,m,d)
	for(j in 1:m){
	  for(i in 1:d){
	    XY[i,1] <- sum(hh[[i]] * A[[j]])
	  }
	  lambda[j,] <- solve(XX) %*% XY
	}
	
	objpre <- obj
	obj <- cptobj_all(A,lambda,h)
	#print(obj)
	
	objv<-c(objv,obj)
	iter<-iter+1
	}
	
	
	
	
  result$objective <- objv
  result$h <- h
  result$lambda <-lambda
  result$iter <- iter
	return(result)
}

multidembed_uc <- function(A,d,maxiter=20,h0=NULL) {
  m<-length(A)
  n<-dim(A[[1]])[1]
  result<-list("objective" =rep(0,d),"lambda" =matrix(0,m,d),"h" =matrix(0,n,d),"iter"=rep(0,d))
  
  
  ###Innitialize obj, h, lambda, gradient
  Abar<-matrix(0,n,n)
  for(i in 1:m){
    Abar<-Abar+A[[i]]
  }
  
  s<-svd(Abar,d,d)
  h<-s$u
  if (! is.null(h0)) {h <- h0}
  
  
  XX <- matrix(0,d,d)
  hh <- list()
  for(i in 1:d){
    hh[[i]] <- h[,i]%*% t(h[,i])
  }
  
  for(i in 1:d){
    for(j in 1:d){
      XX[i,j] <-  sum(hh[[i]] * hh[[j]])
    }
  }
  
  XY <- matrix(0,d,1)
  lambda <- matrix(0,m,d)
  for(j in 1:m){
    for(i in 1:d){
      XY[i,1] <- sum(hh[[i]] * A[[j]])
    }
    lambda[j,] <- solve(XX) %*% XY
  }
  
  obj<-cptobj_all(A,lambda,h)
  objstart<-obj
  objpre<-Inf
  
  iter<-0
  objv<-c(obj)
  
  while(iter<maxiter & ((objpre-obj)>(10^-7)*(1+objstart))){
    Gradiant<-matrix(0,n,d)
    for(i in 1:m){
      Gradiant<-Gradiant-A[[i]]%*%h%*%diag(lambda[i,])*4+diag(diag(h%*%t(h)))%*%h%*%diag(lambda[i,]^2)*4
    }
    Gradiant<-Gradiant/m
    #Gnorm <- apply(Gradiant, 2, function(x){sqrt(sum(x^2))})
    #angle <- diag((t(Gradiant) %*% h)/Gnorm)
    #Gradiant<-Gradiant %*% diag(1 - angle^2 )
    Gnorm<-norm(Gradiant,type='F')
    print(Gnorm)
    
    
    objtmp<-Inf
    stepsize<-1
    while(objtmp>obj-stepsize*2*Gnorm^2*0.01 & stepsize>10^-7 ){
      ht<-h-Gradiant*stepsize
      #htnorm<-apply(ht, 2, function(x){sqrt(sum(x^2))})
      #ht<-ht%*% diag(1/htnorm)
      objtmp<-cptobj_all(A,lambda,ht)
      stepsize<-stepsize*0.5
    }
    
    h<-ht
    for(i in 1:d){
      hh[[i]] <- h[,i]%*% t(h[,i])
    }
    
    for(i in 1:d){
      for(j in 1:d){
        XX[i,j] <-  sum(hh[[i]] * hh[[j]])
      }
    }
    
    lambda <- matrix(0,m,d)
    for(j in 1:m){
      for(i in 1:d){
        XY[i,1] <- sum(hh[[i]] * A[[j]])
      }
      lambda[j,] <- solve(XX) %*% XY
    }
    
    objpre <- obj
    obj <- cptobj_all(A,lambda,h)
    #print(obj)
    
    objv<-c(objv,obj)
    iter<-iter+1
  }
  
  
  
  
  result$objective <- objv
  result$h <- h
  result$lambda <-lambda
  result$iter <- iter
  return(result)
}



