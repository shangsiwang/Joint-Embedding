#n<-20
#m<-10
#A<-list()
#for (i in 1:m){
#Ai<-matrix(runif(n^2,0,1),n,n)
#A[[i]]<-(Ai+t(Ai))/2
#}
#result1<-onedembed(A,innitialize=1)
#result2<-onedembed(A,innitialize=2)
#result3<-onedembed(A,innitialize=3)
#result4<-multidembed(A,d=2,innitialize=1)
                     


# return (Ai- \sum lambdaprevi hprevi hprevi^T) %*% h
resiMultiply<-function(Ai,h,hprevi=NULL,lambdaprevi=NULL){
if(is.null(hprevi)) {
	return(Ai%*%h)
} else {
      n<-dim(hprevi)[1]
	d<-dim(hprevi)[2]
	Aih<-matrix(0,n,1)
	for(i in 1:d){
		Aih<-Aih + lambdaprevi[i]* hprevi[,i] * ( t(hprevi[,i]) %*% h)
      }
	Aih<-Ai%*%h-Aih
	return(Aih)
}
}






onedembed <- function(A,maxiter=20,hprev=NULL,lambdaprev=NULL,innitialize=2) {
	##Innitialize
	m<-length(A)
	n<-dim(A[[1]])[1]
	result<-list("objective" = 0,"lambda" =0,"h" =0,"iter"=0)

	if(innitialize==1){
	Abar<-matrix(0,n,n)
	# You should use residual matrix.
	for(i in 1:m){
		Abar<-Abar+A[[i]]
	}
	s<-svd(Abar,1,1)
	h<-s$u
	}	else if(innitialize==2) {
	s<-svd(A[[1]],1,1)
	h<-s$u
	} else if(innitialize==3) {
	  h<-matrix(rnorm(n),n,1)
	}
	
	h<-h/norm(h,type="F")
	lambda<-rep(0,m)
	for(i in 1:m){
		lambda[i]<-t(h)%*% resiMultiply(A[[i]],h,hprevi=hprev,lambdaprevi=lambdaprev[i,])
	}

	obj<-cptobjprev(A,lambda,h,hprev,lambdaprev)
	objstart<-obj
	objpre<-Inf


	iter<-0
	while(iter<maxiter & ((objpre-obj)>(10^-7)*(1+objstart))){
		iter<-iter+1
		##Iterate update
		Gradiant<-matrix(0,n,1)
		for(i in 1:m){
			Gradiant<-Gradiant-resiMultiply(A[[i]],h,hprevi=hprev,lambdaprevi=lambdaprev[i,])*lambda[i]*4+h*lambda[i]^2*4
		}
		Gradiant<-Gradiant/m
		Gnorm<-norm(Gradiant,type="F")

		obj<-cptobjprev(A,lambda,h,hprev,lambdaprev)
		stepsize<-1
		objtmp<-Inf
		while(objtmp>obj-stepsize*Gnorm^2*0.01 & stepsize>10^-7 ){
			ht<-h-Gradiant*stepsize
			htnorm<-norm(ht,type="F")
			ht<-ht/htnorm
			objtmp<-cptobjprev(A,lambda,ht,hprev,lambdaprev)
			stepsize<-stepsize*0.5
		}

		if(stepsize<=10^-7 ) {
			objpre<-obj
		} else {
			objpre<-obj
			obj<-objtmp
			h<-ht
			for(i in 1:m){
				lambda[i]<-t(h)%*% resiMultiply(A[[i]],h,hprevi=hprev,lambdaprevi=lambdaprev[i,])

			}
		}
	}
	obj<-cptobjprev(A,lambda,h,hprev,lambdaprev)
	result$objective<-obj
	result$lambda<-lambda
	result$h<-h
	result$iter<-iter
	return(result)
}





multidembed <- function(A,d,maxiter=20,innitialize=2) {
	m<-length(A)
	n<-dim(A[[1]])[1]
	result<-list("objective" =rep(0,d),"lambda" =matrix(0,m,d),"h" =matrix(0,n,d),"iter"=rep(0,d))
	Resid<-A
	for(k in 1:d){
	  if(k==1){
		resultd<-onedembed(A,maxiter=20,hprev=NULL,lambdaprev=NULL,innitialize=innitialize)
	  } else {
	  resultd<-onedembed(A,maxiter=20,as.matrix(result$h[,1:(k-1)]),as.matrix(result$lambda[,1:(k-1)]),innitialize)
	  }
		result$objective[k]=resultd$objective
		result$iter[k]=resultd$iter
		result$lambda[,k]=resultd$lambda
		result$h[,k]=resultd$h
		for(i in 1:m){
			Resid[[i]]<-Resid[[i]]-resultd$lambda[i]*resultd$h%*%t(resultd$h)
		}
	}
	return(result)
}



cptobjprev <- function(A,lambda,h,hprev=NULL,lambdaprev=NULL) {
  m<-length(A)
  d<-0
  if(!is.null(hprev)){
  d<-dim(hprev)[2]
  }
  obj<-0
  for(i in 1:m){
    obj<-obj+sum(A[[i]]*A[[i]]) 
    obj<-obj-2*lambda[i]*t(h)%*% A[[i]]%*% h
    obj<-obj+lambda[i]^2
    if(d>0) {
    for (j in 1:d){
      obj<-obj-2*lambdaprev[i,j]*t(hprev[,j])%*% A[[i]]%*% hprev[,j]
    }
    for(j in 1:d){
      for(k in 1:d){
        obj<-obj+lambdaprev[i,j]*lambdaprev[i,k]*(t(hprev[,j])%*% hprev[,k])^2
      }
      obj<-obj+lambdaprev[i,j]*lambda[i]*(t(hprev[,j])%*% h)^2
    }
    }
    }
  
  return(obj)
}

