n<-100
m<-10
A<-list()
for (i in 1:m){
Ai<-matrix(runif(n^2,0,1),n,n)
A[[i]]<-(Ai+t(Ai))/2
}
result<-multidembed(A,2)


onedembed <- function(A,maxiter=20) {
	##Innitialize
	m<-length(A)
	n<-dim(A[[1]])[1]
	result<-list("objective" = 0,"lambda" =0,"h" =0,"iter"=0)
	Abar<-matrix(0,n,n)
	for(i in 1:m){
		Abar<-Abar+A[[i]]
	}
	s<-svd(Abar,1,1)
	h<-sqrt(s$d[1])*s$u
	h<-h/norm(h,type="F")
	H<-h%*%t(h)
	Hnorm<-norm(H,type="F")^2
	lambda<-rep(0,m)
	for(i in 1:m){
		lambda[i]<-1/(Hnorm)*sum(sum(H*A[[i]]))
	}

	obj<-cptobj(A,lambda,h)
	objstart<-obj
	objpre<-Inf


	iter<-0
	while(iter<maxiter & ((objpre-obj)>(10^-7)*(1+objstart))){
		iter<-iter+1
		##Iterate update
		Gradiant<-matrix(0,n,1)
		for(i in 1:m){
			Gradiant<-Gradiant-A[[i]]%*%h*lambda[i]*4+h*lambda[i]^2*4
		}
		Gradiant<-Gradiant/m
		Gnorm<-norm(Gradiant)

		obj<-cptobj(A,lambda,h)
		stepsize<-1
		objtmp<-Inf
		while(objtmp>obj-stepsize*Gnorm^2*0.01 & stepsize>10^-7 ){
			ht<-h-Gradiant*stepsize
			htnorm<-norm(ht,type="F")
			ht<-ht/htnorm
			Ht<-ht%*%t(ht)
			objtmp<-cptobj(A,lambda,ht)
			stepsize<-stepsize*0.5
		}

		if(stepsize<=10^-7 ) {
			objpre<-obj
		} else {
			objpre<-obj
			obj<-objtmp
			h<-ht
			H<-Ht
			Hnorm<-norm(H,type="F")^2
			for(i in 1:m){
				lambda[i]<-1/(Hnorm)*sum(sum(H*A[[i]]))
			}
		}
	}
	obj<-cptobj(A,lambda,h)
	result$objective<-obj
	result$lambda<-lambda
	result$h<-h
	result$iter<-iter
	return(result)
}




multidembed <- function(A,d,maxiter=20) {
	m<-length(A)
	n<-dim(A[[1]])[1]
	result<-list("objective" =rep(0,d),"lambda" =matrix(0,m,d),"h" =matrix(0,n,d),"iter"=rep(0,d))
	Resid<-A
	for(k in 1:d){
		resultd<-onedembed(Resid)
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




cptobj <- function(A,lambda,h) {
	m<-length(A)
	obj<-0
	for(i in 1:m){
		obj<-obj+norm(A[[i]]-lambda[i]*h%*%t(h),type="F")^2
	}
	return(obj)
}
