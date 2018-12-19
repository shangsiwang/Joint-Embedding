# Simulate data ###################
h1 <- rep(1,200) / sqrt(200)
h2 <- rep(c(1,-1),100)/sqrt(200)
h3 <- rep(c(1,1,-1,-1), 50) / sqrt(200)

m <- 2^10
lambda1 <- runif(n = m, 8,16)
lambda2 <- runif(n = m,  0,2)
lambda3 <- runif(n = m, 0,1)

H1 <- tcrossprod(h1)
H2 <- tcrossprod(h2)
H3 <- tcrossprod(h3)

l1 = lambda1[1]
l2 = lambda2[1]
l3 = lambda3[1]
generate_graph <- function(l1,l2,l3) {
  P <- l1*H1 + l2*H2 + l3*H3
  P <- ifelse(P>1, 1, ifelse(P<0, 0, P))
  P.updi <- P[upper.tri(P)]
  A.updi <- 1*(runif(length(P.updi))<P.updi)
  A <- matrix(0, 200, 200)
  A[upper.tri(A)] <- A.updi
  A <- A + t(A)
  A
}

A.list <- mapply(generate_graph, lambda1, lambda2, lambda3,   SIMPLIFY = FALSE)

#### joint embedding
source("joint_embedding.R")
je.A1 <- multidembed(A = A.list, d = 3, Innitialize = 1, maxiter = 100, large.and.sparse = F)

levelplot(tcrossprod(je.A$h[,1]))
levelplot(tcrossprod(je.A$h[,2]))
levelplot(tcrossprod(je.A$h[,3]))

m <- 2^10
lambda1.t <- runif(n = m, 8,16)
lambda2.t <- runif(n = m,  0,2)
lambda3.t <- runif(n = m, 0,1)
A.list.test <- mapply(generate_graph, lambda1.t, lambda2.t, lambda3.t,   SIMPLIFY = FALSE)
L.test <- multiembed_test(A.list.test, je.A$h)
