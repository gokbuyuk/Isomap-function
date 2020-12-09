## Isomap function
## Takes a data set, a localization parameter eps, desired dimension d, 
##    returns: coordinates of a representation  in d-dim, 

library(usedist) #dist_make(data matrix, fcn)

isomap = function(X, eps=NULL, k=NULL, d){
  n = nrow(X)
  Delta_E = as.matrix(dist(X, diag = TRUE, upper = TRUE)) #Find pairwise distances 
  if (eps != 'NULL'){
    Delta_E[Delta_E  > eps] = Inf # Replace values greater than eps with infinity
  }
  Delta_G = Delta_E
  ## Find the shortest path between the points i&j in the epsilon graph
  E = matrix(0,nrow=n,ncol=n) 
  m = 1
  while (m < n-1){
    for (i in 1:n){
      for (j in 1:n){
        E[i,j] = min(Delta_G[i,]+Delta_G[,j])
      }
    }
    Delta_G = E
    m = 2*m
  }
  R0 = cmdscale(Delta_G, d)
  for (c in 1:20) {
    #R0 = mds.guttman.eq(R0,Delta_G)
    
    #
    # Guttman transform with equal weights
    #
    D <- as.vector(as.dist(Delta_G))/as.vector(dist(R0))
    
    #Delta_G <- matrix(1:n,nrow=n,ncol=n)
    i <- as.vector(as.dist(Delta_G))
    j <- as.vector(as.dist(t(Delta_G)))
    Y <- matrix(0,nrow=n,ncol=ncol(R0))
    for (k in 1:length(D)){
      v <- R0[i[k],]-R0[j[k],]
      v <- D[k]*v
      Y[i[k],] <- Y[i[k],]+v
      Y[j[k],] <- Y[j[k],]-v
    }
    R0 = Y/n
  }
  R1 = R0[, 1:d]
  return(R1)
}

### Examples ###

# Rectangular annulus
n1 = 100
R <- matrix(0,nrow=1,ncol=2)
for (i in 1:n1) {
  x1 <- runif(1,min=0,max=3)
  y1 <- runif(1,min=0,max=1)
  x2 <- runif(1,min=0,max=1)
  y2 <- runif(1,min=1,max=3)
  x3 <- runif(1,min=2,max=3)
  y3 <- runif(1,min=1,max=3)
  x4 <- runif(1,min=0,max=3)
  y4 <- runif(1,min=3,max=4)
  R <- rbind(R,c(x1,y1),c(x2,y2),c(x3,y3),c(x4,y4))
}
plot(R)

R1 = isomap(R, eps=0.5, 2)
plot(R1)

X=R
eps=0.5
Delta_E = as.matrix(dist(X, diag = TRUE, upper = TRUE))
Delta.eps = graph.eps(Delta_E,eps)
Delta_G = as.matrix(graph.short(Delta.eps))

# plot the graph
plot(X)
for (i in 1:401){for (j in 1:401){
  if (Delta_G[i,j] <= eps){
    segments(X[i,1],X[i,2], X[j,1],X[j,2],
             col = "blue", lty = par("lty"), lwd = par("lwd"))
  }}}


# Spiral

n2 = 200
S = data.frame(x = numeric(n2),
               y = numeric(n2),
               z = numeric(n2))
for (i in 1:n2){
  t = runif(1, 0, 1)
  S$x[i] = cos(8*t)
  S$y[i] = sin(8*t)
  S$z[i] = t
}

plot(S)
scatterplot3d(S)

S2 = isomap(S, eps=0.3, d=2)
plot(S2)


##Diagnostics
X=S
eps=0.3
Delta_E = as.matrix(dist(X, diag = TRUE, upper = TRUE))
Delta.eps = graph.eps(Delta_E,eps)
Delta_G = as.matrix(graph.short(Delta.eps))

plot(X$y,X$z)

for (i in 1:n2){for (j in 1:n2){
  if (Delta_G[i,j] <= eps){
    segments(X$y[i], X$z[i], X$y[j], X$z[j],
             col = "blue", lty = par("lty"), lwd = par("lwd"))
  }}}

plot(X$x,X$z)

for (i in 1:n2){for (j in 1:n2){
  if (Delta_G[i,j] <= eps){
    segments(X$x[i], X$z[i], X$x[j], X$z[j],
             col = "blue", lty = par("lty"), lwd = par("lwd"))
  }}}


plot(X$y,X$x)

for (i in 1:n2){for (j in 1:n2){
  if (Delta_G[i,j] <= eps){
    segments(X$y[i], X$x[i], X$y[j], X$x[j],
             col = "blue", lty = par("lty"), lwd = par("lwd"))
  }}}



# On Torus
R=1
r=0.5
n3 = 200
Torus = data.frame(x = numeric(n3),
                   y = numeric(n3),
                   z = numeric(n3))
for (i in 1:n3){
  u = runif(1, 0, pi)
  v = runif(1, 0, pi)
  Torus$x[i] = (R+r*cos(u))*cos(v)
  Torus$y[i] = (R+r*cos(u))*sin(v)
  Torus$z[i] = r*sin(u)
}

scatterplot3d(Torus)
plot(Torus)

T2 = isomap(Torus, eps=0.5, 2)
plot(T2)
