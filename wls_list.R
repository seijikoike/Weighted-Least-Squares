wls.list <- function(linear.mod, iterations = 1){

X <- model.matrix(linear.mod)
H <- solve(t(X) %*% X) %*% t(X)
Y.obs <- linear.mod$model[,1]
size <- nrow(linear.mod$model)
betas.ols <- H %*%  Y.obs #The estimated coefficients using ordinary least squares
resid <- residuals(linear.mod)
W.hat <- diag(size) * as.vector(abs(resid)^2)

beta.wls <- solve(t(X) %*% solve(W.hat) %*% X) %*% t(X) %*% solve(W.hat) %*%  Y.obs
#The estimated coefficients with weight least squares method (1st iteration)

beta.list <- list(betas.ols, beta.wls)

i <- 3
while(i <= iterations){
  
  Y.hat <- X %*% beta.wls
  
  resid <- Y.obs - Y.hat
  
  W.hat <- diag(size) * as.vector(abs(resid)^2)
  
  if(abs(min(eigen(W.hat)$values)) < 1e-15){break} 
  #Stops when the estimated variance/covariance matrix (W.hat) becomes computationally
  #singular
  
  beta.wls <- solve(t(X) %*% solve(W.hat) %*% X) %*% t(X) %*% solve(W.hat) %*% Y.obs  
  
  beta.list[[i]] <- beta.wls
  i <- i+1
}

return(beta.list)
}
