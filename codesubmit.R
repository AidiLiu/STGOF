library(iterators)
library(parallel)
library(foreach)
library(doParallel)
library(mixsmsn)


####密度函数
dST <- function(x, loc = 0, sigma2 = 1,shape=1, nu = 4){
  d <- (x - loc)/sqrt(sigma2)
  dens <- 2*dt(d, df = nu)*pt(sqrt((1+nu)/(d^2+nu))*d*shape,1+nu)/sqrt(sigma2)
  return(dens)
}

pST <- function(y, mu, sigma2,lambda, nu){
  dens  <- sn::pst(y, mu, sqrt(sigma2),lambda, nu)
  return(dens)
}

pSN <- function(y, mu, sigma2,lambda){
  dens  <- sn::psn(y, mu, sqrt(sigma2),lambda)
  return(dens)
}

##########检验所用函数和统计量##
AD_value <- function(y, mu, sigma2,lambda, nu){
  z <- sort(y)
  n <- length(y)
  x       <- rep(0,n)
  for(i in 1:n){
    x[i]  <- (2*i-1)*(log(pST(z[i],mu,sigma2,lambda,nu))+log(1-pST(z[n+1-i],mu,sigma2,lambda,nu)))
  }
  return( -n -1/n*sum(x))
}


cdfW <- function(x, nu) {
  result <- pbeta((1 + x / sqrt(x^2 + nu)) / 2, nu / 2, nu / 2)
  result[result <= .Machine$double.xmin] <- .Machine$double.xmin
  result[result >= 1] <- 0.99999  
  return(result)
}




AW_value <- function(y,nu){
  z <- sort(y)
  n <- length(y)
  x       <- rep(0,n)
  for(i in 1:n){
    x[i]  <- (2*i-1)*(log(cdfW(z[i],nu))+log(1-cdfW(z[n+1-i],nu)))
  }
  
  return( -n -1/n*sum(x))
}


cdfH <- function(x, nu) {
  result <- pf(x, 1, nu)
  result[result <= .Machine$double.xmin] <- .Machine$double.xmin
  result[result >= 1] <- 0.99999  
  return(result)
}




AH_value <- function(y,nu){
  z <- sort(y)
  n <- length(y)
  x       <- rep(0,n)
  for(i in 1:n){
    x[i]  <- (2*i-1)*(log(cdfH(z[i],nu))+log(1-cdfH(z[n+1-i],nu)))
  }
  return( -n -1/n*sum(x))
}



cdfRW <- function(x, nu) {
  result <- pbeta((1 + x / sqrt(x^2 + nu)) / 2, nu / 2, nu / 2)
  result[result <= .Machine$double.xmin] <- .Machine$double.xmin
  result[result >= 1] <- 0.99999  
  return(result)
}


cdfRH  <- function(x,nu){
  result <- pf(x, 1, nu)
  result[result <= .Machine$double.xmin] <- .Machine$double.xmin
  result[result >= 1] <- 0.99999  
  return(result)
}


Fnw       <- function(x){
  z <- c()
  for(i in 1:length(x)){
    z[i] <- (sum(x <= x[i])-0.25)/(length(x)+0.25)                                          
  }
  return(z)
}

gen.Skew.normal <- function(n, mu, sigma2, shape, nu=NULL){

  delta <- shape / sqrt(1 + shape^2)
  y <- mu*rep(1,n) + sqrt(sigma2)*(delta*abs(rnorm(n)) + (1 - delta^2)^(1/2)*rnorm(n))
  return(y)
}

gen.Skew.t <- function(n, mu, sigma2, shape, nu ){
  
  y <- mu + (rgamma(n, nu/2, nu/2))^(-1/2)*gen.Skew.normal(n, 0, sigma2, shape)
}


AWHRRvalue  <- function(n,mu,sigma2,lambda,nu){
  y.boost <- gen.Skew.t(n,mu,sigma2,lambda,nu)
  result <- sn::st.mple(y = y.boost, opt.method = "nlminb")
  para1 <- result$dp
  mu.star <- para1[1]
  sigma2.star <- para1[2]^2
  lambda.star <- para1[3]
  nu.star     <- pmax(para1[4],0.1)
  
  W <- abs(y.boost - mu.star) / sqrt(sigma2.star) * sign(rnorm(n))
  AW <- AW_value(W, nu.star)
  
  H <- (y.boost - mu.star)^2 / sigma2.star
  AH <- AH_value(H, nu.star)
  
  v <- abs(y.boost - mu.star)* sign(rnorm(n))
  RW <- v / sqrt(sigma2.star)
  u <- qt(Fnw(RW), nu.star)
  r <- cor(u, v)
  
  vstar <- (y.boost - mu.star)^2
  RH <- vstar / sigma2.star
  ustar <- qf(Fnw(RH), 1, nu.star)
  rstar <- cor(ustar, vstar)
  return(c(AW, AH, r, rstar))
}

nreject1          <- matrix(NA,q,4)
nreject2          <- matrix(NA,q,4)
nreject3          <- matrix(NA,q,4)
nreject4          <- matrix(NA,q,4)
nreject5          <- matrix(NA,q,4)
nreject6          <- matrix(NA,q,4)

AWH_vec           <- matrix(NA,4,B)


  for (i in 1:q) {
      yst    <- gen.Skew.t(n,mu,sigma2,lambda,nu)
      result <- try(sn::st.mple(y = yst, opt.method = "nlminb"), silent = TRUE)
      para   <- result$dp
      mu.hat <- para[1]
      sigma2.hat <- para[2]^2
      lambda.hat <- para[3]
      nu.hat     <- para[4]
      
  for(j in 1:B){
    AWH_vec[,j] <- AWHRRvalue(n,mu.hat, sigma2.hat, lambda.hat, nu.hat)
  }
      AW_vec <- AWH_vec[1,]
      AH_vec <- AWH_vec[2,]
      r_vec <- AWH_vec[3,]
      rstar_vec <- AWH_vec[4,]
      
      criti1 <- quantile(as.numeric(AW_vec), c(0.90,0.95,0.99,0.999))
      criti2 <- quantile(as.numeric(AH_vec), c(0.90,0.95,0.99,0.999))
      criti3 <- quantile(as.numeric(r_vec), c(0.10,0.05,0.01,0.001))
      criti4 <- quantile(as.numeric(rstar_vec),c(0.10,0.05,0.01,0.001))
      
      # 如果 criti1 和 criti2 没有 Inf，则继续存储数据
      nreject1[i,] <- criti1
      nreject2[i,] <- criti2
      nreject3[i,] <- criti3
      nreject4[i,] <- criti4
      print(i)
  }
  
  AW         <- apply(nreject1, 2, function(x) mean(x, na.rm = TRUE))
  AH         <- apply(nreject2, 2, function(x) mean(x, na.rm = TRUE))
  RW         <- apply(nreject3, 2, function(x) mean(x, na.rm = TRUE))
  RH         <- apply(nreject4, 2, function(x) mean(x, na.rm = TRUE))
  
  df <- data.frame(AW = AW, AH = AH, RW = RW, RH = RH)
  rownames(df) <- c("alpha=0.1", "alpha=0.05", "alpha=0.01", "alpha=0.001")
  df