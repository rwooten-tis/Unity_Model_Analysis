# =====================================================================
# TEMPORAL IMPLICIT STRUCTURES (TIS)
# Quadratic Manifold Regression Unity Model Analysis
# =====================================================================

# --- 1. Functional Code

windowsFonts("A" = windowsFont("Times New Roman"))
setwd("E:\\Research\\TIS")
create_file <- function(string){ return(paste(getwd(),"\\",string,".jpg",sep=""))}
publishJPG <- function(string, code){
  file_name <- create_file(string)
  jpeg(file_name,width=4000,height=2500,res=300)
  code
  dev.off()
}

# --- 2. Initialization and Sample Simulation
set.seed(5)

n <- 100
time <- seq(0,1,length.out = n)
x3 <- 100 + 3 * time + rnorm(n,0,0.2)
x2 <- 100 + 2 * time + rnorm(n,0,0.1)
x1 <- 120 + x2 + x2 * x3   + rnorm(n, 0, 1)


data <- data.frame(x1, x2, x3)

angle_C <- function(a2,b2,c2){
  a <- sqrt(a2)
  b <- sqrt(b2)
  C <- acos((a2 + b2 - c2)/ (2*a*b))
  return(C*180/pi)
}

lift_h <- function(a2,b2,c2){
  a <- sqrt(a2)
  B <- angle_C(a2,c2,b2)*pi/180
  return(a*sin(B))
}

R_sq_modified <- function(a2,b2,c2){
  a <- sqrt(a2)
  b <- sqrt(b2)
  C <- angle_C(a2,b2,c2)
  sst <- c2 + 2*a*b*cos(C*pi/180)
  r_sq_modified <- b2/sst
  return(r_sq_modified)
}

qmr_uma <- function(data){
  n <- dim(data)[1]
  m <- dim(data)[2]
  for(i in 1:m){
    if(!is.numeric(data[,i])){
      print("ERROR: Data frame must contain numeric data only")
    }
  }
  k <- m*(m + 1)/2 + m
  unity <- rep(1,n)
  N <- names(data)
  df <- as.matrix(data)
  df_extended <- matrix(nrow = n, ncol = k)
  tracker <- matrix(rep(0,m*k), nrow = m, ncol = k)
  Names <- character(k)
  
  # Extend Data Frame to Include Interaction Terms
  index = 0
  for(i in 1:m){
    index <- index + 1
    df_extended[,index] <- data[,i]
    Names[index] = N[i]
    tracker[i,index] = 1
  }
  
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      index <- index + 1
      df_extended[,index] <- df[,i] * df[,j]
      Names[index] <- paste(N[i],N[j],sep = "")
      tracker[i,index] = 1
      tracker[j,index] = 1
    }
  }
  
  for(i in 1:m){
   index <- index + 1
   df_extended[,index] <- df[,i] * df[,i]
   Names[index] <- paste(N[i],N[i],sep = "")
   tracker[i,index] = 1
  }
  
  DF <- as.data.frame(df_extended)
  names(DF) <- Names
  Data <- as.matrix(DF)
  model <- lm(unity ~ -1 + Data)
  print(summary(model))
  alpha <- coef(model)
  df_est <- matrix(nrow = n, ncol = m)
  
  # Corrected Estimation Block
  for(i in 1:m){
    # Clean element-wise filtering to keep objects as vectors of length k
    alpha_x <- alpha * tracker[i,]
    alpha_c <- alpha * (1 - tracker[i,])
    
    # Compress directly into vectors of length n
    p <- m*(m+1)/2 + i
    df_x2 <- alpha[p]
    df_x <- df_extended[,-p] %*% alpha_x[-p]
    df_c <- df_extended %*% alpha_c
    
    for(j in 1:n){
      A <- df_x2
      B <- df_x[j] / df[j,i]
      C <- (df_c[j] - 1)
      D <- B^2 -4*A*C
      print(length(D[D<0]))
      D[D<0] <- 0
      lower <- (-B - sqrt(D))/(2*A)
      upper <- (-B + sqrt(D))/(2*A)
      # Bypasses matrix swelling; isolates x_i cleanly using its linear factor properties
      df_est[j,i] <- ifelse(abs(df[j,i] - lower) < abs(df[j,i] - upper),lower,upper)
    }
  }
  
  df_error <- df
  df_model <- df
  df_total <- df
  metrics <- matrix(nrow = 8, ncol = (m+1), dimnames = list(c("SSE","SSM","SST","R","R_sq","R_sq_modified","Theta","Lift_h"),c(N,"Overall")))
  
  for(i in 1:m){
    for(j in 1:n){
      df_total[j,i] <- df[j,i] - mean(df[,i])
      df_model[j,i] <- df_est[j,i] - mean(df[,i])
      df_error[j,i] <- df[j,i] -  df_est[j,i]
    }
    metrics[1,i] <- sum((df_error[,i])^2)
    metrics[2,i] <- sum((df_model[,i])^2)
    metrics[3,i] <- sum((df_total[,i])^2)
    metrics[4,i] <- cor(df[,i],df_est[,i])
    metrics[5,i] <- (metrics[4,i])^2
    metrics[6,i] <- R_sq_modified(metrics[1,i],metrics[2,i],metrics[3,i])
    metrics[7,i] <- angle_C(metrics[1,i],metrics[2,i],metrics[3,i])
    metrics[8,i] <- lift_h(metrics[1,i],metrics[2,i],metrics[3,i])
  }
  
  metrics[1,(m+1)] <- sum(metrics[1,1:m])
  metrics[2,(m+1)] <- sum(metrics[2,1:m])
  metrics[3,(m+1)] <- sum(metrics[3,1:m])
  metrics[4,(m+1)] <- sqrt(summary(model)$r.square)
  metrics[5,(m+1)] <- summary(model)$r.square
  metrics[6,(m+1)] <- R_sq_modified(metrics[1,(m+1)],metrics[2,(m+1)],metrics[3,(m+1)])
  metrics[7,(m+1)] <- angle_C(metrics[1,(m+1)],metrics[2,(m+1)],metrics[3,(m+1)])
  metrics[8,(m+1)] <- lift_h(metrics[1,(m+1)],metrics[2,(m+1)],metrics[3,(m+1)])
  
  print(summary(model))
  print(metrics)  
  
}

qmr_uma(data)
