# =====================================================================
# TEMPORAL IMPLICIT STRUCTURES (TIS)
# Multiple Linear Regression Unity Model Analysis
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
x2 <- 100 + 2 * exp(time) + rnorm(n,0,0.1)
x1 <- 120 + time  + rnorm(n, 0, 1)

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

mlr_uma <- function(data){
  n <- dim(data)[1]
  m <- dim(data)[2]
  for(i in 1:m){
    if(!is.numeric(data[,i])){
      print("ERROR: Data frame must contain numeric data only")
    }
  }
  unity <- rep(1,n)
  N <- names(data)
  df <- as.matrix(data)
  model <- lm(unity ~ -1 + df)
  summary(model)
  alpha <- coef(model)
  df_scaled <- df
  for(i in 1:m){
    df_scaled[,i] <- alpha[i] * df[,i]
  }
  df_est <- df
  df_error <- df
  df_model <- df
  df_total <- df
  metrics <- matrix(nrow = 8, ncol = (m+1), dimnames = list(c("SSE","SSM","SST","R","R_sq","R_sq_modified","Theta","Lift_h"),c(N,"Overall")))
  metrics_slr <- matrix(nrow = 5, ncol = (m), dimnames = list(c("SSE","SSM","SST","R","R_sq"),N))
  
  for(i in 1:m){
    for(j in 1:n){
      df_est[j,i] <- (1 - sum(df_scaled[j,-i])) / alpha[i]
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
    model_slr <- lm(df[,i] ~ df[,-i])
    metrics_slr[1,i] <- sum((df[,i] - predict(model_slr))^2)
    metrics_slr[2,i] <- sum((predict(model_slr) - mean(df[,i]))^2)
    metrics_slr[3,i] <- sum((df[,i] - mean(df[,i]))^2)
    metrics_slr[4,i] <- cor(df[,i],predict(model_slr))
    metrics_slr[5,i] <- summary(model_slr)$r.square
    
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
  print(metrics_slr)
}

mlr_uma(data)

