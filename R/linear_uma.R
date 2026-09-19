# =====================================================================
# TEMPORAL IMPLICIT STRUCTURES (TIS)
# Multiple Linear UMA
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
common_sequencer <- seq(0, 1, length.out = n)
x_true <- 5 + 2 * exp(-0.5 * common_sequencer)
y_true <- 7 + 10 * sqrt(common_sequencer)

data <- data.frame(x_true,y_true)

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
unity_model_analysis <- function(data){
  n <- dim(data)[1]
  m <- dim(data)[2]
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
    metrics[7,i] <- angle_C(metrics[1,i],metrics[3,i],metrics[3,i])
    metrics[8,i] <- lift_h(metrics[1,i],metrics[3,i],metrics[3,i])
  }

  metrics[1,(m+1)] <- sum(metrics[1,1:m])
  metrics[2,(m+1)] <- sum(metrics[2,1:m])
  metrics[3,(m+1)] <- sum(metrics[3,1:m])
  metrics[4,(m+1)] <- cor(unity,predict(model))
  metrics[5,(m+1)] <- (metrics[4,(m+1)])^2
  metrics[6,(m+1)] <- R_sq_modified(metrics[1,(m+1)],metrics[2,(m+1)],metrics[3,(m+1)])
  metrics[7,(m+1)] <- angle_C(metrics[1,(m+1)],metrics[2,(m+1)],metrics[3,(m+1)])
  metrics[8,(m+1)] <- lift_h(metrics[1,(m+1)],metrics[3,(m+1)],metrics[3,(m+1)])

  print(metrics)
}

unity_model_analysis(data)

