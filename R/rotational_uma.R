library(shiny)
library(scatterplot3d)

# GUI Environment
windowsFonts("A" = windowsFont("Times New Roman"))

# Functional Code
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

rotational_uma <- function(data){
  n <- dim(data)[1]
  m <- dim(data)[2]
  for(i in 1:m){
    if(!is.numeric(data[,i])){
      print("ERROR: Data frame must contain numeric data only")
    }
  }
  k <- m + m*(m - 1)/2 + m
  unity <- rep(1,n)
  N <- names(data)
  df <- as.matrix(data)
  df_extended <- matrix(nrow = n, ncol = k)
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
    }
  }
  for(i in 1:m){
    index <- index + 1
    df_extended[,index] <- df[,i] * df[,i]
    Names[index] <- paste(N[i],N[i],sep = "")
  }

  DF <- as.data.frame(df_extended)
  names(DF) <- Names
  Data <- as.matrix(DF)
  df_est <- matrix(nrow = n, ncol = k)


  df_error <- df_extended
  df_model <- df_extended
  df_total <- df_extended
  metrics <- matrix(nrow = 8, ncol = k, dimnames = list(c("SSE","SSM","SST","R","R_sq","R_sq_modified","Theta","Lift_h"),Names))

  for(i in 1:k){
    df_response <- Data[,i]
    df_explanatory <- Data[,-i]

    model <- lm(df_response ~ -1 + df_explanatory)
    df_est[,i] <- predict(model)

    for(j in 1:n){
      df_total[j,i] <- df_extended[j,i] - mean(df_extended[,i])
      df_model[j,i] <- df_est[j,i] - mean(df_extended[,i])
      df_error[j,i] <- df_extended[j,i] -  df_est[j,i]
    }
    metrics[1,i] <- sum((df_error[,i])^2)
    metrics[2,i] <- sum((df_model[,i])^2)
    metrics[3,i] <- sum((df_total[,i])^2)
    metrics[4,i] <- cor(df_extended[,i],df_est[,i])
    metrics[5,i] <- (metrics[4,i])^2
    metrics[6,i] <- R_sq_modified(metrics[1,i],metrics[2,i],metrics[3,i])
    metrics[7,i] <- angle_C(metrics[1,i],metrics[2,i],metrics[3,i])
    metrics[8,i] <- lift_h(metrics[1,i],metrics[2,i],metrics[3,i])
  }

  print(metrics)

  unity_model <- lm(unity ~ -1 + df_extended)
  alpha <- coef(unity_model)
  manifold_est <- predict(unity_model, data = df_est)

  plot(unity - manifold_est)
  hist(unity - manifold_est)

}

