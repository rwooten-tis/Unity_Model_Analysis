univariate_analysis <- function(data){
  n <- dim(data)[1]
  m <- dim(data)[2]
  
  # Logically grouped, descriptive row names
  M <- c("Standard Mean", "Median", "Geometric Mean", "Self-Weighted Mean",
         "Minimum", "First Quartile", "Third Quartile", "Maximum",
         "Range", "Variance", "Standard Deviation", "Coeff of Variation",
         "Consistence", "Skewness Index")
  
  N <- names(data)
  metrics <- matrix(nrow = 14, ncol = m, dimnames = list(M, N))
  
  for(i in 1:m){
    x <- data[,i]
    
    # 1. Location
    metrics[1,i]  <- mean(x)
    metrics[2,i]  <- median(x)
    metrics[3,i]  <- exp(mean(log(x))) # Geometric Mean
    metrics[4,i]  <- sum(x^2) / sum(x) # Self-Weighted Mean
    
    # 2. Distribution
    metrics[5,i]  <- min(x)
    metrics[6,i]  <- quantile(x, 0.25)
    metrics[7,i]  <- quantile(x, 0.75)
    metrics[8,i]  <- max(x)
    
    # 3. Spread
    metrics[9,i]  <- max(x) - min(x)
    metrics[10,i] <- var(x)
    metrics[11,i] <- sd(x)
    metrics[12,i] <- sd(x) / mean(x)   # CV
    
    # 4. Shape & UMA Invariants
    metrics[13,i] <- n * (mean(x)^2) / sum(x^2) # Consistence
    metrics[14,i] <- 3 * (mean(x) - median(x)) / sd(x) # Pearson Median Skewness
  }
  
  print(metrics)
}

