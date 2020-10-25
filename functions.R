# считаем количество инфицированных до времени t.

count <- function(infectedtimes, removedtimes, t) {
  sum(infectedtimes < t) - sum(removedtimes < t)
}

# вычисляем интеграл \int_{t} S_t I_t dt используя двойную сумму.

integr <- function(data, snotinf) {
  N <- nrow(data);  
  total <- 0;
  
  for (i in 1:snotinf) {
    total <- total + (N - snotinf) * data$removal[i]
  }
  
  #a <- foreach(i = 1:snotinf, .combine = "+") %do%
  #  (N - snotinf) * data$removal[i]
  
  for (i in 1:snotinf) {
    total <- total + (i + 1 - N) * data$infection[i]
  }
  
  #b <- foreach(i = 1:snotinf, .combine = "+") %do%
  #  (i + 1 - N) * data$infection[i]
  
  for (k in 1:snotinf-2) {
    for (j in k+2:snotinf) {
      total <- total + min(data$removal[k], data$infection[j])
    }
  }
  
  #c <- foreach(k = 1:snotinf - 2, .combine = "+") %:%
  #  foreach(j = k+2:snotinf, .combine = "+") %do%
  #  min(data$removal[k], data$infection[j])
  
  #a + b + c
  total
}  

# вычисляем логарифм произведения \prod_{j!=1} I_{i_j}

logpr <- function(data, infectedtimes, removedtimes, snotinf) {

  time_init <- min(data$infection[1:snotinf])
  index_init <- which(data$infection[1:snotinf]==time_init)

  log_prod <- rep(NA, snotinf);
  
  for (i in 1:snotinf) {
    if (i!=index_init) {
      log_prod[i] <- log(count(infectedtimes, removedtimes, infectedtimes[i]))
    }
    else {
      log_prod[i] <- 0
    }
  }
  sum(log_prod)
}
