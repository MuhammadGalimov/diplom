mcmc <- function(data, iter) {

  # получаем численность населения и количество удаленных индивидуумов
  snotinf <- sum(data$removal != Inf)
  N <- nrow(data)
  
  dd <- 26

  # получаем вектор времен удаления
  removedtimes <- data$removal[1:snotinf]

  # задаем начальные значения времен инфицирования
  #infectedtimes <- seq(0, min(removedtimes), len=snotinf)
  infectedtimes <- rep(NaN, snotinf)
  for (i in 1:snotinf) {
    infectedtimes[i] = removedtimes - dd
  }

  # дополняем таблицу временами заражения
  data$infection <- c(infectedtimes, rep(Inf, N-snotinf))
  data <- data[,c(1,3,2)]

  # начальные значения параметров гамма распределения
  lambda_beta <- 1.0
  nu_beta <- 10^(-3)
  lambda_gamma <- 1.0  
  nu_gamma <- 10^(-3)

  double_sum_cur <- integr(data, snotinf);
  log_cur <- logpr(data, infectedtimes, removedtimes, snotinf);
  sum_minus_cur <- sum(removedtimes-infectedtimes)

  # начальные значения для бетта и гамма
  #beta_cur <- rgamma(1, snotinf - 1 + lambda_beta, 1.0)/(nu_beta + double_sum_cur/N)
  #gamma_cur <- rgamma(1, snotinf + lambda_gamma, 1.0)/(nu_gamma + sum_minus_cur);
  beta_cur <- 0.386
  gamma_cur <- 0.384
  
  res <- matrix(NA, nrow=iter, ncol=3)
  res[1,] <- c(beta_cur, gamma_cur, sum(infectedtimes))
  

  for (i in 2:iter) {

    str(i)

    # выбираем индивидуума, время заражения которого будет обновлено
    choosed_ind <- sample(1:snotinf, 1)

    # получаем предложенное значение
    inftime_cur <- infectedtimes[choosed_ind]
    inftime_can <- removedtimes[choosed_ind] - rexp(1, gamma_cur);

    data$infection[choosed_ind] <- inftime_can;
    infectedtimes[choosed_ind] <- inftime_can;
    
    log_can <- logpr(data, infectedtimes, removedtimes, snotinf);

    if (log_can != -Inf) { 
      
      double_sum_can <- integr(data, snotinf);
      sum_minus_can <- sum(removedtimes-infectedtimes);

      log_q_ratio <- log(dexp(removedtimes[choosed_ind] - inftime_cur, gamma_cur)) - log(dexp(removedtimes[choosed_ind] - inftime_can, gamma_cur))

      log_pi_can <- log_can - (beta_cur/N)*double_sum_can - gamma_cur*sum_minus_can;
      log_pi_cur <- log_cur - (beta_cur/N)*double_sum_cur - gamma_cur*sum_minus_cur;
      
      
      # принимаем новое значение с вероятностью альфа
      u <- runif(1);
      
      if (log(u) < log_pi_can - log_pi_cur + log_q_ratio) {
        log_cur <- log_can;
        double_sum_cur <- double_sum_can;
        sum_minus_cur <- sum_minus_can;
      }
      else {
        data$infection[choosed_ind] <- inftime_cur;
        infectedtimes[choosed_ind] <- inftime_cur;
      }      
    }
    else {
      data$infection[choosed_ind] <- inftime_cur;
      infectedtimes[choosed_ind] <- inftime_cur;
    }
    
    # обновляем бетта
    beta_cur <- rgamma(1, snotinf - 1 + lambda_beta, 1.0)/(nu_beta + double_sum_cur/N)
    
    # обновляем гамма  
    gamma_cur <- rgamma(1, snotinf + lambda_gamma, 1.0)/(nu_gamma + sum_minus_cur)

    # сохраняем значения
    res[i,] <- c(beta_cur, gamma_cur, sum(infectedtimes))
    
  }
  
  res
  
}
