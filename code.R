source("functions.R")
source("mcmc.R")

# загрузка данных
data <- read.csv("covid_data.csv")

# численность населения
N <- 4000000

# количество выздоровевших и умерших
n <- data[5,92] + data[3,92]


times <- 1:length(data[1,-1])
nofdead <- as.numeric(as.vector(data[5,-1]) + as.vector(data[3,-1]))

removal <- c()
for (i_ in 1:(length(nofdead)-1)) {
  difference <- nofdead[i_+1] - nofdead[i_]
  if (difference != 0) {
    for (j_ in 1:difference) {
      removal <- append(removal, times[i_+1])
    }
  }
}

removal <- append(removal, rep(Inf, N - n))
label <- 1:N

# собираем DataFrame для алгоритма
data <- data.frame(label = label, removal = removal)

res <- mcmc(data, 2000)

plot(1:3000, res[,3])
plot(as.vector(res[,3]), type = "l", xlab = "Iteration", ylab = "sum_of_inf_times")

mean(res[2900:3000,1]) / mean(res[2900:3000,2])

rite.csv(res, file = "result2.csv")
res1000 <- read.csv("result1.csv")
plot(as.vector(res1000[,4]), type = "l", xlab = "Iteration", ylab = "sum_of_inf_times")

sum(removal[1:n]) / n - 50000 / n

rexp(1, 0.384)
