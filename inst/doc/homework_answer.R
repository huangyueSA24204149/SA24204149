## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(boot)
library(bootstrap)
library(coda)
library(doParallel)
library(dplyr)
library(foreach)
library(ggplot2)
library(lpSolve)
library(magrittr)
library(microbenchmark)
library(stats)
library(datasets)
library(parallel)
library(SA24204149)

## -----------------------------------------------------------------------------
#Example 1
x<-c(1,2,3,4,5,6,7,8,9)
result<-median(x)#计算向量x的中位数
result

## -----------------------------------------------------------------------------
#Example 2
matrix1<-matrix(c(1,2,3,4),ncol=2,nrow=2)
matrix2<-matrix(c(5,6,7,8),ncol=2,nrow=2)
output<-matrix1%*%matrix2#计算两个矩阵的乘积
output

## -----------------------------------------------------------------------------
#Example 3
data(ToothGrowth)#导入R自带的数据集（VC与牙齿成长的关系）
model <- lm(len ~ supp * dose, data = ToothGrowth)#建立线性回归模型
summary(model)#R方为0.7296拟合效果良好

## -----------------------------------------------------------------------------
#3.4
#首先由密度函数求出分布函数F,再由F计算出F^(-1)=sqrt(-2*σ^2*ln(1-U))
set.seed(3)
n<-1000
U<-runif(1000)
a=1
b=2
c=3
#σ=1
X1<-sqrt(-2*a^2*log(1-U))
hist(X1,prob=TRUE)
y1<-seq(0,6, .01)
lines(y1,(y1/a^2)*exp(-y1^2/(2*a^2)))
#σ=2
X2<-sqrt(-2*b^2*log(1-U))
hist(X2,prob=TRUE)
y2<-seq(0,8, .01)
lines(y2,(y2/b^2)*exp(-y2^2/(2*b^2)))
#σ=3
X3<-sqrt(-2*c^2*log(1-U))
hist(X3,prob=TRUE)
y3<-seq(0,12, .01)
lines(y3,(y3/c^2)*exp(-y3^2/(2*c^2)))
#由三个图我们可以看出这个算法生成的随机数满足原始分布

## -----------------------------------------------------------------------------
#3.11
# 生成随机样本
set.seed(33)
n<-1000
p1 <- 0.75
p2 <- 1 - p1
n1 <- rnorm(n*p1, mean = 0, sd = 1)
n2 <- rnorm(n*p2, mean = 3, sd = 1)
mixture1<- c(n1, n2)
hist(mixture1,prob=TRUE)
y1<-seq(-4,6, .01)
lines(y1,p1*(1/sqrt(2*pi)*exp(-y1^2/2))+p2*(1/sqrt(2*pi)*exp(-(y1-3)^2/2)))
p3 <- 0.5
p4 <- 1 - p3
n3 <- rnorm(n*p3, mean = 0, sd = 1)
n4 <- rnorm(n*p4, mean = 3, sd = 1)
mixture2<- c(n3, n4)
hist(mixture2,prob=TRUE)
y2<-seq(-4,6, .01)
lines(y2,p3*(1/sqrt(2*pi)*exp(-y2^2/2))+p4*(1/sqrt(2*pi)*exp(-(y2-3)^2/2)))
#由图可以看出密度曲线呈双峰分布，我们可以看到两峰的均值为两个分布的均值，而且哪个p大些对应的峰就高些，当p1=0.5时两峰一样高

## -----------------------------------------------------------------------------
#3.20
set.seed(2)
lambda <-3  # 泊松过程的平均发生率
shape<-2
scale<-1
n <- 1000  # 模拟次数
t <- 10      # 时间长度
simulation<-function(lambda,shape,scale,t){
  N<-rpois(1,lambda*t)
  Y<-rgamma(N,shape=shape,scale=scale)
  return(sum(Y))
}
results<-c()
for(i in 1:n){
  results[i]<-simulation(lambda,shape,scale,t)
}
mean<-mean(results)
var<-var(results)
mean#估计值
var#估计值
mean1 <- lambda * t * (shape/ scale)
var1 <- lambda * t * (shape *(shape+1)*scale^2)
mean1#理论值
var1#理论值

## -----------------------------------------------------------------------------
##5.4
set.seed(2)

monte_carlo_beta <- function(x, n = 100000) {
  samples <- rbeta(n, 3,3)
  c <- sum(samples <= x)
  return(c/n)
}
estimated <- sapply(seq(0.1, 0.9, by = 0.1), monte_carlo_beta)##返回结果向量
true<- pbeta(seq(0.1, 0.9, by = 0.1), 3, 3)
# 输出估计值和实际值进行比较
data.frame(Estimated = estimated, Actual = true)

## -----------------------------------------------------------------------------
##5.9
set.seed(1)
s <- 1
n <- 1000
antithetic <- function(s, n) {
  a <- runif(n)
  X <- s * sqrt(-2 * log(a))
  X_ <- s * sqrt(-2 * log(1 - a))
  return((X + X_) / 2)
}
Sample_antithetic <- antithetic(s, n)
var_antithetic <- var(Sample_antithetic)
rayleigh <- function(s, n) {
  U <- runif(n)
  return(s * sqrt(-2 * log(U)))
}
# 独立采样
X_1 <- rayleigh(s, n)
X_2 <- rayleigh(s, n)
Sample <- (X_1 + X_2) / 2
var_independent <- var(Sample)
reduced <- (var_independent - var_antithetic) / var_independent * 100
cat("独立采样的样本方差为", var_independent, "\n")
cat("对偶变量法采样的样本方差为", var_antithetic, "\n")
cat("方差减少的百分比为：", reduced, "%\n")

## -----------------------------------------------------------------------------
##5.13
g <- function(x) {
  ifelse(x > 1, x ^ 2 / sqrt(2 * pi) * exp(-x ^ 2 / 2), 0)
}
f1 <- function(x) {
  ifelse(x > 1, exp(-(x - 1)), 0)
}
f2 <- function(x) {
  ifelse(x > 1, 2 * exp(-2 * (x - 1)), 0)
}


# 生成数据
set.seed(1)
x <- seq(1.01, 5, by = 0.1)
data <- data.frame(
  x = x,
  g = g(x),
  f1 = f1(x),
  f2 = f2(x)
)

# 绘制图像
ggplot(data, aes(x = x)) +
  geom_line(aes(y = g, color = "g(x)")) +
  geom_line(aes(y = f1, color = "f1")) +
  geom_line(aes(y = f2, color = "f2")) +
  labs(x = "x", y = "函数值", color = "函数") +
  scale_color_manual(values = c(
    "g(x)" = "green",
    "f1" = "red",
    "f2" = "blue"
  ))+
theme_set(theme(plot.margin = margin(0.2, 0.2, 0.2, 02, "cm")))
m <- 1000
samples1 <- rexp(m, 1) + 1
gf1 <- g(samples1) / f1(samples1)
var1 <- var(gf1)
samples2 <- rexp(m, 2) + 1
gf2 <- g(samples2) / f1(samples2)
var2 <- var(gf2)
cat("f1产生的方差为：", var1, "\n")
cat("f2产生的方差为：", var2, "\n")
##f2产生的方差更小

## -----------------------------------------------------------------------------
##question
##定义排序函数
sort <- function(x) {
  if (length(x) <= 1)
    return(x)
  a <- x[sample(length(x), 1)]
  left <- x[x < a]
  right <- x[x > a]
  return(c(sort(left), a, sort(right)))
}
set.seed(1)
values <- c(1 * 10 ^ 4, 2 * 10 ^ 4, 4 * 10 ^ 4, 6 * 10 ^ 4, 8 * 10 ^ 4)
m <- 100 
times <- numeric(length(values))
for (i in seq_along(values)) {
  n <- values[i]
  total <- 0
  for (j in 1:m) {
    x <- 1:n
    start <- Sys.time()
    sort(x)
    end <- Sys.time()
    total <- total + as.numeric(difftime(end, start, units = "secs"))
  }
  times[i] <- total / m
}
t <- values * log(values)
time <- lm(times ~ values)
summary(time)

## -----------------------------------------------------------------------------
set.seed(12345)
skewness <- function(x) {
  n <- length(x)
  mean_x <- mean(x)
  std_dev <- sd(x)
  m3 <- sum((x - mean_x) ^ 3) / n
  return(m3 / std_dev ^ 3)
}
n <- 1000
m <- 1e4

skewness_values <- numeric(m)
for (i in 1:m) {
  x <- rnorm(n)
  skewness_values[i] <- skewness(x)
}

# 计算分位数和标准误差
quantiles <- c(0.025, 0.05, 0.95, 0.975)
estimated_quantiles <- quantile(skewness_values, quantiles)

density <- dnorm(estimated_quantiles)
var <- quantiles * (1 - quantiles) / (n * density ^ 2)
se <- sqrt(var)
large_sample_quantiles <- qnorm(quantiles, mean = 0, sd = sqrt(6 / n))
print(quantiles)#分位数
print(se)#标准差
print(estimated_quantiles)#估计分位数
print(large_sample_quantiles)#大样本分位数

## -----------------------------------------------------------------------------
bivariate_normal <- function(n, m) {
  p_Pearson <- rep(NA, m)
  p_spearman <- rep(NA, m)
  p_kendall <- rep(NA, m)
  for (i in 1:m) {
    x <- rnorm(n)
    y <- rnorm(n)
    cor <- cor.test(x, y)
    p_Pearson[i] <- cor$p.value
    spearman<- cor.test(x, y, method = 'spearman')
    p_spearman[i] <- spearman$p.value
    kendall <- cor.test(x, y, method = 'kendall')
    p_kendall[i] <- kendall$p.value
  }
  list(
    p_Pearson = p_Pearson,
    p_spearman = p_spearman,
    p_kendall = p_kendall
  )
}
##功效函数
compute_power <- function(p_values) {
  mean(p_values < 0.05)
}

power_bivariate_normal <- function(n, m) {
 results <- bivariate_normal(n, m)
  power_Pearson <- compute_power(results$p_Pearson)
  power_spearman <- compute_power(results$p_spearman)
  power_kendall <- compute_power(results$p_kendall)
  cat('Pearson相关检验的功效：', power_Pearson, '\n')
  cat('Spearman检验的功效：', power_spearman, '\n')
  cat('Kendall检验的功效：', power_kendall, '\n')
}

# 生成指数分布的数据并进行检验
exponential <- function(n, m) {
  p_Pearson <- rep(NA, m)
  p_spearman <- rep(NA, m)
  p_kendall <- rep(NA, m)
  for (i in 1:m) {
    x <- rexp(n)
    y <- rexp(n)
    cor<- cor.test(x, y)
    p_Pearson[i] <- cor$p.value
    spearman<- cor.test(x, y, method = 'spearman')
    p_spearman[i] <- spearman$p.value
    kendall<- cor.test(x, y, method = 'kendall')
    p_kendall[i] <- kendall$p.value
  }
  list(
    p_Pearson = p_Pearson,
    p_spearman = p_spearman,
    p_kendall = p_kendall
  )
}
power_exponential <- function(n, m) {
  results <- exponential(n, m)
  power_Pearson <- compute_power(results$p_Pearson)
  power_spearman <- compute_power(results$p_spearman)
  power_kendall <- compute_power(results$p_kendall)
  cat('Pearson相关检验的功效：', power_Pearson, '\n')
  cat('Spearman检验的功效：', power_spearman, '\n')
  cat('Kendall检验的功效：', power_kendall, '\n')
}
set.seed(1009)
n <- 100
m <- 1000
power_bivariate_normal(n, m)# 二元正态分布
cat('\n')
power_exponential(n, m)# 指数分布

## -----------------------------------------------------------------------------
#  Z 统计量
p1 <- 0.651
p2 <- 0.676
n <- 10000
s1 <- sqrt(p1 * (1 - p1) / n)
s2 <- sqrt(p2 * (1 - p2) / n)
s <- sqrt(s1 ^ 2 + s2 ^ 2)
z <- (p2 - p1) / s
#  P 值
p<- 2 * pnorm(-abs(z))
cat("p值为：",p)
if (p < 0.05) {
  cat("\n在 0.05显著水平下，两种方法功效有显著差异。\n")
} else {
  cat("在 0.05显著水平下，两种方法功效没有显著差异。\n")
}

## -----------------------------------------------------------------------------
# 设置参数
N <- 1000
null_hypotheses <- 950
alternative_hypotheses <- 50
alpha <- 0.1
m <- 10000

# 初始化计数器
bonferroni_fwer <- 0
bonferroni_fdr <- 0
bonferroni_tpr <- 0
bh_fwer <- 0
bh_fdr <- 0
bh_tpr <- 0

# 进行m次模拟
for (i in 1:m) {
  # 生成p值
  null_p_values <- runif(null_hypotheses, min = 0, max = 1)
  alternative_p_values <- rbeta(alternative_hypotheses, shape1 = 0.1, shape2 = 1)
  all_p_values <- c(null_p_values, alternative_p_values)
  
  # 计算原始显著性水平下的拒绝情况
  reject_null <- null_p_values <= alpha
  reject_alternative <- alternative_p_values <= alpha
  
  # Bonferroni校正
  bonferroni_threshold <- alpha / N
  reject_null_bonferroni <- null_p_values <= bonferroni_threshold
  reject_alternative_bonferroni <- alternative_p_values <= bonferroni_threshold
  
  # B-H校正
  bh_threshold <- p.adjust(all_p_values, method = "BH", n = N)
  N_k<-seq(1,1000,1)*alpha/N
  reject_null_bh<-bh_threshold<=N_k
  reject_alternative_bh <- alternative_p_values <= N_k
  bonferroni_fwer <- bonferroni_fwer + sum(reject_null_bonferroni) / null_hypotheses
  bonferroni_fdr <- bonferroni_fdr + sum(reject_null_bonferroni) / length(reject_null_bonferroni)
  bonferroni_tpr <- bonferroni_tpr + sum(reject_alternative_bonferroni) / alternative_hypotheses
  bh_fwer <- bh_fwer + sum(reject_null_bh) / null_hypotheses
  bh_fdr <- bh_fdr + sum(reject_null_bh) / length(reject_null_bh)
  bh_tpr <- bh_tpr + sum(reject_alternative_bh) / alternative_hypotheses
}
bonferroni_fwer <- bonferroni_fwer / m
bonferroni_fdr <- bonferroni_fdr / m
bonferroni_tpr <- bonferroni_tpr / m
bh_fwer <- bh_fwer / m
bh_fdr <- bh_fdr / m
bh_tpr <- bh_tpr / m
result_table <- matrix(c(bonferroni_fwer, bonferroni_fdr, bonferroni_tpr,bh_fwer, bh_fdr, bh_tpr), nrow = 3, byrow = TRUE)
colnames(result_table) <- c("Bonferroni ", "B-H ")
rownames(result_table) <- c("FWER", "FDR", "TPR")
print(result_table)

## -----------------------------------------------------------------------------
times <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
# 危险率lambda的极大似然估计值
n <- length(times)
lambda_hat <- n/sum(times)
cat("MLE of the hazard rate lambda：", lambda_hat, "\n")
set.seed(1)
# Bootstrap
B <- 1000 # 设置自助抽样次数
lambda_boot <- numeric(B)
for (i in 1:B) {
  boot_sample <- sample(times, replace = TRUE)
  lambda_boot[i] <- length(boot_sample)/sum(boot_sample)
}

bias <- mean(lambda_boot) - lambda_hat
se <- sd(lambda_boot)

cat("Bias:", bias, "\n")
cat("Standard error:", se)

## -----------------------------------------------------------------------------
set.seed(1)
# 计算平均故障间隔时间
mean_time <- 1/lambda_hat

# 计算 95%置信区间
## 标准正态方法
z <- qnorm(0.975)
normal_interval <- mean_time + c(-z * se, z * se) * (1/lambda_hat^2)

## 基本方法
basic_interval <- quantile(lambda_boot, c(0.025, 0.975))
basic_interval <- 1/basic_interval

## 百分位数方法
percentile_interval <- quantile(1/lambda_boot, c(0.025, 0.975))

## BCa方法
boot_obj <- boot(data = times, statistic = function(x, i) {
  n <- length(x[i])
  return(n/sum(x[i]))
}, R = B)
bca_interval <- boot.ci(boot_obj, type = "bca")$bca[4:5]
bca_interval <- 1/bca_interval

cat("标准正态法 95%置信区间：[", normal_interval[1], ",", normal_interval[2], "]\n")
cat("基本方法 95%置信区间：[", basic_interval[2], ",", basic_interval[1], "]\n")
cat("分位数法 95%置信区间：[", percentile_interval[1], ",", percentile_interval[2], "]\n")
cat("BCa法 95%置信区间：[", bca_interval[2], ",", bca_interval[1], "]\n")

## -----------------------------------------------------------------------------
##7.8
set.seed(1)
sc<-scor
# 计算theta的函数
estimate_theta <- function(x) {
    cov(x) %>%
    prcomp() %>%
    {.$sdev[1]/sum(.$sdev)}
}
n<-nrow(sc)
theta.j<- numeric(n)
for (i in 1:n){    
  theta.j[i]<-estimate_theta(sc[-i,])
}
theta.hat<-estimate_theta(sc)
bias<-(n-1)*(mean(theta.j)-theta.hat) #BIAS
se<-sqrt((n-1)*var(theta.j)) #SE
cat("Bias:", round(bias,3), "\n")
cat("standard error:", round(se,3),"\n")

## -----------------------------------------------------------------------------
#8.2
# 置换检验函数
permutation_test_spearman <- function(x, y, B = 1000) {
  observed_stat <- cor(x, y, method = "spearman")
  n <- length(x)
  stats <- numeric(B)
  for (i in 1:B) {
    permuted_y <- sample(y)
    stats[i] <- cor(x, permuted_y, method = "spearman")
  }
  p_value <- mean(abs(stats) >= abs(observed_stat))
  return(list(observed_statistic = observed_stat, p_value = p_value))
}

# 生成随机数据
set.seed(1)
x <- rnorm(1e3)
y <- rnorm(1e3)

perm_result <- permutation_test_spearman(x, y)
cor_test_result <- cor.test(x, y, method = "spearman")

# 结果打印
round(perm_result$observed_statistic, 3) %>%
  cat("Observed Spearman correlation statistic from permutation test:",
      .,
      "\n")
round(perm_result$p_value, 3) %>%
  cat("P-value from permutation test:", ., "\n")
round(cor_test_result$p.value, 3) %>%
  cat("P-value from cor.test:", ., "\n")

## -----------------------------------------------------------------------------
##9.3
set.seed(123)
# 设置参数
theta <- 1
a <- 0
n_samples <- 10000
delete<- 1000

# Metropolis-Hastings采样器函数
metropolis_hastings <- function(n_samples, theta, a, delete) {
  samples <- numeric(n_samples)
  current_sample <- 0
  for (i in 1:n_samples) {
    proposal <- rnorm(1, mean = current_sample, sd = 1)
    acceptance_ratio <- dcauchy(proposal, location = a, scale = theta) / dcauchy(current_sample, location = a, scale = theta)
    if (runif(1) < acceptance_ratio) {
      current_sample <- proposal
    }
    samples[i] <- current_sample
  }
  return(samples[-(1:delete)])
}

# 生成样本
generated_samples <- metropolis_hastings(n_samples, theta, a, delete)

# 计算并比较分位数
true_quantiles <- qcauchy(seq(0.1, 0.9, by = 0.1), location = a, scale = theta)
generated_quantiles <- quantile(generated_samples, probs = seq(0.1, 0.9, by = 0.1))

# 输出结果
print("True Cauchy Quantiles:")
print(true_quantiles)
print("Generated Cauchy Quantiles:")
print(generated_quantiles)


## -----------------------------------------------------------------------------
##9.8
# 设置参数
n <- 10
a <- 2
b <- 3
num_samples <- 10000

# 初始化x和y
x <- sample(0:n, 1)
y <- runif(1)

# 容器用于存储样本
samples_x <- numeric(num_samples)
samples_y <- numeric(num_samples)

# Gibbs采样循环
for (i in 1:num_samples) {
  # 从条件分布中采样x给定y
  x <- rbinom(1, size = n, prob = dbinom(0:n, size = n, prob = y))
  samples_x[i] <- x
  
  # 从条件分布中采样y给定x
  y <- rbeta(1, shape1 = x + a, shape2 = n - x + b)
  samples_y[i] <- y
}
print(mean(samples_x))
print(mean(samples_y))##输出平均值

## -----------------------------------------------------------------------------
##For each of the above exercise, use the Gelman-Rubin method to monitor convergence of the chain, and run the chain until it converges approximately to the target distribution according to Rˆ < 1.2


# 假设我们已经有一个MCMC样本数据，这里用一个模拟的例子
set.seed(123)
num_chains <- 5
chain_length <- 1000
mcmc_chains <- list()
for (i in 1:num_chains) {
  mcmc_chains[[i]] <- rnorm(chain_length, mean = 0, sd = 1)
}

# 将链转换为矩阵形式
mcmc_matrix <- do.call(cbind, mcmc_chains)

# 定义一个函数来计算Gelman-Rubin统计量
gelman_rubin <- function(chains) {
  n <- nrow(chains)
  m <- ncol(chains)
  B <- var(apply(chains, 2, mean))
  W <- apply(chains, 2, var)
  mean_W <- mean(W)
  var_plus <- (m - 1) / m * W + (n / m) * B
  R_hat <- sqrt(var_plus / W)
  return(R_hat)
}

# 初始计算Gelman-Rubin统计量
R_hat <- gelman_rubin(mcmc_matrix)
cat("Initial Gelman-Rubin statistic:", R_hat, "\n")

# 循环采样直到收敛
while (any(R_hat > 1.2)) {
  # 增加新的样本点到每个链中
  for (i in 1:num_chains) {
    new_samples <- rnorm(chain_length, mean = 0, sd = 1)
    mcmc_chains[[i]] <- c(mcmc_chains[[i]], new_samples)
  }
  
  # 更新矩阵
  mcmc_matrix <- do.call(cbind, mcmc_chains)
  
  # 重新计算Gelman-Rubin统计量
  R_hat <- gelman_rubin(mcmc_matrix)
  cat("Updated Gelman-Rubin statistic:", R_hat, "\n")
}

cat("少于 1.2\n")


## -----------------------------------------------------------------------------
##11.3

compute_kth_term <- function(k, a) {
  d <- length(a)
  norm_a <- sqrt(sum(a^2))  # Euclidean norm
  term <- (-1)^k / (factorial(k) * 2^k) * (norm_a^(2*k+2)) / ((2*k+1)*(2*k+2)) *
    gamma((d+1)/2) * gamma(k + 3/2) / gamma(k + d/2 + 1)
  return(term)
}

# function to compute sum
compute_sum <- function(a, tol = 1e-10) {
  norm_a <- sqrt(sum(a^2))  # Euclidean norm
  sum_series <- 0
  k <- 0
  term <- compute_kth_term(k, a)
  while (abs(term) > tol) {
    sum_series <- sum_series + term
    k <- k + 1
    term <- compute_kth_term(k, a)
  }
  return(sum_series)
}

# Evaluate the sum when a=(1,2)
a <- c(1, 2)
result <- compute_sum(a)
cat('sum when a = (1,2) :', result)

## -----------------------------------------------------------------------------
##11.5
define_equation <- function(k, a){
  int_func <- function(u, n){(1 + u^2/(n - 1))^(-n/2)}
  get_c <- function(n, a){sqrt(a^2 * n / (n + 1 - a^2))}
  expr <- function(n, a) {
    this_int_func <- function(u){
      int_func(u, n)}
    c_val <- get_c(n - 1, a)
    2/sqrt(pi*(n - 1)) * exp(lgamma(n/2)-lgamma((n - 1)/2)) * 
      integrate(this_int_func, lower = 0, upper = c_val)$value}
  
  lhs <- expr(k, a)
  rhs <- expr(k + 1, a)
  return(lhs - rhs)
}
solve_eq <- function(k) {
  f <- function(a){return(define_equation(k, a))}
  eps <- 1e-2
    if (f(eps) < 0 && f(sqrt(k) - eps) > 0 || f(eps) > 0 && f(sqrt(k) - eps) < 0) {
    return(uniroot(f, interval = c(eps, sqrt(k) - eps))$root)
  } else {
    return(NA)
  }
}
r11.5 <- sapply(c(4:25, 100, 500, 1000), function(k) {
  solve_eq(k) %>% 
    return()
})
r11.5

## -----------------------------------------------------------------------------
#  11.4
findIntersection = function (k) {
  s_k.minus.one = function (a) {1 - pt(sqrt(a ^ 2 * (k - 1) / (k - a ^ 2)), df = k - 1)}
  s_k = function (a) {1 - pt(sqrt(a ^ 2 * k / (k + 1 - a ^ 2)), df = k)}
  f = function (a) {s_k(a) - s_k.minus.one(a)}
  
  eps = 1e-2
  return(uniroot(f, interval = c(eps, sqrt(k) - eps))$root)
}

r11.4 <- sapply(c(4:25, 100, 500, 1000), function (k) {
  findIntersection(k)
})

r11.4
##发现在k较小时二者都有非常相近的估计，但11.4相对更稳健。k很大时，EX11.5解的估计相对较差

## -----------------------------------------------------------------------------
##

# 定义观测数据
observed_data <- c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
tau <- 1
n <- length(observed_data)

# 初始化参数
lambda_init <- 1

# E-step: 计算E[Zi | Yi, λ]
e_step <- function(lambda, data, tau) {
  return(pexp(data, rate = lambda))
}

# M-step: 更新λ
m_step <- function(lambda, data, tau) {
  weights <- e_step(lambda, data, tau)
  weighted_sum <- sum(weights * data)
  total_weight <- sum(weights)
  new_lambda <- total_weight / weighted_sum
  return(new_lambda)
}

# E-M算法迭代
em_algorithm <- function(lambda_init, data, tau, tol = 1e-6, max_iter = 100) {
  lambda <- lambda_init
  for (i in 1:max_iter) {
    new_lambda <- m_step(lambda, data, tau)
    if (abs(new_lambda - lambda) < tol) {
      break
    }
    lambda <- new_lambda
  }
  return(lambda)
}

# 运行E-M算法
estimated_lambda <- em_algorithm(lambda_init, observed_data, tau)
cat("Estimated λ using E-M algorithm:", estimated_lambda, "\n")

# 计算MLE
mle_lambda <- function(data, tau) {
  n <- length(data)
  return(n / sum(data))
}

# 计算MLE
mle_estimated_lambda <- mle_lambda(observed_data, tau)
cat("Estimated λ using MLE:", mle_estimated_lambda, "\n")

## -----------------------------------------------------------------------------
##11.7
# 定义目标函数的系数(lpSolve默认是求最大化所以取个负号)
objective <- c(-4, -2, -9)
# 引入s1，s2两个松弛变量，定义约束矩阵
c <- matrix(c(2, 1, 1, 1, 
                        1, -1, 3, 1), 
nrow = 2, byrow = TRUE)
# 定义约束的方向（"<="表示小于等于）
d <- c("<=", "<=")
# 定义约束的右侧值
r<- c(2, 3)
# 用lp函数求解
result <- lp("min", objective, c, d, r)
# 输出结果
print(result$solution)
print(result$objval)
##最终结果再取个负号，所以结果为2

## -----------------------------------------------------------------------------
##EX3

# 定义公式列表
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
# 用 for 循环拟合线性模型
f<- list()
for (i in seq_along(formulas)) {
  f[[i]] <- lm(formulas[[i]], data = mtcars)
}
print(f)
#使用 lapply() 拟合线性模型
l <- lapply(formulas, function(frm) lm(frm, data = mtcars))
print(l)

## -----------------------------------------------------------------------------
##EX4
# 定义公式
formula <- mpg ~ disp
# 创建引导样本列表
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})
# 用 for 循环来拟合线性模型
f<- list()
for (i in seq_along(bootstraps)) {
  f[[i]] <- lm(formula, data = bootstraps[[i]])
}
print(f)
# 用 lapply() 来拟合线性模型
l <- lapply(bootstraps, function(data) lm(formula, data = data))
print(l)

## -----------------------------------------------------------------------------
##EX5
##对前两问得到的结果取个summary$r.squared，再用sapply对每个结果提取即可

## -----------------------------------------------------------------------------
##EX3
# 模拟 t-test 的性能
trials <- replicate(
  5,
  t.test(rpois(9, 9), rpois(6, 9)),
  simplify = FALSE
)
# 匿名函数提取 p-value
p_values1 <- sapply(trials, function(trial) trial$p.value)
print(p_values1)
#  [[ 提取 p-value
p_values2 <- sapply(trials, "[[", "p.value")
print(p_values2)

## -----------------------------------------------------------------------------

# 定义一个并行版本的 lapply 函数
parallel_lapply <- function(X, FUN, ...) {
  # 获取可用的核心数
  num_cores <- detectCores() - 7
  
  # 使用 mclapply 进行并行计算
  result <- mclapply(X, FUN, ..., mc.cores = num_cores)
  
  return(result)
}

# 示例用法
input_vector <- 1:10
square_function <- function(x) {
  Sys.sleep(1)  # 模拟耗时操作
  return(x^2)
}

# 调用并行版本的 lapply
output_vector <- parallel_lapply(input_vector, square_function)
print(output_vector)

## -----------------------------------------------------------------------------
##EX4
fast <- function(x, y) {
  table_xy <- table(x, y)
  observed <- as.vector(table_xy)
  row_sums <- colSums(table_xy)
  col_sums <- rowSums(table_xy)
  total_sum <- sum(observed)
  expected <- outer(row_sums, col_sums, "*") / total_sum
  chisq <- sum((observed - expected)^2 / expected)
  return(chisq)
}

x <- c(1, 2, 3, 4, 5)
y <- c(1, 2, 3, 4, 5)
result <- fast(x, y)
print(result)

## -----------------------------------------------------------------------------
##ex5
fast <- function(x, y) {
  # 获取唯一值及其索引
  x_levels <- unique(x)
  y_levels <- unique(y)
  
  # 创建结果矩阵
  result <- matrix(0, nrow = length(x_levels), ncol = length(y_levels))
  rownames(result) <- as.character(x_levels)
  colnames(result) <- as.character(y_levels)
  
  # 填充结果矩阵
  for (i in seq_along(x)) {
    result[as.character(x[i]), as.character(y[i])] <- result[as.character(x[i]), as.character(y[i])] + 1
  }
  return(result)
}
# 示例用法
x <- c(1, 2, 3, 1, 2, 3, 1, 2, 3)
y <- c(3, 2, 1, 3, 2, 1, 3, 2, 1)
# 使用自定义的fast_table函数
result <- fast(x, y)
print(result)


## -----------------------------------------------------------------------------
# Number of samples to generate
n_samples <- 1000

# Parameters for the Gibbs sampler
n <- 10
a <- 2
b <- 3

# Call the C++ function to generate the samples
result <- Gibbs(n_samples, n, a, b)

# Extract the samples
x_samples <- result$x
y_samples <- result$y
head(x_samples)
head(y_samples)

set.seed(123)

## 使用R自带函数生成随机数
Gibbs_R <- function(n, a, b, n_samples) {
  x_samples_r <- numeric(n_samples)
  y_samples_r <- numeric(n_samples)

  # 初始化y
  y_temp <- runif(1)

  for (i in 1:n_samples) {
    x_samples_r[i] <- rbinom(1, n, y_temp)
    y_temp <- rbeta(1, x_samples_r[i] + a, n - x_samples_r[i] + b)
    y_samples_r[i] <- y_temp
  }

  return(list(x_samples = x_samples_r, y_samples = y_samples_r))
}

## 调用Gibbs_R生成随机数
result_R <- Gibbs_R(n, a, b, n_samples)
x_r_samples <- result_R$x_samples
y_r_samples <- result_R$y_samples
x_r_samples
y_r_samples
## 比较
compare<-microbenchmark(
    Rcpp_function = Gibbs(n, a, b, n_samples),
    R_function = Gibbs_R(n, a, b, n_samples)
) 
  summary(compare)
  ##可以看出Cppfunction时间更短

## -----------------------------------------------------------------------------
set.seed(1)
d1 <- c(-2.961, 0.478, -0.391, -0.869, -0.460, 
            -0.937, 0.779, -1.409, 0.027, -1.569);
d2  <- c(1.608, 1.009,  0.878,  1.600, -0.263,  
             0.680, 2.280,  2.390, 1.793,  8.091, 1.468)
bootstrap<- function(sample1, sample2, R) {
  n1 <- length(sample1)
  n2 <- length(sample2)
  mean_diff <- numeric(R)
  for (i in 1:R) {
    re1 <- sample(sample1, n1, replace = TRUE)
    re2 <- sample(sample2, n2, replace = TRUE)
    mean_diff[i] <- mean(re1) - mean(re2)
  }
  return(sd(mean_diff))
}
result <- bootstrap(d1, d2, R = 10000)
print(result)

