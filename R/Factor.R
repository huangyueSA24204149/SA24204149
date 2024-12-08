#' @importFrom Rcpp sourceCpp
#' @useDynLib SA24204149

NULL
##########################################################################
# Function factor to estimate number of factors
#     rhat0 <- factor(X,rmax)
# X TxN, maxit maximum number of iteration, defalut is four
##########################################################################
#' Estimate number of factors
#' 
#' return the number of factors
#' @param X a given matrix
#' @param rmax the largest number of factor
#' @return the number of factor 
#' @export
factor <- function(X, rmax){
  # 设置最大迭代次数为4
  maxiteration <- 4
  
  # 获取矩阵的列数和行数
  N <- NCOL(X)
  T <- NROW(X)
  
  # 中心化数据矩阵
  xd <- X - matrix(1,T,1)%*%colMeans(X)
  
  # 计算协方差矩阵
  v <- t(xd)%*%xd/T
  
  # 计算特征值和特征向量
  temp <- eigen(v)
  
  # 初始化j为rmax + 1
  j <- rmax + 1
  
  # 提取第j到j+4个特征值
  lamb <- temp$values[j:(j+4)]
  
  # 计算相邻特征值之差
  dl <- temp$values[1:rmax] - temp$values[2:(rmax+1)]
  
  # 生成趋势项
  trend <-  c((j-1):(j+3))^(2/3)
  
  # 创建全1向量
  x1vec <- matrix(1,5,1)
  
  # 构造设计矩阵
  xx <- cbind(trend,x1vec)
  
  # 计算X'X
  XpX <- t(xx)%*%xx
  
  # 计算(X'X)的逆矩阵
  invXpX <- solve(XpX)
  
  # 计算最小二乘估计向量
  bhatvec <- invXpX %*% t(xx) %*% lamb
  
  # 计算delta值
  delta <- 2*abs(bhatvec[1])
  
  # 如果所有dl都小于delta，则rhat0设为0
  if (all(dl<delta)) {rhat0 <- 0
  } else {
    # 否则找到第一个大于delta的dl对应的索引
    rhat0 <- max((which(dl>delta)))
    if (rhat0<rmax){
      iter <- 1
      while (iter < (maxiteration)){
        j <- rhat0 + 1
        lamb <- temp$values[j:(j+4)]
        dl <- temp$values[1:rhat0] - temp$values[2:(rhat0+1)]
        trend <-  c((j-1):(j+3))^2/3
        xx <- cbind(trend,x1vec)
        XpX <- t(xx)%*%xx
        invXpX <- solve(XpX)
        bhatvec <- invXpX %*% t(xx) %*% lamb
        delta <- 2*abs(bhatvec[1])
        rhat <- max((which(dl>delta)))
        if (rhat == rhat0) {iter <- iter + (maxiteration+1)}
        iter <- iter +1
      }	
      rhat0 <- rhat
    }
  }
  return(rhat0)
}