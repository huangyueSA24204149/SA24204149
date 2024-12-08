//' @useDynLib SA24204149
//' @import Rcpp

#include <Rcpp.h>
using namespace Rcpp;
//' @title Loss Function
//' @description calculate loss funtion
//' @param ytrue true vector
//' @param ypred predicted vector
//' @param alpha weight 
//' @return loss
//' @export
// [[Rcpp::export]]
double LossFunction(NumericVector ytrue, NumericVector ypred, double alpha) {
  // 获取向量长度
  int n = ytrue.size();
  
  // 检查输入向量的长度是否相同
  if (n != ypred.size()) {
    stop("The length of the true and predicted vectors must be the same.");
  }
  
  // 初始化损失值
  double loss = 0;
  
  // 计算每个元素的平方误差，并应用权重函数
  for (int i = 0; i < n; i++) {
    double error = ytrue[i] - ypred[i];
    double weight = exp(-alpha * abs(error)); // 使用指数衰减的权重函数
    loss += weight * error * error;
  }
  
  return loss;
}