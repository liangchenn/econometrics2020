#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector CDF(NumericVector x1, NumericVector x2, float beta1, float beta2){
  
  int len = x1.size();
  NumericVector res(len);
  
  for(int i = 0; i < len; i++){
    
    float u = x1[i]*beta1 - x2[i]*beta2;
    
    res[i] = exp(u)/(1 + exp(u));
    
  }
  
  return(res);
}

// [[Rcpp::export]]
float C_F(NumericVector x1, NumericVector x2, NumericVector y, float beta1, float beta2){
  
  int len = x1.size();
  //NumericVector res(len);
  float summation = 0;
  
  for(int i = 0; i < len; i++){
    
    float u = x1[i]*beta1 - x2[i]*beta2;
    
    float temp = exp(u)/(1 + exp(u));
    
    
    summation += y[i]*log(temp) + (1-y[i])*log(1-temp);
    
    
    
  }
  
  return(summation);
  
}


// [[Rcpp::export()]]
NumericVector CppLoglikelihood(NumericVector beta1, NumericVector beta2){
  
  Environment env = Environment::global_env();
  NumericVector x1 = env["x1"];
  NumericVector x2 = env["x2"];
  NumericVector y  = env["y"];
  
  int len = beta1.size();
  
  NumericVector res(len);
  
  for(int i = 0; i < len; i++){
    
    
    res[i] = C_F(x1, x2, y, beta1[i], beta2[i]);
    
    
  }
  
  return(res);
  
}








