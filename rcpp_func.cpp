#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

//' Identify for which index the cumulative sum of prob is smaller than value
//' 
//' Identify for which index the cumulative sum of probabilities (prob) is smaller than
//' a given value (value). This function helps with multinomial draws
//' 
//' @param value a real number between 0 and 1
//' @param prob a vector of probabilities that sum to one
//' @return this function returns the integer res
//' @export

int whichLessDVPresence(double value, NumericVector prob) {
  int res=-1;
  double probcum = 0;
  
  for (int i = 0; i < prob.length(); i++) {
    probcum = probcum + prob(i); //calculate the cumulative sum of the vector prob
    if (value < probcum) {
      res = i;
      break;
    }
  }
  return res;
}

//' Summarize data according to the different location and species groups
//' 
//' This function determines how many observations are in each location group l 
//' and species group s for two cases: 
//' - when dat(i,j)=1: these results are stored in matrix nql1
//' - when dat(i,j)=0: these results are stored in matrix nql0
//' 
//' @param z integer vector with cluster assignments of locations
//' @param w integer vector with cluster assignments of species
//' @param dat this matrix has L rows (e.g., locations) and S columns (e.g., species)
//'        and contains the presence-absence data
//' @param ngrloc maximum number of location groups
//' @param ngrspp maximum number of species groups
//' @return this function returns a list with two matrices: nql1 and nql0
//' @export

// [[Rcpp::export]]
Rcpp::List getql(IntegerVector z, IntegerVector w, IntegerMatrix dat,
                 int ngrloc, int ngrspp) {
  IntegerMatrix nql1(ngrloc,ngrspp);
  IntegerMatrix nql0(ngrloc,ngrspp);
  
  for(int i=0; i<dat.nrow();i++){
    for (int j=0; j<dat.ncol(); j++){
      if (dat(i,j)==1) {
        nql1(z[i],w[j])=nql1(z[i],w[j])+1;
      }
      if (dat(i,j)==0){
        nql0(z[i],w[j])=nql0(z[i],w[j])+1;
      }
    }
  }
  
  Rcpp::List resTemp = Rcpp::List::create(Rcpp::Named("nql1") = nql1,
                                          Rcpp::Named("nql0") = nql0);
  return(resTemp);
}

//' Convert the vector v into the vector res using the stick-breaking formula
//' 
//' This function converts the input vector v into another vector of probabilities 
//' that sum to one based on the stick-breaking equation
//' 
//' @param v vector of probabilities (these do not have to sum to 1)
//' @return this function returns a vector res with probabilities that sum to one
//' @export
//' 
// [[Rcpp::export]]
NumericVector convertSBtoNormal(NumericVector v) {
  NumericVector res(v.size());
  
  res[0]=v[0];
  double prod=1-v[0];
  for(int j=1; j<v.size();j++){
    res(j)=v[j]*prod;    
    prod=prod*(1-v[j]);
  }
  
  return (res);
}

//' Sample the cluster assignment z for each location 
//' 
//' This function samples the cluster assignment z for each location
//' 
//' @param ltheta this is log(theta)
//' @param dat this matrix has L rows (e.g., locations) and S columns (e.g., species)
//'        and contains the presence-absence data
//' @param dat1m this matrix has L rows (e.g., locations) and S columns (e.g., species)
//'        and is calculated as 1-dat
//' @param lpsi this is log(psi)
//' @param l1mpsi this is log(1-psi)
//' @param w this is the cluster assignment for each species
//' @param runi this is a vector of uniform random variables
//' @return this function returns a vector z with the cluster assignment of each locations
//' @export
//' 
// [[Rcpp::export]]
IntegerVector samplez(NumericVector ltheta,
                      IntegerMatrix dat,
                      IntegerMatrix dat1m,
                      NumericMatrix lpsi,
                      NumericMatrix l1mpsi,
                      IntegerVector w,
                      NumericVector runi) {
  IntegerVector z(dat.nrow());
  NumericVector prob(ltheta.size());
  double tmp=0;
  
  for(int i=0; i<dat.nrow();i++){
    //calculate log-probabilities
    for (int k=0; k<ltheta.size(); k++){
      tmp=0;
      for (int j=0; j<dat.ncol(); j++){ //sum over all species
        tmp=tmp+dat(i,j)*lpsi(k,w[j])+dat1m(i,j)*l1mpsi(k,w[j]);      
      }        
      prob[k]=tmp+ltheta[k]; //add prior log-probability for that group
    }
    prob=prob-max(prob);  //for numerical stability
    prob=exp(prob);       //exponentiate log-probability
    prob=prob/sum(prob); //normalize to sum to 1
    
    //multinomial draw
    z[i]=whichLessDVPresence(runi[i],prob)+1;
  }
  return (z);
}

//' Sample the cluster assignment w for each species
//' 
//' This function samples the cluster assignment w for each species
//' 
//' @param lphi this is log(phi)
//' @param dat this matrix has L rows (e.g., locations) and S columns (e.g., species)
//'        and contains the presence-absence data
//' @param dat1m this matrix has L rows (e.g., locations) and S columns (e.g., species)
//'        and is calculated as 1-dat
//' @param lpsi this is log(psi)
//' @param l1mpsi this is log(1-psi)
//' @param z this is the cluster assignment for each location
//' @param runi this is a vector of uniform random variables
//' @return this function returns a vector w with the cluster assignment of each species
//' @export
// [[Rcpp::export]]
IntegerVector samplew(NumericVector lphi,
                      IntegerMatrix dat,
                      IntegerMatrix dat1m,
                      NumericMatrix lpsi,
                      NumericMatrix l1mpsi,
                      IntegerVector z,
                      NumericVector runi) {
  IntegerVector w(dat.ncol());
  NumericVector prob(lphi.size());
  double tmp=0;
  
  for(int j=0; j<dat.ncol();j++){
    //calculate log-probabilities
    for (int k=0; k<lphi.size(); k++){
      tmp=0;
      for (int i=0; i<dat.nrow(); i++){ //sum over all locations
        tmp=tmp+dat(i,j)*lpsi(z[i],k)+dat1m(i,j)*l1mpsi(z[i],k);      
      }        
      prob[k]=tmp+lphi[k]; //add log prior-probability
    }
    prob=prob-max(prob); //for numerical stability
    prob=exp(prob);      //exponentiate log-probabilities
    prob=prob/sum(prob); //normalize results to sum to 1
    
    //multinomial draw
    w[j]=whichLessDVPresence(runi[j],prob)+1;
  }
  return (w);
}

