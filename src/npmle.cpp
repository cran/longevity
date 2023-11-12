// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

//' Turnbull's sets
//'
//' Given truncation and censoring sets,
//' construct disjoint increasing intervals whose
//' left and right endpoints lie in L and R
//' and which contain no other members of L and R
//'
//' @param Lcens set of left censoring limits
//' @param Rcens vector of right censoring limits
//' @param Ltrunc vector of left truncation limits
//' @param Rtrunc vector of right truncation limits
//' @param status integer vector giving status of censoring set
//' @export
//' @return a matrix containing limits of intervals for EM
//' @keywords internal
// [[Rcpp::export(.turnbull_intervals)]]
arma::dmat turnbull_intervals(
    arma::dvec Lcens,
    arma::dvec Rcens,
    arma::dvec Ltrunc,
    arma::dvec Rtrunc,
    arma::uvec status
){
  arma::uword n = Lcens.n_elem;
  if(n != Rcens.n_elem){
    Rcpp::stop("Left and right censoring vectors must be of the same length.");
  }
  // if(Ltrunc.n_elem != Rtrunc.n_elem){
  //   Rcpp::stop("Left and right truncation vectors must be of the same length.");
  // }
  double eps = sqrt(arma::datum::eps);
  // Sort and get unique elements of left and right sets
  // Ltrunc = arma::unique(Ltrunc);
  // Rtrunc = arma::unique(Rtrunc);

  // Make intervals open to the left by adding epsilon
  arma::dvec Lcens_o = Lcens;
  for(arma::uword i = 0; i < n; i++){
   if(status(i) != 1){
     Lcens_o(i) += eps;
   }
  }

  arma::dvec left = arma::sort(arma::unique(arma::join_cols(Lcens_o, Rtrunc)));
  arma::dvec right = arma::sort(arma::unique(arma::join_cols(Rcens, Ltrunc))) + eps/2;
  // Create container, reduce size latter
  arma::uword nmax = std::min(left.n_elem, right.n_elem);
  arma::dmat turnset(nmax, 2);
  arma::uword Lind = 0;
  arma::uword Rind = 0;
  arma::uword i = 0;
  bool keepgoing = true;
  // FIRST INTERVAL ONLY
  // Start with Lind = 0 since min(Rset) > min(Lset)
  while(keepgoing){
  for(arma::uword j = Rind; j < right.n_elem; ++j){
    Rind = j;
    // find the next value of R > L
    // stop at first occurrence
    if(right(j) >= left(Lind)){
      // Increment left set if there was a gap
      for(arma::uword k = Lind; k < left.n_elem; ++k){
        if(left(k) <= right(j)){
          Lind = k;
        } else{
          // Stop as soon as we exceed current upper bound
          break;
        }
      }
      break;
    }
  }
  // Save interval
  turnset(i,0) = left(Lind);
  turnset(i,1) = right(Rind) - eps/2;
  i++;
  // Increment the index set for the left
  Lind += 1;
  // If no more left truncation point or we have already
  // included the rightmost point of Rind, abort
  if((Lind >= left.n_elem) | ((Rind + 1) >= right.n_elem)){
    keepgoing = false;
  }
  } // end of while loop
  return turnset.head_rows(i);
}

//' Identification sets
//'
//' @param tsets Turnbull's sets
//' @param lcens numeric vector of left censoring
//' @param rcens numeric vector of right censoring
//' @param ltrunc numeric vector of left truncation
//' @param rtrunc numeric vector of right truncation
//' @param trunc logical are observation truncated?
//' @export
//' @return a matrix with the bounds of the intervals for Turnbull sets
//' @keywords internal
// [[Rcpp::export(.censTruncLimits)]]
arma::umat censTruncLimits(
  arma::dmat tsets,
  arma::dvec lcens,
  arma::dvec rcens,
  arma::dvec ltrunc,
  arma::dvec rtrunc,
  bool trunc,
  bool cens
){
  arma::uword n = lcens.n_elem;
  if((lcens.n_elem != n) || (rcens.n_elem != n)){
    Rcpp::stop("All vectors of censoring intervals should be of the same length.");
  }
  arma::uword J = tsets.n_rows;
  if(tsets.n_cols != 2){
    Rcpp::stop("\"tsets\" should be a matrix with two columns.");
  }
  double eps = sqrt(arma::datum::eps);
  // Initialize containers
  arma::umat bounds(n, 4, arma::fill::value(J));
  // Columns correspond to censLow, censUpp, truncLow, truncUpp
for(arma::uword i = 0; i < n; ++i){
  for(arma::uword j = 0; j < J; ++j){
   // fully observed
    if(fabs(rcens(i) - lcens(i)) < eps){ // check equality status == 1
      // Compute I_{ij} for censoring - the description in Turnbull 1976 says the contrary...
      // L_i, R_i should cover [q_j, p_j]
      if((tsets(j,0) >= lcens(i)) && (tsets(j,1) <= rcens(i))){
        // Initialized to J, so pick the smallest integer
        if(j < bounds(i,0)){ // initialized to J, so true for the first
          bounds(i,0) = j;
        }
        // Since sets are increasing, pick the largest
        bounds(i,1) = j;
      }
    } else{ // Interval data, consider (left, right)
      if((tsets(j,0) >= (lcens(i) + eps)) && (tsets(j,1) <= rcens(i))){
        // Initialized to J, so pick the smallest integer
        if(j < bounds(i,0)){ // initialized to J, so true for the first
          bounds(i,0) = j;
        }
        // Since sets are increasing, pick the largest
        bounds(i,1) = j;
      }
    }
  }
}
if(!trunc){
  for(arma::uword i = 0; i < n; ++i){
    bounds(i,2) = 0;
    bounds(i,3) = J - 1;
  }
} else{
  if((ltrunc.n_elem != n) | (rtrunc.n_elem != n)){
    Rcpp::stop("All vectors of censoring intervals should be of the same length.");
  }
  for(arma::uword i = 0; i < n; ++i){
    for(arma::uword j = 0; j < J; ++j){
    if((ltrunc(i) <= tsets(j,0)) && (rtrunc(i) >= tsets(j,1))){
        // Initialized to J, so pick the smallest integer
        if(j < bounds(i,2)){
          bounds(i,2) = j;
        }
        // Since sets are increasing, pick the largest
        bounds(i,3) = j;
      }
    }
  }
}
// Rcpp::Rcout << bounds << std::endl;
// Intervals; if the value is J, then intervals are empty
return bounds;

}


//' Turnbull EM algorithm (low storage implementation)
//'
//' @param tsets Turnbull's sets
//' @param lcens numeric vector of left censoring
//' @param rcens numeric vector of right censoring
//' @param ltrunc numeric vector of left truncation
//' @param rtrunc numeric vector of right truncation
//' @param cens logical; if \code{FALSE}, then \code{censUpp = censLow} and a particular update can be avoided in the EM algorithm
//' @param tol tolerance level for terminating the EM algorithm
//' @param maxiter maximum number of iteration for the EM algorithm
//' @param weights vector of weights for observations
//' @return a list with the probabilities and the standard errors
//' @keywords internal
// [[Rcpp::export(.turnbull_em)]]
Rcpp::List turnbullem(
    arma::dmat tsets,
    arma::dvec lcens,
    arma::dvec rcens,
    arma::dvec ltrunc,
    arma::dvec rtrunc,
    arma::dvec weights,
    bool cens = true,
    bool trunc = true,
    double tol = 1e-12,
    double zerotol = 1e-10,
    arma::uword maxiter = 1e5){
  arma::uword n = lcens.n_elem;
  // Create containers
  if(weights.n_elem != n){
    Rcpp::stop("Invalid weight vector");
  }
  if((lcens.n_elem != n) | (rcens.n_elem != n)){
    Rcpp::stop("All vectors of censoring intervals should be of the same length.");
  }
  if(any(lcens > rcens)){
      Rcpp::stop("Invalid arguments: \"lcens\" > \"rcens\" for some elements.");
  }
  arma::uword J = tsets.n_rows;
  if(tsets.n_cols != 2){
    Rcpp::stop("\"tsets\" should be a matrix with two columns.");
  }
  arma::umat limits = censTruncLimits(
    tsets,
    lcens,
    rcens,
    ltrunc,
    rtrunc,
    trunc,
    cens);
  arma::uvec censLow = limits.col(0);
  arma::uvec censUpp = limits.col(1);
  arma::uvec truncLow = limits.col(2);
  arma::uvec truncUpp = limits.col(3);

  // Create containers
  arma::dvec pCur(J);
  arma::dvec grad(J);
  arma::dvec pNew(J, arma::fill::zeros);
  // Equiprobable initial value
  pNew.fill(1/(double)J);
  arma::dvec uiCum(J, arma::fill::zeros);
  double sum_p = 0;
  double abstol = 0;
  int niter = 0;
  int nviolation = 0;
  bool convergence = false;
  bool mle = true;
  double N = arma::sum(weights);
  // Lagrange multiplier of equality constraint
  double mu0 = 0;
  if(cens & !trunc){
    mu0 = N;
  } else if (!cens & trunc){
    mu0 = -N;
  }
  arma::dvec mu(J);
  // Rcpp::Rcout << limits << std::endl;
  while(!convergence && niter < maxiter){
    abstol = arma::max(arma::abs(pCur - pNew));
    // Check convergence of the method
    if((abstol < tol) && (niter > 1)){
      // If difference is negligible, check the KKT conditions
      nviolation = 0;
      for(arma::uword j = 0; j < J; ++j){
        if(pCur(j) > 0){
          mu(j) = 0;
        } else if(pCur(j) == 0){
          mu(j) = mu0 - grad(j);
          if(mu(j) <= 0){
            nviolation += 1;
            // Zero-ed component shouldn't be zero...
            pCur(j) = 1.2*zerotol;
          }
        }
      }
      if(nviolation > 0){
        mle = false;
        convergence = false;
        pCur = pCur/arma::sum(pCur);
      } else{
      convergence = true;
      }
    }

    grad.zeros();
    pCur = pNew;
    Rcpp::checkUserInterrupt();
    uiCum.zeros();
    for(arma::uword i = 0; i < n; ++i){
      // Censoring step - only
      if(censLow(i) < J){
        // This clause is not triggered unless there are censored observations
        if(cens){
          sum_p = arma::sum(pCur(arma::span(censLow(i), censUpp(i))));
          for(arma::uword j = censLow(i); j <= censUpp(i); ++j){
            uiCum(j) += weights(i) * pCur(j)/sum_p;
            grad(j) += weights(i) / sum_p; //KKT conditions discussed in 2.2 of Gentleman and Geyer
          }
        } else {
           uiCum(censLow(i)) +=  weights(i);
           grad(censLow(i)) +=  weights(i); // since alpha_{ij}=1 for all, denominator is 1
        }
      }
      if(trunc & (truncLow(i) < J)){
        // Truncation step
        sum_p = arma::sum(pCur(arma::span(truncLow(i), truncUpp(i))));
        for(arma::uword j = 0; j < J; ++j){
          if(j < truncLow(i) || j > truncUpp(i)){
            uiCum(j) += weights(i) * pCur(j) / sum_p;
          } else {
            grad(j) -= weights(i) / sum_p; //NEW
          }
        }
      }
    }
    pNew = uiCum / arma::sum(uiCum);
    // Rcpp::Rcout << pNew << std::endl;
    // Zero-out probabilities that are negligible
    for(arma::uword j = 0; j < J; ++j){
      if(pNew(j) < zerotol){
        pNew(j) = 0;
      }
    }
    //NEW Check KKT conditions with the gradient after zeroing entries

    pNew = pNew / arma::sum(pNew);
    niter ++;
  }
  return Rcpp::List::create(Rcpp::Named("p") = pNew,
                            Rcpp::Named("rgrad") = mu,
                            Rcpp::Named("conv") = convergence,
                            Rcpp::Named("abstol") = abstol,
                            Rcpp::Named("niter") = niter,
                            Rcpp::Named("KKT") = mle);
}
