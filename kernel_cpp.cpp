#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
#include <stdlib.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;



vec K_sob_cpp2(subview_col<double>  s, double t){
    vec k1s = s-0.5;
    double k1t = t-0.5;
    vec k1abst = abs(s-t)-0.5;

    vec k2s = (k1s % k1s - 1.0/12)/2;
    double k2t = (k1t * k1t - 1.0/12)/2;

    vec k4abst = (pow(k1abst,4) - k1abst%k1abst/2 + 7.0/240)/24;

    return 1.0 + k1s * k1t + k2s * k2t - k4abst;
}



mat getK_sob(const vec & t){
    int n = t.n_elem, i;
    mat out(n,n);
    for (i=0; i<n; i++){
        out.submat(i, i, n-1, i) = K_sob_cpp2(t.subvec(i, n-1), t[i]);
    }
    return symmatl(out);
}

// [[Rcpp::export]]
mat getK_I(const vec & t){
    int n = t.n_elem, i, j;
    mat out(n,n);
    for (i=0; i<n; i++){
      for (j=0; j<i+1; j++){
        if(t(i) == t(j)){out(i,j) = 1;}else{
          out(i,j) = 0;
        }
      }
    }
    return symmatl(out);
}


// [[Rcpp::export]]
mat getK_sob_prod(const mat & X){
  mat K = ones<mat>(X.n_rows, X.n_rows);
  int i;
  for (i=0; i<X.n_cols; i++){
    K = K % getK_sob(X.col(i));
  }
  return K;
}

// [[Rcpp::export]]
mat getK_I_prod(const mat & X){
  mat K = ones<mat>(X.n_rows, X.n_rows);
  int i;
  for (i=0; i<X.n_cols; i++){
    K = K % getK_I(X.col(i));
  }
  return K;
}

// [[Rcpp::export]]
vec K_gaussian_cpp(const vec &s, double t, double h){
  vec temp = square(s-t);
  return exp(-0.5 * temp /(h*h));
}

// [[Rcpp::export]]
mat getG_grid(const vec & t1, const vec & t2, double h){
  // t1 is V;
  // t2 is Vg
    int n1 = t1.n_elem, n2 = t2.n_elem, i;
    mat out(n1,n2);
    for (i=0; i<n1; i++){
        out.row(i) = K_gaussian_cpp(t2, t1(i), h);
    }
    return out;
}



//[[Rcpp::export]]
List eval_obj_grad2_cpp(const vec &w, int N,  const vec &ind, const mat&nP1, const vec & nq1, double lam2, const mat & GG, const vec & Vnw_w){
  vec z = -1.0 * ones(N);
  uvec  tind = find(ind > 0);
  z.elem(tind) = w - 1.0;
  double V = sum(Vnw_w.elem(tind) % w % w)/N;
  mat tZnP1 = (diagmat(z) * nP1).t();
  mat Mtemp = tZnP1 * GG * tZnP1.t();
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, Mtemp + diagmat(nq1));
  uword maxidx = index_max(eigval);
  vec v = eigvec.col(maxidx);
  List out;
  double obj = eigval(maxidx) + lam2 * V;
  vec tzv = tZnP1.t() * v;
  vec gra = 2 *  (tzv.elem(tind)) % (GG.rows(tind) * tzv) + 2.0/N * lam2 * Vnw_w.elem(tind) % w;
  NumericVector gra_r = wrap(gra);
  gra_r.attr("dim") = R_NilValue;
  out["objective"] = obj;
  out["gradient"] = gra_r;
  return out;
}


//[[Rcpp::export]]
vec eval_obj_cpp(const vec &w, int N,  const vec &ind, const mat&nP1, const vec & nq1, double lam2, const mat & GG, const vec & Vnw_w){
  vec z = -1.0 * ones(N);
  uvec  tind = find(ind > 0);
  z.elem(tind) = w - 1.0;
  double V = sum(Vnw_w.elem(tind) % w % w)/N;
  mat tZnP1 = (diagmat(z) * nP1).t();
  mat Mtemp = tZnP1 * GG * tZnP1.t();
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, Mtemp + diagmat(nq1));
  uword maxidx = index_max(eigval);
  vec obj(1);
  obj(0) = eigval(maxidx) + lam2 * V;
  return obj;
}



//[[Rcpp::export]]
vec eval_grad_cpp(const vec &w, int N,  const vec &ind, const mat&nP1, const vec & nq1, double lam2, const mat & GG, const vec & Vnw_w){
  vec z = -1.0 * ones(N);
  uvec  tind = find(ind > 0);
  z.elem(tind) = w - 1.0;
  double V = sum(Vnw_w.elem(tind) % w % w)/N;
  mat tZnP1 = (diagmat(z) * nP1).t();
  mat Mtemp = tZnP1 * GG * tZnP1.t();
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, Mtemp + diagmat(nq1));
  uword maxidx = index_max(eigval);
  vec v = eigvec.col(maxidx);
  vec tzv = tZnP1.t() * v;
  vec gra = 2 *  (tzv.elem(tind)) % (GG.rows(tind) * tzv) + 2.0/N * lam2 * Vnw_w.elem(tind) % w;
  //std::cout << (sign(gra)% gra).max() << (sign(gra)% gra).min() << '\n';
  return gra;
}
