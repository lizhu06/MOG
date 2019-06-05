#include <RcppEigen.h>
#include <R.h>
#include <Rmath.h>

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace Rcpp;
using namespace std;
using Eigen::Map;                 
using Eigen::MatrixXd;                  
using Eigen::VectorXd;                  

// [[Rcpp::export]]
void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}

// [[Rcpp::export]]
void update_gamma1_sc(double pi1, double sigma2, int G, int P,  
  Eigen::VectorXi& g_index, Eigen::VectorXd& gamma_1, 
  Eigen::VectorXd& gamma_1_rep_per_feature,
  Eigen::VectorXd& gamma_2, Eigen::VectorXd& b, 
  Eigen::VectorXd& Y, Eigen::MatrixXd& X) {
  
  Eigen::VectorXd tailvec = gamma_2.array() * b.array();
  Eigen::ArrayXd beta_old = gamma_1_rep_per_feature.array() * tailvec.array();
  Eigen::VectorXd res_old = Y - X * beta_old.matrix();
  double mul = 1.0 / (2.0 * sigma2);
  double L0_L1;

  for(int g = 0; g < G; g++) {
    int fg = g+1;
    Eigen::VectorXd res_change = res_old;
    if(gamma_1[g] == 0){
      for(int p = 0; p < P; p++){
        if(g_index[p] == fg){
          res_change = res_change - X.col(p) * tailvec[p];
        }
      }
      L0_L1 = exp(mul * (res_change.dot(res_change) - 
        res_old.dot(res_old)));
    }else{
      for(int p = 0; p < P; p++){
        if(g_index[p] == fg){
          res_change = res_change + X.col(p) * tailvec[p];
        }
      }
      L0_L1 = exp(mul * (res_old.dot(res_old) - 
        res_change.dot(res_change)));
    }
    double p_gamma_1 = pi1 / (pi1 + (1.0-pi1) * L0_L1);
    double new_gamma_1 = R::rbinom(1, p_gamma_1);
    if(new_gamma_1 != gamma_1[g]){
      res_old = res_change;
      gamma_1[g] = new_gamma_1;
    }else{
      res_change = res_old;
    }
  }
}


// [[Rcpp::export]]
void update_gamma2_sc(Eigen::VectorXd& pi2_rep_per_feature, double sigma2, int P, 
  Eigen::VectorXi& g_index, Eigen::VectorXi& m_fg,
  Eigen::VectorXd& gamma_1_rep_per_feature, 
  Eigen::VectorXd& gamma_2, 
  Eigen::VectorXd& b, Eigen::VectorXd& beta, Eigen::VectorXd& Y, Eigen::MatrixXd& X, 
  Eigen::VectorXd& v_weight, int BernoulliWeighted_int) {

  Eigen::VectorXd tailvec = gamma_1_rep_per_feature.array() * b.array();
  beta = gamma_2.array() * tailvec.array();
  Eigen::VectorXd res_old = Y - X * beta.matrix();
  double mul = 1.0 / (2.0 * sigma2);
  double L0_L1;

  for(int p = 0; p < P; p++) {
    //if(m_fg[g_index[p] - 1] > 1){
      Eigen::VectorXd res_change = res_old;
      if(gamma_2[p] == 0){
        res_change = res_change - X.col(p) * tailvec[p];
        L0_L1 = exp(mul * (res_change.dot(res_change) - 
          res_old.dot(res_old)));
      }else{
        res_change = res_change + X.col(p) * tailvec[p];
        L0_L1 = exp(mul * (res_old.dot(res_old) - 
          res_change.dot(res_change)));
      }

      double p_gamma_2 = 0.0;
      if(BernoulliWeighted_int == 1){
        double pi2_weighted = pi2_rep_per_feature[p]*v_weight[p];
        p_gamma_2 = pi2_weighted/(pi2_weighted 
          + (1.0-pi2_weighted) * L0_L1);
      }else{
        p_gamma_2 = pi2_rep_per_feature[p] / (pi2_rep_per_feature[p] 
          + (1.0-pi2_rep_per_feature[p]) * L0_L1);
      }
      double new_gamma_2 = R::rbinom(1, p_gamma_2);

      if(new_gamma_2 != gamma_2[p]){
        res_old = res_change;
        gamma_2[p] = new_gamma_2;
        beta[p] = gamma_1_rep_per_feature[p] * gamma_2[p] * b[p];
      }else{
        res_change = res_old;
      }
    //}
  }
}

// [[Rcpp::export]]
void update_b_sc(Eigen::VectorXd& pi2_rep_per_feature, double sigma2, 
  Eigen::VectorXd& s2_vec, int P,
  Eigen::VectorXd& gamma_1_rep_per_feature, 
  Eigen::VectorXd& gamma_2, 
  Eigen::VectorXd& b, Eigen::VectorXd& beta, 
  Eigen::VectorXd& Y, Eigen::MatrixXd& X, Eigen::VectorXd& w_weight) {

  Eigen::VectorXd res = Y - X * beta.matrix();

  for(int p = 0; p < P; p++) {
    if(gamma_1_rep_per_feature[p] * gamma_2[p] == 0){
      b[p] = R::rnorm(0.0, sqrt(s2_vec[p] * w_weight[p]));
    }else{
      double sigma2_b = 1.0 / ((X.col(p).dot(X.col(p))) / sigma2 + 1.0 / (w_weight[p]*s2_vec[p]));
      double mu_b = X.col(p).dot((res + X.col(p) * beta[p])) / sigma2 * sigma2_b;
      b[p] = R::rnorm(mu_b, sqrt(sigma2_b));
    }
    double betap_new = gamma_1_rep_per_feature[p] * 
      gamma_2[p] * b[p];
    res = res - X.col(p) * (betap_new - beta[p]);
    beta[p] = betap_new;
  }
}

// [[Rcpp::export]]
double update_pi1_sc(Eigen::VectorXd& gamma_1, 
  double alpha1, double beta1, int G) {
  double sum_g = gamma_1.sum();
  double pi1 = R::rbeta(sum_g + alpha1, G - sum_g + beta1);
  return (pi1);
}


// [[Rcpp::export]]
void update_pi2_sc(Eigen::VectorXd& gamma_2, Eigen::VectorXd& pi2, 
  Eigen::VectorXi& g_index, 
  Eigen::VectorXi& m_fg, double alpha2, double beta2, int G, int P, 
  Eigen::VectorXd& v_weight,  
  Eigen::VectorXd& pi2_loglikeli, double pi2_prop_n, int MH_ind) {

  for(int g = 0; g < G; g++){
    //if(m_fg[g] > 1){
      double fg = g + 1;
      if(MH_ind == 1){ // MH
        double loglikeli_new=0.0;
        //double lower = max(0.0, pi2[g] - pi2_prop_sigma);
        //double upper = min(pi2[g]+pi2_prop_sigma, 1.0);
        //double pi2_new = R::runif(lower, upper);
        double pi2_new = R::rbeta(pi2_prop_n * pi2[g], 
          pi2_prop_n * (1.0 - pi2[g]));
        

        //if(pi2_new <= 0){
        //  pi2_new = 0.0001;
        //}
        //if(pi2_new >= 1){
        //  pi2_new = 0.9999;
        //}

        for(int p = 0; p < P; p++){
          if(g_index[p] == fg){
            if(gamma_2[p]==1){
              loglikeli_new = loglikeli_new + 
                log(pi2_new * v_weight[p]);
            }else{
              loglikeli_new = loglikeli_new + 
                log(1-(pi2_new * v_weight[p]));
            }
          }
        }

        double acc_p = min(0.0, loglikeli_new - pi2_loglikeli[g] + 
          R::dbeta(pi2[g], pi2_prop_n*pi2_new, pi2_prop_n*(1-pi2_new), 1) -
          R::dbeta(pi2_new, pi2_prop_n*pi2[g], pi2_prop_n*(1-pi2[g]), 1));
        double u_uni = R::runif(0, 1);
        if(log(u_uni) < acc_p){
          pi2[g] = pi2_new;
          pi2_loglikeli[g] = loglikeli_new;
        }

      } else {
        //double fg = g + 1;
        double sum_p = 0.0;
        for(int p = 0; p < P; p++){
          if(g_index[p] == fg){
            sum_p = sum_p + gamma_2[p];
          }
        }
        pi2[g] = R::rbeta(sum_p + alpha2, m_fg[g] - sum_p + beta2);
      }
    //}
  }
  //return (pi2);
}

// [[Rcpp::export]]
double update_sigma2_sc(double N, double alpha_sigma, double beta_sigma,
  Eigen::VectorXd& gamma_1_rep_per_feature, Eigen::VectorXd& gamma_2,
  Eigen::VectorXd& beta, Eigen::VectorXd& Y, Eigen::MatrixXd& X) {
    
    Eigen::VectorXd res = Y - X * beta.matrix();
    double arg1 = (N*1.0) / 2.0 + alpha_sigma;
    double arg2 = res.dot(res) / 2.0 + beta_sigma;
    double sigma2 = 1.0 / R::rgamma(arg1, 1.0 /arg2);
    return (sigma2);
}

// [[Rcpp::export]]
void update_s2_sc(int P, double alpha_s, double beta_s,
  Eigen::VectorXi& uni_types, Eigen::VectorXi& types, Eigen::VectorXd& b, 
  Eigen::VectorXd& w_weight, Eigen::VectorXd& s2) {

  int T = uni_types.size(); 
  //Eigen::VectorXd s2(T);
  for(int t = 0; t < T; t++){
    double sum_t = 0.0;
    double sum_b2 = 0.0;
    for(int p = 0; p < P; p++){
      if(types[p] == uni_types[t]){
        sum_t = sum_t + 1;
        sum_b2 = sum_b2 + pow(b[p], 2)/w_weight[p];
      }
    }
    double arg1 = sum_t / 2.0 + alpha_s;
    double arg2 = sum_b2 / 2.0 + beta_s;
    s2[t] = 1.0 / R::rgamma(arg1, 1.0/arg2);
  }
}

// [[Rcpp::export]]
Rcpp::List MCMC_sc(int seed, int burnInIter, int keepIter, 
  int G, int P, int printInt, double N, 
  double alpha1, double alpha2, double alpha_sigma,
  double alpha_s, double beta1, double beta2, 
  double beta_sigma, double beta_s,
  double pi1, double sigma2, 
  Eigen::VectorXi& uni_types, Eigen::VectorXi& types,
  Eigen::VectorXi& m_fg, Eigen::VectorXi& g_index, 
  Eigen::VectorXd& gamma_1, Eigen::VectorXd& gamma_1_rep_per_feature, 
  Eigen::VectorXd& gamma_2,
  Eigen::VectorXd& b, Eigen::VectorXd& pi2, Eigen::VectorXd& pi2_rep_per_feature, 
  Eigen::VectorXd& s2,
  Eigen::VectorXd& Y, Eigen::MatrixXd& X, 
  Eigen::VectorXd& w_weight, Eigen::VectorXd& v_weight, int s2Weighted_int, 
  int BernoulliWeighted_int, 
  Eigen::VectorXd& pi2_loglikeli, double pi2_prop_n, int MH_ind){

  int num_t = uni_types.size();

  Eigen::MatrixXd GAMMA_1(G, keepIter);
  Eigen::MatrixXd GAMMA_2(P, keepIter);
  Eigen::MatrixXd B(P, keepIter);
  Eigen::MatrixXd PI2(G, keepIter);
  Eigen::MatrixXd S2(num_t, keepIter);
  Eigen::VectorXd SIGMA2(keepIter); 
  Eigen::VectorXd PI1(keepIter);
  Eigen::MatrixXd PI2_LOGL_store(G, keepIter);

  set_seed(seed);

  Eigen::VectorXd beta = gamma_1_rep_per_feature.array() * 
    gamma_2.array() * b.array();
  Eigen::VectorXd s2_vec(P);
  for(int p = 0; p < P; p++){
    s2_vec[p] = s2[types[p] - 1];
  }
  int totalIter = burnInIter + keepIter;
  // MCMC iterations
  for(int i = 0; i < totalIter; i++){
    
    update_gamma1_sc(pi1, sigma2, G, P, g_index, 
      gamma_1, gamma_1_rep_per_feature, gamma_2, 
      b, Y, X); // beta is not updated after this step
    for(int p = 0; p < P; p++){
      gamma_1_rep_per_feature[p] = gamma_1[g_index[p]-1];
    }

    update_gamma2_sc(pi2_rep_per_feature, sigma2,
      P, g_index, m_fg, gamma_1_rep_per_feature, gamma_2, 
      b, beta, Y, X, v_weight, BernoulliWeighted_int); // beta is updated
    
    update_b_sc(pi2_rep_per_feature, sigma2, s2_vec, P, 
      gamma_1_rep_per_feature, gamma_2, 
      b, beta, Y, X, w_weight);  // beta is updated as well
    
    pi1 = update_pi1_sc(gamma_1, alpha1, beta1, G);

    update_pi2_sc(gamma_2, pi2, 
      g_index, m_fg, alpha2, beta2, G, P, 
      v_weight, pi2_loglikeli, pi2_prop_n, MH_ind);
    
    for(int p = 0; p < P; p++){
      pi2_rep_per_feature[p] = pi2[g_index[p] - 1];
    }

    sigma2 = update_sigma2_sc(N, alpha_sigma, beta_sigma, 
      gamma_1_rep_per_feature, gamma_2, beta, Y, X);
    
    update_s2_sc(P, alpha_s, beta_s, uni_types, types, b, 
      w_weight, s2);

    for(int p = 0; p < P; p++){
      s2_vec[p] = s2[types[p] - 1];
    }

    // store results
    if(i > burnInIter - 1){
      int store_index = i - burnInIter;
      GAMMA_1.col(store_index) = gamma_1;
      GAMMA_2.col(store_index) = gamma_2;
      B.col(store_index) = b;
      PI1[store_index] = pi1;
      PI2.col(store_index) = pi2;
      PI2_LOGL_store.col(store_index) = pi2_loglikeli;
      SIGMA2[store_index] = sigma2;
      S2.col(store_index) = s2;
    }

    if(i%printInt ==  0){
      cout<< "MCMC iteration: " << i << ".\n";
    }
  }

  return Rcpp::List::create(Rcpp::Named("GAMMA_1") = GAMMA_1,
                            Rcpp::Named("GAMMA_2") = GAMMA_2,
                            Rcpp::Named("B") = B,
                            Rcpp::Named("PI1") = PI1,
                            Rcpp::Named("PI2") = PI2,
                            Rcpp::Named("S2") = S2,
                            Rcpp::Named("SIGMA2") = SIGMA2,
                            Rcpp::Named("PI2_LOGL_store") = PI2_LOGL_store
                            ); 
}































