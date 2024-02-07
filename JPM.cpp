#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data and design variables
  DATA_VECTOR(SUVR); // Outcome
  
  // Fixed effect design
  DATA_ARRAY(X_tracer); // Tracer design matrix
  DATA_VECTOR(x_centaur); // CenTauR values for groups
  DATA_ARRAY(X_h2h_subj); // CenTauR design matrix
  
  // Random effect design
  DATA_ARRAY(Z_subj); // Subject-level random-effect variance design matrix
  
  // Variance design
  DATA_ARRAY(V); // Variance design matrix for residual error
  DATA_ARRAY(V_subj); // Variance design matrix for random effect 
  
  // Parameters
  PARAMETER_VECTOR(slope_tracer);  
  PARAMETER_VECTOR(intercept_tracer);
  PARAMETER_VECTOR(centaur_h2h_subj);  
  PARAMETER_VECTOR(centaur_subj);  
  PARAMETER_VECTOR(log_sigma_tracer);
  PARAMETER_VECTOR(log_sigma_subj_group);  
  
  // Data sizes  
  int n = SUVR.size(); // Number of observations
  int n_tracer = X_tracer.dim[1]; // Number of tracers
  int n_h2h = X_h2h_subj.dim[1]; // Number of CenTauR groups
  int n_var = V.dim[1]; // Number of variance parameters
  int n_anchor = Z_subj.dim[1]; // Number of subjects with random effects
  int n_anchor_var = V_subj.dim[1]; // Number of subject-level variance parameters
  
  
  // Variance parameters
  vector<Type> sigma_tracer(log_sigma_tracer.size());
  sigma_tracer = exp(log_sigma_tracer);
  vector<Type> sigma_subj_group(log_sigma_subj_group.size());
  sigma_subj_group = exp(log_sigma_subj_group);
  
  // Intermediate vectors for picking out right parameters
  vector<Type> slope(n);
  vector<Type> intercept(n);
  vector<Type> centaur(n);  
  vector<Type> random(n);  
  vector<Type> sd(n);
  vector<Type> subject_sd(n_anchor);
  
  // Conditional fit
  vector<Type> yfit(n);
  
  // Negative log likelihood
  Type nll = Type(0.0);
  
  // Loop over anchor point subjects for random effect 
  // contribution to log likelihood
  for (int j = 0; j < n_anchor; j++) { 
    subject_sd(j) = 0;
    for (int k = 0; k < n_anchor_var; k++) {
      subject_sd(j) += V_subj(j, k) * sigma_subj_group(k);
    }
    // Add random effect contribution to negative log likelihood
    nll += -dnorm(centaur_subj(j), Type(0.0), subject_sd(j), true);
  }
  
  // Loop over all observations
  for (int i = 0; i < n; i++) { 
    // Calculate tracer CenTauR slope and intercept for observation
    slope(i) = 0.0;
    intercept(i) = 0.0;
    for (int j = 0; j < n_tracer; j++) {
      slope(i) += X_tracer(i, j) * slope_tracer(j);
      intercept(i) += X_tracer(i, j) * intercept_tracer(j);
    }
    
    // Initialize CenTauR group (0 for h2h)
    centaur(i) = x_centaur(i);
    
    // Add h2h CenTauR fixed effects
    for (int j = 0; j < n_h2h; j++) {
      centaur(i) += X_h2h_subj(i, j) * centaur_h2h_subj(j);
    }
    
    // Random effect contribution (anchor point groups only)
    random(i) = 0.0;
    for (int j = 0; j < n_anchor; j++) {
      random(i) += Z_subj(i, j) * centaur_subj(j);
    }
    
    // Calculate the conditional fit
    yfit(i) = slope(i) * (centaur(i) + random(i)) + intercept(i);
    
    // Calculate residual variance of observation
    sd(i) = 0.0;
    
    // Multiply sd with slope to ensure 
    // equal variance on CenTauR scale
    for (int j = 0; j < n_var; j++) {
      sd(i) += slope(i) * V(i, j) * sigma_tracer(j);
    }
    

    // Add contribution to negative log likelihood
    nll += -dnorm(SUVR(i), yfit(i), sd(i), true);
  }
  
  return nll;
}
