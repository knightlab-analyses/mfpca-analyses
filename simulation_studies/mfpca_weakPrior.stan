functions {
  // Doesn't quite take the cholesky factor of `A`.
  // It replaces the diagonal elements with
  // exp(0.5 * o[n] + 2), to guarantee positive definitiveness
  // `o` should be unconstrained parameters.
  matrix pos_chol_offset(matrix A, vector o){
    int M = rows(A);
    matrix[M,M] L = rep_matrix(0, M, M);
    for (n in 1:M){
      real L_n_n;
      real Linv_n_n;
      for (m in (n+1):M){
  real L_m_n = A[m,n];
  for (k in 1:(n-1)){
    L_m_n -= L[m,k] * L[n,k];
  }
  L[m,n] = L_m_n;
      }
      L_n_n = exp(0.5 * o[n] + 2.0);
      L[n,n] = L_n_n;
      Linv_n_n = 1/L_n_n;
      for (m in (n+1):M){
  L[m,n] *= Linv_n_n;
      }
    }
    return L;
  }

  // `x` is a vector of unconstrained parameters
  // `sizes` is a vector of block sizes.
  // `M == sum(sizes)`; figure it's better to pass in as an arg rather than recalculate
  matrix constrain(vector x, int[] sizes, int M){
    int N = num_elements(sizes);
    int extracted_elements = sizes[1];
    matrix[M,M] L = rep_matrix(0, M, M);
    int row_pos = extracted_elements; // row position with respect to L
    for (r in 1:extracted_elements){ // do this one first because we can skip pos_chol_offset
      L[r,r] = exp(0.5*x[r] + 2.0);
    }
    // Alternatively, just create a matrix with undefined entries and then fill structural zeros with `0` in the loops.
    for (nr in 2:N){ // c for column block of L
      int rbs = sizes[nr]; // row block size
      matrix[rbs,rbs] D = rep_matrix(0, rbs, rbs);
      int col_pos = 0; // col position with respect to L
      for (nc in 1:nr-1){ // r for row
  int cbs = sizes[nc]; // column block size
  matrix[rbs,cbs] A;
  for (c in 1:cbs){
    A[:,c] = x[(extracted_elements+1):(extracted_elements+rbs)];
    extracted_elements += rbs;
  }
  L[(row_pos+1):(row_pos+rbs),(col_pos+1):(col_pos+cbs)] = A;
  D -= A * A';
  col_pos += cbs;
      }
      L[(row_pos+1):(row_pos+rbs),(row_pos+1):(row_pos+rbs)] = pos_chol_offset(D, x[(extracted_elements+1):(extracted_elements+rbs)]);
      extracted_elements += rbs;
      row_pos += rbs;
    }
    return L;
  }

  int unconstrained_len(int[] sizes){
    int s = 0;
    int N = num_elements(sizes);
    for (nc in 1:N){
      int sz = sizes[nc];
      s += sz;
      for (nr in (nc+1):N){
  s += sz * sizes[nr];
      }
    }
    return s;
  }

  real triangle_logdet(matrix L){
    real ld = 0.0;
    int N = rows(L);
    for (n in 1:N){
      ld += log(L[n,n]);
    }
    return ld;
  }
  
}

data {
  int <lower=1> N; // number of total unique subjects
  int <lower=1> P; // number of blocks
  int <lower=1> M; // number of total visits across all subjects
  int <lower=1> K[P]; // number of principal components for P blocks
  int <lower=1> Q[P]; // number of basis for the cubic spline basis for P blocks

  int <lower=1> subject_starts[N]; // starting index for each subjects' visit
  int <lower=1> subject_stops[N]; // ending index for each subjects' visit
  int <lower=1> subject_starts_each[N*P]; // starting index for each subjects' visit in each block
  int <lower=1> subject_stops_each[N*P]; // ending index for each subjects' visit in each block

  int V[N, P]; // visits matrix for each subject in each block 
  matrix[M, sum(Q)] B; // transpose of sparse spline basis 
  vector[M] Y; // record response at each available time point for each subject 

  // below are indexing needed for Theta estimation
  int <lower=1> len_Theta; // number of non-zero values in Theta
  int Theta_idx[P];
  int Q_idx[P];
  int K_idx[P];
}

transformed data{
  int sum_K = sum(K);
  int nu_len = unconstrained_len(K);     
}

parameters { 
  vector[sum(Q)] theta_mu;
  vector[len_Theta] vals_Theta; // non-zero values in Theta
  real<lower=0> sigma_eps[P];
  vector[nu_len] nu;
  matrix[sum(K),N]  alpha;
}

transformed parameters{
  vector[M] cov_y_vector; // to reformat sigma_eps;
  matrix[sum(Q), sum(K)] Theta;
  matrix[sum_K,sum_K] VCOV_factor = constrain(nu, K, sum_K);

  {
    int idx_eps; // for reformatting sigma_eps; 

    // reformatting sigma_eps;
    idx_eps = 0;
    for (nn in 1:N){
      for (p in 1:P){
        idx_eps = idx_eps + 1;
        cov_y_vector[subject_starts_each[idx_eps]:subject_stops_each[idx_eps]] = rep_vector(sigma_eps[p], 
                                                                                            V[nn,p]);
      }
    }

    // estimate non-zero values in Theta only
    Theta = rep_matrix(0, sum(Q), sum(K));
    Theta[1:Q_idx[1], 1:K_idx[1]] = to_matrix(vals_Theta[1:Theta_idx[1]], Q[1], K[1]);
    for (p in 2:P){
      Theta[Q_idx[p-1]+1:Q_idx[p], K_idx[p-1]+1:K_idx[p]] = 
              to_matrix(vals_Theta[Theta_idx[p-1]+1: Theta_idx[p]], Q[p], K[p]);
    }

  } // close modeling
}

model {
  theta_mu ~ normal(0, 1);
  vals_Theta ~ normal(0, 1);
  sigma_eps ~ cauchy(0, 1);

  // nu ~ normal(0,1);
  // likelihood
  // quadratic of multivariate normal
  target += -0.5*sum(rows_dot_self(mdivide_left_tri_low(VCOV_factor, alpha)));
  // -N * log(|VCOV_factor|)
  target += -N * triangle_logdet(VCOV_factor);

  // posterior draws
  for(n in 1:N){
    Y[subject_starts[n]: subject_stops[n]] ~ normal(B[subject_starts[n]: subject_stops[n]]*
                     theta_mu + B[subject_starts[n]: subject_stops[n]] * Theta * alpha[, n], cov_y_vector[subject_starts[n]: subject_stops[n]]);
  } 

}

generated quantities {
  cov_matrix[sum_K] cov_alpha = tcrossprod(VCOV_factor);
}



