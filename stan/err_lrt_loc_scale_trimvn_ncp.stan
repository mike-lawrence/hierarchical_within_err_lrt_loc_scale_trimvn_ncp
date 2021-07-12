//aria: compile=1
//aria: compile_debug=1
//aria: run_debug=0 #because the auto-generated debug data doesn't work
//aria: syntax_ignore = c('has no priors.')
functions{
	// flatten_lower_tri: function that returns the lower-tri of a matrix, flattened to a vector
	vector flatten_lower_tri(matrix mat) {
		int n_cols = cols(mat);
		int n_uniq = (n_cols * (n_cols - 1)) %/% 2;
		vector[n_uniq] out ;
		int i = 1;
		for(c in 1:(n_cols-1)){
			for(r in (c+1):n_cols){
				out[i] = mat[r,c];
				i += 1;
			}
		}
		return(out) ;
	}
}
data{

	// n: number of subj
	int<lower=1> n ;

	// k: number of within-subj predictors
	int<lower=1> k ;

	// m_lrt: number of lrt observations
	int<lower=1> m_lrt ;

	// lrt: lrt on each trial
	vector[m_lrt] lrt ;

	// m_err: number of error summary cells
	int<lower=1> m_err ;

	// err_cell_count_total: number of error observations per cell of the summary
	int err_cell_count_total[m_err] ;

	// err_cell_count_errors: number of errors per cell of the summary
	int err_cell_count_errors[m_err] ;

	// w: unique entries in the within predictor matrix
	matrix[k,k] w ;

	// lrt_row: index of each lrt in flattened subject-by-condition value matrix
	array[m_lrt] int lrt_row ;

	// err_row: index of each err summary in flattened subject-by-condition value matrix
	array[m_err] int err_row ;

	// is_intercept: binary indicator of columns that reflect intercept parameters
	array[k,1] int<lower=0,upper=1> is_intercept ;

}
transformed data{

	// obs_lrt_mean: mean lrt value
	real obs_lrt_mean = mean(lrt) ;

	// obs_lrt_sd: sd of lrts
	real obs_lrt_sd = sd(lrt) ;

	// lrt_: lrt scaled to have zero mean and unit variance
	vector[m_lrt] lrt_ = (lrt-obs_lrt_mean)/obs_lrt_sd ;

	// compute observed intercept for the error data
	real obs_err_intercept = logit(mean(to_vector(err_cell_count_errors)./to_vector(err_cell_count_total))) ;

	//wn: an array where each element is a row from w repeated for n rows
	array[k] matrix[n,k] wn ;
	for(i_k in 1:k){
		wn[i_k] = rep_matrix(w[i_k],n) ;
	}

	// which_r: matrix-like array to facilitate extracting the pairwise correlations
	//    from the lower-tri-vector parameter representation
	array[k,k] int which_r ;
	int which_r_num = 0 ;
	for(i_k in 1:(k-1)){
		for(j_k in (i_k+1):k){
			which_r_num = which_r_num + 1;
			which_r[i_k,j_k] = which_r_num ;
		}
	}

}
parameters{

	//for parameters below, trailing underscore denotes that they need to be un-scaled in generated quantities

	// z_m_: mean (across subj) for each coefficient
	matrix[k,3] z_m_ ;

	// z_s_: sd (across subj) for each coefficient
	matrix<lower=0>[3,k] z_s_ ;

	// z_z: a helper variable for implementing non-centered parameterization of z_
	array[3] matrix[k,n] z_z ;

	// r_: population-level correlations (on cholesky factor scale) amongst within-subject predictors
	array[k] cholesky_factor_corr[3] r_ ;

}
model{
	////
	// Priors ----
	////

	// normal(0,1) priors on all means
	to_vector(z_m_) ~ std_normal() ;

	// weibull priors on all sds
	to_vector(z_s_) ~ weibull(2,1) ;

	for(i_k in 1:k){
		// flat lkj prior on correlations
		r_[i_k] ~ lkj_corr_cholesky(1) ;
	}

	////
	// Mid-hierarchy structure ----
	////

	for(i in 1:3){
		to_vector(z_z[i]) ~ std_normal() ; // with TP implies multi_normal(z_m_,z_s_,r_)
	}

	{// local environment to avoid saving intermediate variables

		// compute z_ (non-centered parameterization)
		array[k] matrix[n,3] z_ ;
		matrix[n,k] err_z_ ;
		matrix[n,k] mrt_z_ ;
		matrix[n,k] srt_z_ ;
		for(i_k in 1:k){
			z_[i_k] = (
				rep_matrix(z_m_[i_k],n)
				+ transpose(
					diag_pre_multiply(
						z_s_[,i_k]
						, r_[i_k]
					)
					* append_row(
						append_row(
							z_z[1][i_k]
							, z_z[2][i_k]
						)
						, z_z[3][i_k]
					)
				)
			) ;
			err_z_[,i_k] = z_[i_k][,1];
			mrt_z_[,i_k] = z_[i_k][,2];
			srt_z_[,i_k] = z_[i_k][,3];
		}

		// condition values implied by by-subject coefficients & contrasts
		matrix[n,k] err_z_dot_w ;
		matrix[n,k] mrt_z_dot_w ;
		matrix[n,k] srt_z_dot_w ;
		for(i_k in 1:k){
			err_z_dot_w[,i_k] = rows_dot_product(err_z_ , wn[i_k]) ;
			mrt_z_dot_w[,i_k] = rows_dot_product(mrt_z_ , wn[i_k]) ;
			srt_z_dot_w[,i_k] = rows_dot_product(srt_z_ , wn[i_k]) ;
		}

		////
		// Observation-level structure ----
		////

		// lrt
		lrt_ ~ normal(
			to_vector(mrt_z_dot_w)[lrt_row]
			, sqrt(exp(to_vector(srt_z_dot_w)))[lrt_row]
		) ;

		// err
		err_cell_count_errors ~ binomial_logit(
			err_cell_count_total
			, (obs_err_intercept+to_vector(err_z_dot_w))[err_row]
		) ;
	}
}
generated quantities{

	// group-level quantities ----

	// extract and unscale/unshift
	array[k] vector[3] r ;
	vector[k] err_m ;
	vector[k] mrt_m ;
	vector[k] srt_m ;
	vector[k] err_s ;
	vector[k] mrt_s ;
	vector[k] srt_s ;
	for(i_k in 1:k){
		r[i_k] = flatten_lower_tri( multiply_lower_tri_self_transpose(r_[i_k]) ) ;
		err_m[i_k] = z_m_[i_k,1] ;
		mrt_m[i_k] = z_m_[i_k,2]*obs_lrt_sd ;
		srt_m[i_k] = z_m_[i_k,3]*obs_lrt_sd ;
		err_s[i_k] = z_s_[1,i_k] ;
		mrt_s[i_k] = z_s_[2,i_k]*obs_lrt_sd ;
		srt_s[i_k] = z_s_[3,i_k]*obs_lrt_sd ;
		if(is_intercept[i_k,1]==1){
			err_m[i_k] += obs_err_intercept ;
			mrt_m[i_k] += obs_lrt_mean ;
		}
	}

	// condition values implied by mean coefficients & contrasts
	vector[k] err_elogodds ;
	vector[k] lrt_emean ;
	vector[k] lrt_esd ;
	for(i_k in 1:k){
		err_elogodds[i_k] = dot_product(err_m,w[i_k]);
		lrt_emean[i_k] = dot_product(mrt_m,w[i_k]);
		lrt_esd[i_k] = sqrt(exp(dot_product(srt_m,w[i_k])));
	}

	// subject-level quantities ----

	// leave commented-out to avoid excess write during sampling (can compute from samples later)

	// matrix[n,k] err_q_ ;
	// matrix[n,k] mrt_q_ ;
	// matrix[n,k] srt_q_ ;
	// {
	// 	// compute z_ (non-centered parameterization)
	// 	array[k] matrix[n,3] z_ ;
	// 	for(i_k in 1:k){
	// 		z_[i_k] = (
	// 			rep_matrix(z_m_[i_k],n)
	// 			+ transpose(
	// 				diag_pre_multiply(
	// 					z_s_[,i_k]
	// 					, r_[i_k]
	// 				)
	// 				* append_row(
	// 					append_row(
	// 						z_z[1][i_k]
	// 						, z_z[2][i_k]
	// 					)
	// 					, z_z[3][i_k]
	// 				)
	// 			)
	// 		) ;
	// 		err_q_[,i_k] = z_[i_k][,1] ;
	// 		mrt_q_[,i_k] = z_[i_k][,2] * obs_lrt_sd;
	// 		srt_q_[,i_k] = z_[i_k][,3] * obs_lrt_sd;
	// 		if(is_intercept[i_k,1]==1){
	// 			err_q_[,i_k] += obs_err_intercept ;
	// 			mrt_q_[,i_k] += obs_lrt_mean ;
	// 		}
	// 	}
	// }
	// // condition values implied by by-subject coefficients & contrasts
	// matrix[n,k] err_q_dot_w ;
	// matrix[n,k] mrt_q_dot_w ;
	// matrix[n,k] srt_q_dot_w ;
	// for(i_k in 1:k){
	// 	err_q_dot_w[,i_k] = rows_dot_product(err_q_ , wn[i_k]) ;
	// 	mrt_q_dot_w[,i_k] = rows_dot_product(mrt_q_ , wn[i_k]) ;
	// 	srt_q_dot_w[,i_k] = rows_dot_product(srt_q_ , wn[i_k]) ;
	// }


}
