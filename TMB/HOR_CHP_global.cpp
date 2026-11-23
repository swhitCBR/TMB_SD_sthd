  #include <TMB.hpp>
  template<class Type>
  Type objective_function<Type>::operator() (){

  //////////////////
  // LOADING DATA //
  //////////////////

  DATA_MATRIX( XX );			// 
  DATA_MATRIX( XX_s );			// 
  DATA_MATRIX( Dat );			// 
  
  
  DATA_IVECTOR( ind_all );
  DATA_IVECTOR( f1_ind );
  DATA_IVECTOR( f2_ind );
  DATA_IVECTOR( f3_ind );
  
  DATA_IVECTOR( col_Dat1 );
  DATA_IVECTOR( col_Dat2 );
  DATA_IVECTOR( col_Dat3 );
  
  DATA_IVECTOR( col_S );
  DATA_IVECTOR( col_P11 );
  DATA_IVECTOR( col_P12 );
  
  DATA_IVECTOR( y123 );
  DATA_IVECTOR( y12 );
  DATA_IVECTOR( y13 );  

  DATA_IVECTOR( y1to7 );
  DATA_IVECTOR( y1to4 );
  DATA_IVECTOR( y1256 );  
  DATA_IVECTOR( y1357 );
  DATA_IVECTOR( y5to7 );
  DATA_IVECTOR( y347 );  
  DATA_IVECTOR( y246 );  

  DATA_VECTOR( lvec1 );			// 
  DATA_VECTOR( lvec2 );			// 
  DATA_INTEGER( P_lc_n );		// number of p and lambda linear combinations
  DATA_INTEGER( n_relgrps );
  DATA_IVECTOR( relgrps_ind );  

  
  ////////////////////////////////////////////////////
  // INITIALIZING PARAMETERS (Must be real numbers) //
  ////////////////////////////////////////////////////

  PARAMETER_VECTOR (S_pars); 
  // PARAMETER_VECTOR (yr_pars); 
  // PARAMETER (logSD_YR ); 
  PARAMETER (logSD_RELGRP ); 
  PARAMETER_VECTOR(RELGRP_err);
  // PARAMETER_VECTOR(OBS_err);
  // PARAMETER (logSD_OBS ); 
  PARAMETER_MATRIX (P_pars); // columns 1 and 2 (0 and 1 in C++) represent intercept and slope

 // Transforming SDs
  // Type Sigma_YR = exp(logSD_YR);  
  Type Sigma_RELGRP = exp(logSD_RELGRP);  
  // Type Sigma_OBS = exp(logSD_OBS);  
	
  int n_S_pars = S_pars.size();	
  // int n_yrs = yr_pars.size();
  int n_all = ind_all.size();
  
  // vector<Type> S_pars_tmp(yr_pars.size()+n_S_pars); // Result vector must be allocated in advance!
  // S_pars_tmp << S_pars(0),yr_pars,S_pars.segment(1,(n_S_pars-1));              // Syntax for concatination of v1 and v2


  // RANDOM EFFECT LIKELIHOODS
  // Type RE_YR_nll = 0;    // for(int ii=0; ii<n_yrs; ii++){	 				// RE_YR_nll -= dnorm(yr_pars(ii),Type(0),exp(logSD_YR),1); 		// } 

  Type RE_RELGRP_nll = 0;  
  for(int ii=0; ii<n_relgrps; ii++){	
				RE_RELGRP_nll -= dnorm(RELGRP_err(ii),Type(0),exp(logSD_RELGRP),1);
		} 
		
  // Type RE_OBS_nll = 0;  
  // for(int ii=0; ii<n_all; ii++){	
				// RE_OBS_nll -= dnorm(OBS_err(ii),Type(0),exp(logSD_OBS),1);
		// } 

  ADREPORT(Sigma_RELGRP)

  REPORT(RE_RELGRP_nll)
  // REPORT(S_pars_tmp)
  REPORT(n_S_pars)
  
  ////////////////////////////////////////////////////
  // INITIALIZING PARAMETERS (Must be real numbers) //
  ////////////////////////////////////////////////////

	// int f1_n = f1_ind.size();	
	// int f2_n = f2_ind.size();	
	// int f3_n = f3_ind.size();	
	int n = Dat.rows();	


	// Vector of survival linear combinations
	// vector<Type> lp_S = XX_s * S_pars_tmp;                     
	// vector<Type> lp_S_wo_err = XX_s * S_pars_tmp;                     // NOTE:matrix multiplication
	vector<Type> lp_S_wo_err = XX_s * S_pars;                     // NOTE:matrix multiplication

	
	// vector<Type> lp_S = lp_S_wo_err  + RELGRP_err(relgrps_ind) + OBS_err;                     
	vector<Type> lp_S = lp_S_wo_err  + RELGRP_err(relgrps_ind) ;                     
	
	// vector<Type> test = lp_S_wo_err + RELGRP_err(relgrps_ind);
	// REPORT(test)
	
	vector<Type> log_S = lp_S - log(1 + exp(lp_S));
	vector<Type> S = exp(log_S);


	// // Matrix of detection linear combinations
	matrix<Type> lp_P(n,P_lc_n); 
	
	for(int kk=0; kk<P_lc_n; kk++){	lp_P.col(kk) = XX * P_pars(col_P11,kk); } 
	
	// elementwise tranformation:  log_P = lp_P - log(1 + exp(lp_P))
	matrix<Type> log_pls1_exp_P = log(1 + exp(lp_P.array()).array());
	matrix<Type> log_P =lp_P - log_pls1_exp_P; 
	
	// not all of the linear combinations need a 'q' (e.g., lambda)
	matrix<Type> log_q = log(1 - exp(log_P.array()));
	matrix<Type> q = exp(log_q.array());

	// terms in Chi likelihood
	vector<Type> log_q11 = log_q.col(2);
	vector<Type> log_q12 = log_q.col(3);
	vector<Type> log_q13 = log_q12;  //P12 and P13 equated
	// vector<Type> log_q13 = log_q.col(4);

	vector<Type> q11 = exp(log_q11);
	vector<Type> q12 = exp(log_q12);
	vector<Type> q13 = exp(log_q13);


  //////////////////////////
  // LIKELIHOOD SUMMATION //
  //////////////////////////
  
	// common vectors and matrices
	
  	vector<Type> log_lvec1 = log(lvec1);
	vector<Type> log_lvec2 = log(lvec2);
  
  	/// formulation 3 ///

	// Data matrices f3
	vector<Type> f3_D2c1to7 = Dat(f3_ind,col_Dat3(y1to7)).rowwise().sum();
	vector<Type> f3_D2c1to4 = Dat(f3_ind,col_Dat3(y1to4)).rowwise().sum();
	vector<Type> f3_D2c1256 = Dat(f3_ind,col_Dat3(y1256)).rowwise().sum();
	vector<Type> f3_D2c1357 = Dat(f3_ind,col_Dat3(y1357)).rowwise().sum();
	vector<Type> f3_D2c5to7 = Dat(f3_ind,col_Dat3(y5to7)).rowwise().sum();
	vector<Type> f3_D2c347 = Dat(f3_ind,col_Dat3(y347)).rowwise().sum();
	vector<Type> f3_D2c246 = Dat(f3_ind,col_Dat3(y246)).rowwise().sum();

	vector<Type> f3_tmp0_v1 = log_lvec1(f3_ind)*f3_D2c1to7;
  	Type f3_tmp0 = sum(f3_tmp0_v1);
		
	vector<Type> f3_tmp1v1 = log_S(f3_ind)*f3_D2c1to7;
  	Type f3_tmp1 = sum(f3_tmp1v1);
	
	vector<Type> f3_tmp2v1 = log_P.col(2); 
	vector<Type> f3_tmp2v2 = f3_tmp2v1(f3_ind)*f3_D2c1to4;
	Type f3_tmp2 = sum(f3_tmp2v2);

	vector<Type> f3_tmp3v1 = log_P.col(3);
	vector<Type> f3_tmp3v2 = f3_tmp3v1(f3_ind)*f3_D2c1256;
	Type f3_tmp3 = sum(f3_tmp3v2);
	
	vector<Type> f3_tmp4v1 = log_P.col(3); // P13 is equated to P12
	vector<Type> f3_tmp4v2 = f3_tmp4v1(f3_ind)*f3_D2c1357;
	Type f3_tmp4 = sum(f3_tmp4v2);	

	vector<Type> f3_tmp5v1 = log_q11(f3_ind); 
	vector<Type> f3_tmp5v2 = f3_tmp5v1*f3_D2c5to7;
	Type f3_tmp5 = sum(f3_tmp5v2);	
	
	vector<Type> f3_tmp6v1 = log_q12(f3_ind); 
	vector<Type> f3_tmp6v2 = f3_tmp6v1*f3_D2c347;
	Type f3_tmp6 = sum(f3_tmp6v2);	

	vector<Type> f3_tmp7v1 = log_q13(f3_ind); 
	vector<Type> f3_tmp7v2 = f3_tmp7v1*f3_D2c246;
	Type f3_tmp7 = sum(f3_tmp7v2);	

	vector<Type> f3_log_chi0i = log(1 - exp(log_S(f3_ind))*lvec1(f3_ind)*(1-q11(f3_ind)*q12(f3_ind)*q13(f3_ind)));
	vector<Type> f3_tmp8v1 = Dat(f3_ind,col_Dat3(7));
	vector<Type> f3_tmp8v2 = f3_log_chi0i*f3_tmp8v1;
	Type f3_tmp8 = sum(f3_tmp8v2);

	Type f3_nll = -1*(f3_tmp0 + f3_tmp1 + f3_tmp2 + f3_tmp3 + f3_tmp4 + f3_tmp5 + f3_tmp6 + f3_tmp7 + f3_tmp8);	

	// REPORT(f3_tmp0)
	// REPORT(f3_tmp1)
	// REPORT(f3_tmp2)
	// REPORT(f3_tmp3)
	// REPORT(f3_tmp4)
	// REPORT(f3_tmp5)
	// REPORT(f3_tmp6)
	// REPORT(f3_tmp7)
	// REPORT(f3_tmp8)	
	
	/// formulation 1 ///
	
	// Data matrices f1
	vector<Type> f1_D2c12 = Dat(f1_ind,col_Dat1(y12)).rowwise().sum();
	vector<Type> f1_D2c13 = Dat(f1_ind,col_Dat1(y13)).rowwise().sum();
	vector<Type> f1_D2c123 = Dat(f1_ind,col_Dat1(y123)).rowwise().sum();
	
	vector<Type> f1_tmp0_1v1 = log_lvec1((f1_ind))*f1_D2c123;
  	Type f1_tmp0_1 = sum(f1_tmp0_1v1);
	
	vector<Type> f1_tmp0_2v1 = log_lvec2((f1_ind))*f1_D2c13;
  	Type f1_tmp0_2 = sum(f1_tmp0_2v1);
	
	vector<Type> f1_tmp1v1 = log_S(f1_ind)*f1_D2c123;
  	Type f1_tmp1 = sum(f1_tmp1v1);
	
	vector<Type> f1_tmp2v1 = log_P.col(0);
	vector<Type> f1_tmp2v2 = f1_tmp2v1(f1_ind)*f1_D2c12;
	Type f1_tmp2 = sum(f1_tmp2v2);
	
	vector<Type> f1_tmp3v1 = log_q.col(0);
	vector<Type> f1_tmp3v2 = Dat(f1_ind,col_Dat1(2));
	vector<Type> f1_tmp3v3 = f1_tmp3v1(f1_ind)*f1_tmp3v2;
	Type f1_tmp3 = sum(f1_tmp3v3);

	vector<Type> f1_tmp4v1 = log_P.col(1); // log_lambda vector
	vector<Type> f1_tmp4v2 = f1_tmp4v1(f1_ind)*f1_D2c13;
  	Type f1_tmp4 = sum(f1_tmp4v2);
	
	// form 1 specific vectors
	vector<Type> log_1minLamL2i =log(1- exp(f1_tmp4v1(f1_ind))*lvec2(f1_ind)); // lambda vector	
	vector<Type> q_tmp = q.col(0); // q1i
	vector<Type> f1_chi0i = 1 - S(f1_ind)*lvec1(f1_ind) + S(f1_ind)*lvec1(f1_ind)*q_tmp(f1_ind)*exp(log_1minLamL2i); // lambda vector	
	
	vector<Type> f1_tmp5v1 = Dat(f1_ind,col_Dat1(1));
	vector<Type> f1_tmp5v2 = f1_tmp5v1*log_1minLamL2i;
	Type f1_tmp5 = sum(f1_tmp5v2);	
	
	vector<Type> f1_tmp6v1 = Dat(f1_ind,col_Dat1(3));
	vector<Type> f1_tmp6v2 = log(f1_chi0i)*f1_tmp6v1;
	Type f1_tmp6 = sum(f1_tmp6v2);	

	// REPORT(f1_D2c12)
	REPORT(f1_D2c123)
	  
  REPORT(log_S)
	REPORT(f1_chi0i)
	REPORT(f1_tmp0_1)
	REPORT(f1_tmp0_2)
	REPORT(f1_tmp1)
	REPORT(f1_tmp2)
	REPORT(f1_tmp3)
	REPORT(f1_tmp4)
	REPORT(f1_tmp5)
	REPORT(f1_tmp6)
	
	// Form 1 Likelihood
	Type f1_nll = -1*(f1_tmp0_1 + f1_tmp0_2 + f1_tmp1 + f1_tmp2 + f1_tmp3 + f1_tmp4 + f1_tmp5 + f1_tmp6);	
	REPORT(f1_nll)

	/// formulation 2 ///


// form 1 specific vectors

	// vector<Type> f2_chi0i = 1 - exp(log_S(f2_ind))*lvec1(f2_ind)*(1-q11(f2_ind)*q12(f2_ind));
	// vector<Type> f2_log_chi0i = log(f2_chi0i);
	vector<Type> f2_log_chi0i = log(1 - exp(log_S(f2_ind))*lvec1(f2_ind)*(1-q11(f2_ind)*q12(f2_ind)));

	// Data matrices f2
	vector<Type> f2_D2c12 = Dat(f2_ind,col_Dat2(y12)).rowwise().sum();
	vector<Type> f2_D2c13 = Dat(f2_ind,col_Dat2(y13)).rowwise().sum();
	vector<Type> f2_D2c123 = Dat(f2_ind,col_Dat2(y123)).rowwise().sum();

	vector<Type> f2_tmp0v1 = log_lvec1((f2_ind))*f2_D2c123;
  	Type f2_tmp0 = sum(f2_tmp0v1);
	vector<Type> f2_tmp1v1 = log_S(f2_ind)*f2_D2c123;
  	Type f2_tmp1 = sum(f2_tmp1v1);

	vector<Type> f2_tmp2v1 = log_P.col(2);
	vector<Type> f2_tmp2v2 = f2_tmp2v1(f2_ind)*f2_D2c12;
	Type f2_tmp2 = sum(f2_tmp2v2);
	
	vector<Type> f2_tmp3v1 = log_P.col(3);
	vector<Type> f2_tmp3v2 = f2_tmp3v1(f2_ind)*f2_D2c13;
	Type f2_tmp3 = sum(f2_tmp3v2);
	
	vector<Type> f2_tmp4v1 = Dat(f2_ind,col_Dat2(2));
	vector<Type> f2_tmp4v2 = log_q11(f2_ind)*f2_tmp4v1;
	Type f2_tmp4 = sum(f2_tmp4v2);

	vector<Type> f2_tmp5v1 = Dat(f2_ind,col_Dat2(1));
	vector<Type> f2_tmp5v2 = log_q12(f2_ind)*f2_tmp5v1;
	Type f2_tmp5 = sum(f2_tmp5v2);

	vector<Type> f2_tmp6v1 = Dat(f2_ind,col_Dat2(3));
	vector<Type> f2_tmp6v2 = f2_log_chi0i*f2_tmp6v1;
	Type f2_tmp6 = sum(f2_tmp6v2);

	// Form 2 Likelihood
	Type f2_nll = -1*(f2_tmp0 + f2_tmp1 + f2_tmp2 + f2_tmp3 + f2_tmp4 + f2_tmp5 + f2_tmp6);	
	// REPORT(f2_nll)
	// REPORT(f2_tmp0)
	// REPORT(f2_tmp1)
	// REPORT(f2_tmp2)
	// REPORT(f2_tmp3)
	// REPORT(f2_tmp4)
	// REPORT(f2_tmp5)
	// REPORT(f2_tmp6)
	
	REPORT(RELGRP_err)

	// joint log likelihood (objective_function)
	Type nll = f1_nll + f2_nll + f3_nll + RE_RELGRP_nll;
	REPORT(nll)

  return nll;
  
  };
