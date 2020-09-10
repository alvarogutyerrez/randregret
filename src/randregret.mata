

mata:
	void randregret_LL(transmorphic scalar M, real scalar todo,
	real rowvector b, real scalar lnf,
	real rowvector g, real matrix H)
 {	
	// variables declaration	
	real matrix panvar  
	real matrix paninfo 
	real scalar npanels 

	real scalar n 
	real scalar i 
	
	real scalar mu_raw
	real scalar mu
	
	real matrix Y 
	real matrix X 
	
	real matrix r_i 
	real matrix x_n 
	real matrix y_n 
	
	cons_demanded		= st_global("cons_demanded")
	regret_fn			= st_global("regret_fn") 	
	upmu 				= strtoreal(st_global("uppermu")) 	
	
	// variables creation
	Y = moptimize_util_depvar(M, 1) 	
	X= moptimize_init_eq_indepvars(M,1)	
	// Coefficients
	// regret function
	id_regret_eq=moptimize_util_eq_indices(M,1)
	b_regret = b[|id_regret_eq|]
	
	if (cons_demanded=="YES"){
	    
		ASC = moptimize_init_eq_indepvars(M,2) 
		id_ASC_eq=moptimize_util_eq_indices(M,2)
		b_ASC = b[|id_ASC_eq|]
		
		if (regret_fn=="mu"){
		    param_eq_mu  = moptimize_util_eq_indices(M, 3)
			mu_raw = b[|param_eq_mu|]
			mu = (exp(mu_raw) / (1 + exp(mu_raw))) *upmu
		}
		else if (regret_fn=="gene") {
		    param_eq_gamma  = moptimize_util_eq_indices(M, 3)
			gamma_raw = b[|param_eq_gamma|]
			gamma = invlogit(gamma_raw)  
		}
	}   
	else if (cons_demanded=="NO"){
	    if (regret_fn=="mu"){
		    param_eq_mu  = moptimize_util_eq_indices(M, 2)
			mu_raw = b[|param_eq_mu|]
			mu = invlogit(mu_raw) *upmu	
			
		}
		else if (regret_fn=="gene") {
		    param_eq_gamma  = moptimize_util_eq_indices(M, 2)
			gamma_raw = b[|param_eq_gamma|]
			gamma = invlogit(gamma_raw)  
			// (exp(gamma_raw) / (1 + exp(gamma_raw)))
		}		
	}


	// group; choice situations id.
	st_view(panvar = ., ., st_global("group_mata"))

	paninfo = panelsetup(panvar, 1) 	
	npanels = panelstats(paninfo)[1]
	
	st_numscalar("__n_cases", npanels)


	
	lnfj = J(npanels, 1, 0) 

	for(n=1; n <= npanels; ++n) {
		x_n = panelsubmatrix(X, n, paninfo) 
		y_n = panelsubmatrix(Y, n, paninfo) 
		// ------regret function-------//		
		if (regret_fn=="classic"){
			regret_n=RRM_log(x_n,b_regret,1,1)
		}
		else if (regret_fn=="mu"){
		    regret_n=RRM_log(x_n,b_regret,mu,1)
		}
		else if (regret_fn=="gene"){
		    regret_n=RRM_log(x_n,b_regret,1,gamma)
		}
		// ------  ASC  ------- //		
		if (cons_demanded=="YES") { 
			asc_n 	= panelsubmatrix(ASC, n, paninfo) 
			ASC_prod= asc_n*b_ASC'
			regret_i =rowsum(regret_n) + ASC_prod 
		}
		else if (cons_demanded=="NO"){
			regret_i =rowsum(regret_n) 
		}
		R_exp = exp(-regret_i)
		p_i = colsum((R_exp :* y_n))  / colsum(R_exp)
		lnfj[n] = ln(p_i)
	}
	lnf = moptimize_util_sum(M, lnfj)
}
end


mata:
real matrix RRM_log(real matrix 	x_n, 
					real rowvector 	b, 
					real scalar 	mu, 
					real scalar 	gamma)	
{
	real scalar i, j ,r_m
	real matrix regret_n 
	
	regret_n = J(rows(x_n), cols(x_n), 0)
	for(i=1; i <= rows(x_n); ++i) { 
		for(j=1; j <= rows(x_n); ++j) { 
			if (i!=j) { 
				r_i = ln(gamma :+ exp( (b :/ mu) :* ( x_n[j , . ] :-  x_n[i, . ]))) 				
				regret_n[i, . ] = regret_n[i, . ] :+ mu*r_i 		
				} 
			}	  
		}
	return(regret_n)
}
end


// randregretpred function.
mata:
real matrix pbb_pred(real matrix X,real matrix ASC, real colvector panvar )
{
	b_hat  		= st_matrix("b_hat")
 	mu_hat 		= st_matrix("mu_hat")
	gamma_hat	= st_matrix("gamma_hat")	
		
	ASC_hat  	= st_matrix("ASC_hat") 
	ASC_option 	= st_global("e(ASC)") 

	rrmfn 		= st_global("e(rrmfn)") 	
	proba		= st_local("proba") 

	paninfo = panelsetup(panvar, 1) 
	npanels = panelstats(paninfo)[1]
	
	for(n=1; n <= npanels; ++n) { 	
		x_n = panelsubmatrix(X, n, paninfo)
		// Generating regret part
		if (rrmfn != "pure") {
			hat_regret_n = RRM_log(x_n,b_hat,mu_hat,gamma_hat) 
		}
		else {
			hat_regret_n = -b_hat:*x_n
		}
		// Adding ASC fo the regret function
		if (ASC_option=="YES") {
			asc_n = panelsubmatrix(ASC, n, paninfo)  
			ASC_prod= asc_n*ASC_hat'
			regret_i =rowsum(hat_regret_n) + ASC_prod 
		}
		else{
			regret_i =rowsum(hat_regret_n) 
		}
		// Exp(-regret) 
		if (rrmfn != "pure") R_exp_hat = exp(-regret_i)
		else  R_exp_hat = exp(-regret_i) 
		// Computing Prediction
		p_hat_i = R_exp_hat   :/ colsum(R_exp_hat)
		// Probability prediction: proba option [default]		
		if (proba != "") {
			if (n==1) 	p_hat	= p_hat_i	
			else 		p_hat	= p_hat  \ p_hat_i	
			}
		else{ // Linear prediction: xb option
			if (n==1) 	regret_hat	= regret_i	
			else 		regret_hat	= regret_hat  \ regret_i	
			}	
	}	
	// Output prediction: 
	if (proba != "") return(p_hat)
	else return(regret_hat)	
}
end




//  Analytical Gradients utilities

// dR_db_fn: partial R_in w.r.t beta
mata:
real matrix dR_db_fn(real matrix x_n, real rowvector b, real scalar mu, real scalar gamma )
{
	// Equation (18), (19) or (23) on the article (depending on the inputs).
	real scalar i, j 
	real matrix dR_db 

	real matrix num
	real matrix den
	
	rows_x_n =rows(x_n) 
	dR_db = J(rows(x_n), cols(x_n), 0)
	
	for(i=1; i <= rows_x_n; ++i) { 
		for(j=1; j <= rows_x_n; ++j) { 
			if (i!=j) { 
				num = ( exp( (b :/ mu) :* ( x_n[j , . ] :-  x_n[i, . ])) ) :* ( x_n[j , . ] :-  x_n[i, . ])
				den = mu * (gamma :+  exp(  (b :/ mu) :* ( x_n[j , . ] :-  x_n[i, . ]))   ) 		
				dR_db[i, . ] = dR_db[i, . ] :+  (num :/ den)  		
				} 
			}	  
		}
		return(dR_db)
}
end

// dR_db_fn: partial R_in w.r.t mu
mata:
real matrix dR_dmu_fn(real matrix x_n, real rowvector b, real scalar mu)
{
   	// Equation (26) replaced into equation (25) on the article. 
	real scalar i, j 
	real matrix dR_dmu 

	real matrix num
	real matrix den
	rows_x_n =rows(x_n) 	
	dR_dmu = J(rows_x_n, cols(x_n), 0)
	for(i=1; i <= rows_x_n; ++i) { 
		for(j=1; j <= rows_x_n; ++j) { 
			if (i!=j) { 
		num = exp( (b :/ mu) :* ( x_n[j , . ] :-  x_n[i, . ]))  :* ( x_n[j , . ] :-  x_n[i, . ])	:* b 
		den = (mu*mu) :*  (1 :+  exp( (b :/ mu) :* ( x_n[j , . ] :-  x_n[i, . ])) )
		dR_dmu[i, . ] = dR_dmu[i, . ] :+  (num :/ den)  		
				} 
			}	  
		}
		// the derivative consider all the alternatives, 
		// therefore we need a sum over them (rows)
		dR_dmu = rowsum(-dR_dmu)
		return(dR_dmu)
}
end

// dR_db_fn: partial R_in w.r.t gamma
mata:
real matrix dR_dgamma_fn(real matrix x_n, real rowvector b, real scalar gamma)
{
   	// Equation (21) on the article 
	real scalar i, j 
	real matrix dR_dgamma 

	real matrix num
	real matrix den
	
	rows_x_n =rows(x_n) 	
	dR_dgamma = J(rows_x_n, cols(x_n), 0)
	
	for(i=1; i <= rows_x_n; ++i) { 
		for(j=1; j <= rows_x_n; ++j) { 
			if (i!=j) { 
		partial =   (gamma :+  exp( b  :* ( x_n[j , . ] :-  x_n[i, . ]))):^(-1)
		dR_dgamma[i, . ] = dR_dgamma[i, . ] :+  partial		
				} 
			}	  
		}
		// the derivative consider all the alternatives, 
		// therefore we need a sum over them (rows)
		dR_dgamma = rowsum(dR_dgamma)
		return(dR_dgamma)
}
end




// Robust & Cluster Standar Error Computations
mata:
real matrix std_errs(real matrix 	D, 
					 real colvector Y,
					 real matrix 	X,
					 real matrix 	ASC,
					 real colvector panvar )
{
	external ASC_hat
	external b_hat
	external mu_star_hat
	external gamma_star_hat
	external D
	external cluster_var_mata
	
	// st_numscalar("r(unique_value)")
 
	cluster_option	= st_global("cluster_mata")
	ASC_option		= st_global("cons_demanded")
	rrmfn			= st_global("regret_fn") 	
	upmu 			= strtoreal(st_global("uppermu")) 	
	
	

	// id each choice situation
	paninfo = panelsetup(panvar, 1)     
	npanels = panelstats(paninfo)[1] 
	for(n=1; n <= npanels; ++n) { 
		
		// Extract only the submatrix of individual n
		x_n = panelsubmatrix(X, n, paninfo) 
		y_n = panelsubmatrix(Y, n, paninfo) 
		if (ASC_option == "YES"){
			asc_n = panelsubmatrix(ASC, n, paninfo) 			
		} 
		else{
			asc_n = J(1,1,0)
			ASC_hat = J(1,1,0)
		} 
		if (rrmfn=="classic"){
			// If not ASC, this is just zero.
			ASC_prod=  asc_n*ASC_hat' 
			// estimated regret
			regret_n=RRM_log(x_n,b_hat,1,1) 
			R_in =rowsum(regret_n) :+ ASC_prod
			P_in =(exp(-R_in))  :/ colsum(exp(-R_in)) 
			
			// gradient of regret params
			dR_db= dR_db_fn(x_n, b_hat,1,1)
			grad  = colsum(-(y_n - P_in):*dR_db)
			
			if  (ASC_option == "YES"){
				// if ASC then added to the gradients
				grad_dASC = colsum(-(y_n - P_in):*asc_n)
				grad = grad,  grad_dASC  
			} 
				
		}
		else if (rrmfn=="mu"){
			mu_hat = invlogit(mu_star_hat)*upmu
			// If not ASC, this is just zero.
			ASC_prod=  asc_n*ASC_hat' 
			// estimated regret
			regret_n=RRM_log(x_n,b_hat,mu_hat ,1) :+  ASC_prod
			R_in =rowsum(regret_n) 
			P_in =(exp(-R_in))  :/ colsum(exp(-R_in)) 

			// gradient of regret params
			dR_db= dR_db_fn(x_n, b_hat,mu_hat ,1)
			grad_b  = colsum(-(y_n - P_in) :* mu_hat  :* dR_db)
	  
			// gradient of mu
			dR_dmu = dR_dmu_fn(x_n,b_hat,mu_hat )
			scale_factor = ((mu_hat)*(upmu-mu_hat))/(upmu)
			grad_mu=colsum(-(y_n - P_in) :* (mu_hat  :* dR_dmu + R_in :/ mu_hat ))*scale_factor
			
			// full gradient vector
			grad =grad_b,grad_mu

			if  (ASC_option == "YES"){
				// if ASC then added to the gradients
				grad_dASC = colsum(-(y_n - P_in):*asc_n)
				grad = grad_b,  grad_dASC, grad_mu 
			} 		
			
		}
		else if (rrmfn=="gene"){			
			gamma_hat = invlogit(gamma_star_hat)
			// If not ASC, this is just zero.
			ASC_prod=  asc_n*ASC_hat' 
			// estimated regret
			regret_n=RRM_log(x_n,b_hat,1,gamma_hat) :+  ASC_prod
			R_in =rowsum(regret_n) 
			P_in =(exp(-R_in))  :/ colsum(exp(-R_in)) 
			
			// gradient of regret params
			dR_db= dR_db_fn(x_n, b_hat,1,gamma_hat)
			grad_b  = colsum(-(y_n - P_in):*dR_db)

			// gradient of gamma		
			dR_dgamma = dR_dgamma_fn(x_n, b_hat,gamma_hat)
			scale_factor =gamma_hat * (1- gamma_hat)
			grad_gamma  = colsum(-(y_n - P_in):*dR_dgamma*scale_factor)
			
			// full gradient vector
			grad =grad_b,grad_gamma
			
			if  (ASC_option == "YES"){
				// if ASC then added to the gradients
				grad_dASC = colsum(-(y_n - P_in):*asc_n)
				grad = grad_b, grad_dASC,grad_gamma   
			}	 			
		}
		
		if (n==1) U = grad    
		else	  U = U \ grad    
	}

	if (cluster_option!=""){

		cluinfo = panelsetup(cluster_var_mata, 1)     
		nclusters = panelstats(cluinfo)[1] 
		st_numscalar("__n_clusters", nclusters)

		for(c=1; c <= nclusters; ++c) { 
			U_c = colsum(panelsubmatrix(U, c, cluinfo)) 
			if (c==1) U_cluster = U_c    
			else	  U_cluster= U_cluster \ U_c  	
		}
		 V_r = D*cross(U_cluster,U_cluster)*D*(nclusters/(nclusters-1))
	}
	else V_r = D*cross(U,U)*D*(npanels/(npanels-1))

	
	return(V_r)
}
end

/*

capture mata: mata drop randregret_LL() RRM_log() pbb_pred() dR_db_fn() /*
*/ dR_dmu_fn() dR_dgamma_fn()  std_errs()


*/


 