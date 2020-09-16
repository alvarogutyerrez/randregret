cscript "stub of test script"
display "to be filled in"


clear all
set more off
which randregret
which randregret_pure

*Simulated data to recover parameters  

*RRM Models:

* classic 	= 2 attributes 

* gene 		= 2 attributes + gamma 


* mu   		= 2 attributes + mu  

* pure 		= 2 (pure)attributes  


global n_obs = 2500
global n_rep = 100


mata:
void data_gen() 
{
		st_view(X = ., ., "x*")
		st_view(panvar = ., ., "id") 
		betas  = st_matrix("betas_st")
		gamma  = strtoreal(st_global("real_gamma")) 
		mu 	   = strtoreal(st_global("real_mu")) 
		
		
		pure_model 	   = st_global("pure_model") 
		
		
		paninfo = panelsetup(panvar, 1)     
        npanels = panelstats(paninfo)[1] 
		
		for(n=1; n <= npanels; ++n) { 
            x_n = panelsubmatrix(X, n, paninfo) 
			if (pure_model =="YES"){
				regret_n = betas:* x_n   
			} 
			else{
				regret_n =RRM_log_sim_data(x_n, betas, mu, gamma)  
			} 
			regret_sum =rowsum(regret_n) 
			R_exp = exp(-regret_sum)
			p_i = 	R_exp :/ colsum(R_exp)
				
			cum_p_i =runningsum(p_i)	
			rand_draws = J(rows(x_n),1,uniform(1,1)) 
			pbb_balance = rand_draws:<cum_p_i
			cum_pbb_balance = runningsum(pbb_balance)
			choice_n = (cum_pbb_balance:== J(rows(x_n),1,1))

		if (n==1)     Y =choice_n	
		else Y=Y \ choice_n	
		}
resindex = st_addvar("byte","choice")
st_store((1,rows(Y)),resindex,Y) 		
}
end


mata:
function RRM_log_sim_data(	real matrix x_n, 
							real rowvector betas,
							real scalar mu,
							real scalar gamma)
{
real scalar i, j ,r_m
real matrix regret_n 
regret_n = J(rows(x_n), cols(x_n), 0)
for(i=1; i <= rows(x_n); ++i) { 
	for(j=1; j <= rows(x_n); ++j) { 
		if (i!=j) { 
			r_i = ln(gamma :+ exp( (betas :/ mu) :* ( x_n[j , . ] :-  x_n[i, . ]))) 				
			regret_n[i, . ] = regret_n[i, . ] :+ mu*r_i 		
			} 
		}	  
	}
return(regret_n)
}
end

//------------------------------//
//-----Simulate classic RRM-----//
//------------------------------//
capture program drop sim_classicRRM
program define sim_classicRRM, rclass	
    version 12
	drop _all
	set obs ${n_obs}
	gen id = _n
	local n_choices =3
	expand `n_choices'
	bys id : gen alternative = _n
	gen x1 =  runiform(-2,2)
	gen x2 =  runiform(-2,2)
	global real_gamma = 1
	global real_mu = 1
	matrix betas_st = (1,-2)
    mata: data_gen()
	randregret choice  x* , altern(alternative)  gr(id) 	///
							rrm(classic) nocons  tech(bfgs) 
end
simulate  	_b _se   ,  reps(${n_rep})  seed(157)  : sim_classicRRM


graph hbox RRM_b_x1 RRM_b_x2, 				///
ytitle(Estimated Parameters)				///
title(Box-Plot Estimated Coefficients) 		///
subtitle(Classic RRM Model) 				///
caption("True Coefficients; \beta = (1,-2)")

graph export "sim_classicRRM_graph.pdf" , as(pdf) replace

sum RRM_b_x1 RRM_b_x2



//------------------------------//
//---------Simulate GRRM--------//
//------------------------------//
capture program drop sim_GRRM
program define sim_GRRM, rclass	
    version 12
	drop _all
	set obs ${n_obs}
	gen id = _n
	local n_choices =3
	expand `n_choices'
	bys id : gen alternative = _n
	gen x1 =  runiform(-2,2)
	gen x2 =  runiform(-2,2)
	global real_gamma = 0.5
	global real_mu = 1
	matrix betas_st = (1,-2)
    mata: data_gen()
	randregret choice  x* , altern(alternative)  gr(id) ///
							rrm(gene) nocons  tech(bfgs) 
end
simulate  	_b _se  				///
			gamma = e(gamma)  		///
			gamma_se = e(gamma_se) 	/// 
								,  reps(${n_rep})  seed(157)  : sim_GRRM

graph hbox RRM_b_x1 RRM_b_x2 _eq5_gamma,	///
ytitle(Estimated Parameters)				///
subtitle(Box-Plot Estimated Coefficients) 		///
title(Generalized RRM Model) 			///
caption("True Coefficients; \beta = (1,-2); gamma= 0.5")

graph export "sim_GRRM_graph.pdf" , as(pdf) replace

sum RRM_b_x1 RRM_b_x2 _eq5_gamma
			
			
//------------------------------//
//--------Simulate muRRM--------//
//------------------------------//
capture program drop sim_muRRM
program define sim_muRRM, rclass	
    version 12
	drop _all
	set obs ${n_obs}
	gen id = _n
	local n_choices =3
	expand `n_choices'
	bys id : gen alternative = _n
	gen x1 =  runiform(-2,2)
	gen x2 =  runiform(-2,2)
	global real_gamma = 1
	global real_mu = 0.75
	matrix betas_st = (1,-2)
    mata: data_gen()
	randregret choice  x* , altern(alternative)  gr(id) ///
							rrm(mu) nocons  tech(bfgs) 
end
simulate  	_b _se  				///
			mu = e(mu)  			///
			mu_se = e(mu_se) 		/// 
								,  reps(${n_rep})  seed(157)  : sim_muRRM

graph hbox RRM_b_x1 RRM_b_x2 _eq5_mu,	///
ytitle(Estimated Parameters)				///
subtitle(Box-Plot Estimated Coefficients) 		///
title(Mu RRM Model) 			///
caption("True Coefficients; \beta = (1,-2); mu= 0.75")

graph export "sim_muRRM_graph.pdf" , as(pdf) replace

sum RRM_b_x1 RRM_b_x2 _eq5_mu
			
			

//------------------------------//
//--------Simulate PRRM --------//
//------------------------------//
capture program drop sim_PRRM
program define sim_PRRM, rclass	
    version 12
	drop _all
	set obs ${n_obs}
	gen id = _n
	local n_choices =3
	expand `n_choices'
	bys id : gen alternative = _n
	gen _x1 =  runiform(-2,2)
	gen _x2 =  runiform(-2,2)

	*Transformation using randregre_pure
	randregret_pure _x1 , gr(id) sign(pos) prefix(x_)
	randregret_pure _x2 , gr(id) sign(neg) prefix(x_)	
	
	matrix betas_st = (1,-2)
    
	global pure_model = "YES"
	mata: data_gen()
	
	randregret choice , altern(alternative)  gr(id) 	///
						rrm(pure) nocons  tech(bfgs) 	///
						pos(_x1) neg(_x2)
						
end

simulate  _b _se 	,  reps(${n_rep})  seed(157)  : sim_PRRM


graph hbox choice_b_x1 choice_b_x2  ,			///
ytitle(Estimated Parameters)					///
subtitle(Box-Plot Estimated Coefficients) 		///
title(Pure RRM Model) 							///
caption("True Coefficients; \beta = (1,-2)")

graph export "sim_PRRM_graph.pdf" , as(pdf) replace

sum choice_b_x1 choice_b_x2



