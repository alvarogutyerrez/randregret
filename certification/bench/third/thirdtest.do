cscript "stub of test script"
display "to be filled in"


clear all
set more off
which randregret
which randregret_pure
which randregretpred


ssc install estout
* Simulated data to test robust/cluster error correction 

* Checking the variance equivalence of:
*	(*)	non-robust
*	(*) robust
*	(*) robust cluster
* 						...	on a homoscedastic (large) data.


*RRM Models:

* classic 	= 2 attributes -> b = (1,-2)

* gene 		= 2 attributes -> b = (1,-2) & gamma = 0.5 

* mu   		= 2 attributes -> b = (1,-2) & mu 	 = .75  

* pure 		= 2 (pure)attributes -> b = (1,-2)  

global n_obs = 250000
global tol = 1e-4


mata:
void data_gen() 
{
		st_view(X = ., ., "x*")
		st_view(panvar = ., ., "id") 
		betas  = st_matrix("betas_st")
		gamma  = strtoreal(st_global("real_gamma")) 
		mu 	   = strtoreal(st_global("real_mu")) 
			
		pure_model 	 = st_global("pure_model") 
			
		paninfo = panelsetup(panvar, 1)     
        npanels = panelstats(paninfo)[1] 
		
		for(n=1; n <= npanels; ++n) { 
            x_n = panelsubmatrix(X, n, paninfo) 
			if (pure_model =="YES") {
				regret_n = -1*betas:* x_n 
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
version 12
drop _all
set seed 157
set obs ${n_obs}
gen cluster_id = _n
local n_cluster=3
expand `n_cluster'
sort cluster_id
gen id = _n
local n_choices =3
expand `n_choices'
bys id : gen alternative = _n
sort  cluster_id id

gen x1 =  runiform(-2,2)
gen x2 =  runiform(-2,2)

global real_gamma = 1
global real_mu = 1
matrix betas_st = (1,-2)
mata: data_gen()


// Non-robust variance
eststo  c :randregret choice  x* , altern(alternative)  gr(id) 		///
						rrm(classic) nocons  from(betas_st) 

matrix classicRRM_V = e(V)						
						
// Robust variance
eststo  cR :randregret choice  x* , altern(alternative)  gr(id) 	///
						rrm(classic) nocons  from(betas_st)			///
						robust

matrix classicRRM_V_robust = e(V)
						
// Cluster-robust variance
eststo cC : randregret choice  x* , altern(alternative)  gr(id) 	///
						rrm(classic) nocons  from(betas_st)			///
						cluster(cluster_id)
matrix classicRRM_V_cluster = e(V)
						
esttab  , label nodepvar nonumber    b(a6) se(a6)
esttab  using "var_correction_classicRRM"  , nostar label nodepvar nonumber  b(a6) se(a6) replace
estimates clear						


display 			 mreldif(classicRRM_V, classicRRM_V_robust)
rcof "noisily assert mreldif(classicRRM_V, classicRRM_V_robust) < ${tol}" == 0

display				 mreldif(classicRRM_V, classicRRM_V_cluster)
rcof "noisily assert mreldif(classicRRM_V, classicRRM_V_cluster) <${tol}" == 0







//------------------------------//
//---------Simulate GRRM--------//
//------------------------------//
version 12
drop _all
set seed 157
set obs ${n_obs}
gen cluster_id = _n
local n_cluster=3
expand `n_cluster'
sort cluster_id
gen id = _n
local n_choices =3
expand `n_choices'
bys id : gen alternative = _n
sort  cluster_id id

gen x1 =  runiform(-2,2)
gen x2 =  runiform(-2,2)

global real_gamma = 0.5
global real_mu = 1
matrix betas_st = (1,-2)
mata: data_gen()


matrix init_gene  = betas_st,logit(${real_gamma})

// Non-robust variance
eststo  g :randregret choice  x* , altern(alternative)  gr(id) 			///
						rrm(gene) nocons  from(init_gene)  show notlr  	///

matrix GRRM_V = e(V)						
												
// Robust variance
eststo  gR :randregret choice  x* , altern(alternative)  gr(id) 		///
						rrm(gene) nocons  from(init_gene)  show notlr  	///
						robust

matrix GRRM_V_robust = e(V)						
		
// Cluster-robust variance
eststo gC : randregret choice  x* , altern(alternative)  gr(id) 		///
						rrm(gene) nocons  from(init_gene)  show notlr 	///
						cluster(cluster_id)

matrix GRRM_V_cluster = e(V)						
		
esttab  , label nodepvar nonumber    b(a6) se(a6)
esttab  using "var_correction_GRRM"  , nostar label nodepvar nonumber  b(a6) se(a6) replace
estimates clear								


display				 mreldif(GRRM_V, GRRM_V_robust)
rcof "noisily assert mreldif(GRRM_V, GRRM_V_robust) < ${tol}" == 0

display				 mreldif(GRRM_V, GRRM_V_cluster)
rcof "noisily assert mreldif(GRRM_V, GRRM_V_cluster) < ${tol}" == 0


											
//------------------------------//
//--------Simulate muRRM--------//
//------------------------------//

version 12
drop _all
set seed 157
set obs ${n_obs}
gen cluster_id = _n
local n_cluster=3
expand `n_cluster'
sort cluster_id
gen id = _n
local n_choices =3
expand `n_choices'
bys id : gen alternative = _n
sort  cluster_id id

gen x1 =  runiform(-2,2)
gen x2 =  runiform(-2,2)

global real_gamma = 1
global real_mu = 0.75
matrix betas_st = (1,-2)
mata: data_gen()

*invlogit(-1.72)*5 \approx 0.75 in transformed scale
matrix init_mu = betas_st, -1.72

// Non-robust variance
eststo  u :randregret choice  x* , altern(alternative)  gr(id) 		///
						rrm(mu) nocons notlr show  from(init_mu)  	///

						
matrix muRRM_V = e(V)						
// Robust variance
eststo  uR :randregret choice  x* , altern(alternative)  gr(id) 		///
						rrm(mu) nocons notlr show  from(init_mu) 		///
						robust

matrix muRRM_V_robust= e(V)						
						
// Cluster-robust variance
eststo uC : randregret choice  x* , altern(alternative)  gr(id) 		///
						rrm(mu) nocons notlr show  from(init_mu) 		///
						cluster(cluster_id)

matrix muRRM_V_cluster = e(V)						
		
		
esttab  , label nodepvar nonumber    b(a6) se(a6)
esttab  using "var_correction_muRRM"  , nostar label nodepvar nonumber  b(a6) se(a6) replace
estimates clear		


display				 mreldif(muRRM_V, muRRM_V_robust)
rcof "noisily assert mreldif(muRRM_V, muRRM_V_robust) < ${tol}" == 0

display				 mreldif(muRRM_V, muRRM_V_cluster)
rcof "noisily assert mreldif(muRRM_V, muRRM_V_cluster) < ${tol}" == 0


//------------------------------//
//--------Simulate PRRM --------//
//------------------------------//

version 12

drop _all
set seed 157
set obs ${n_obs}
gen cluster_id = _n
local n_cluster=3
expand `n_cluster'
sort cluster_id
gen id = _n
local n_choices =3
expand `n_choices'
bys id : gen alternative = _n
sort  cluster_id id

gen _x1 =  runiform(-2,2)
gen _x2 =  runiform(-2,2)

*Transformation using randregre_pure
randregret_pure _x1 , gr(id) sign(pos) prefix(x_)
randregret_pure _x2 , gr(id) sign(neg) prefix(x_)	

global pure_model = "YES"

global real_gamma 	= 1
global real_mu 		= 1
matrix betas_st 	= (1,-2)
mata: data_gen()


// Non-robust variance
eststo p : randregret choice , altern(alternative)  gr(id) 		///
						rrm(pure) nocons  tech(bfgs) 			///
						pos(_x1) neg(_x2)						

matrix PRRM_V = e(V)						
						
// Robust variance
eststo pR : randregret choice , altern(alternative)  gr(id) 	///
						rrm(pure) nocons  tech(bfgs) 			///
						pos(_x1) neg(_x2)						///
						robust
matrix PRRM_V_robust = e(V)						
// Cluster-robust variance
eststo pC : randregret choice , altern(alternative)  gr(id) 	///
						rrm(pure) nocons  tech(bfgs) 			///
						pos(_x1) neg(_x2)						///
						cluster(cluster_id)

matrix PRRM_V_cluster = e(V)						
								
esttab  , label nodepvar nonumber    b(a6) se(a6)							
esttab  using "var_correction_PRRM"  , nostar label nodepvar nonumber  b(a6) se(a6) replace
estimates clear		

display				 mreldif(PRRM_V, PRRM_V_robust)
rcof "noisily assert mreldif(PRRM_V, PRRM_V_robust) < ${tol}" == 0

display				 mreldif(PRRM_V, PRRM_V_cluster)
rcof "noisily assert mreldif(PRRM_V, PRRM_V_cluster) < ${tol}" == 0

