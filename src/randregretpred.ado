*! version 1.0.0 01Jun2020
program define randregretpred, eclass
        version 12
	
        syntax newvarname  [if] [in] , 	///
		GRoup(varname) 			 		///
		ALTernatives(varname)			///
		[PROBA							///
		XB]

		marksample touse , novarlist
		markout `touse'  `group' `alternatives'

		
		// Pre-checks
		if "`e(cmd)'"=="" 		exit 301 
		if ("`e(rrmfn)'"=="") {
			di as red "randregretpred must be performed after randregret"
			exit 301
			}	
		if ("`proba'"=="" & "`xb'"=="") {
			loc proba = "proba"
			di "(option proba assumed; probability of success given one success within group)"
		}  	
		
		
		// panvar, regarless type of model.
		mata: st_view(panvar = ., ., "`group'") 
		
		// In case of PRRM:
		if ("`e(rrmfn)'"=="pure") {
			// (1) Dealing with transformed regressors
			tempvar prefix
			tempvar _temp
			loc negative= e(negative) 
			loc positive= e(positive) 		
			
			capture  randregret_pure `positive' , gr(`group')  sign(pos) prefix(`prefix')
			capture  randregret_pure `negative' , gr(`group')  sign(neg) prefix(`prefix')

			mata: st_view(X = ., ., "`prefix'*") 	
			
			tempname ASC_hat b_hat b_all
			matrix `b_all' = e(b)
			// (2) Dealing with ASC (if specified)	
			if "`e(ASC)'"=="YES" {
				// Generate ASC as tempvars
				qui levelsof `alternatives', local(levels_altern)
				// tailored made ASC variables
				tempvar ASC_
				foreach i of local levels_altern {
					tempvar ASC_`i'
					qui gen `ASC_'`i' = (`alternatives' == `i')
				}				
				qui tab `alternatives' 
				local n_altern = r(r) 
				if "`e(basealternative)'"!=""{
					drop `ASC_'`e(basealternative)'
				}
				else{ //drop the alternative with lower number
						qui sum  `alternatives' , meanonly 
						local min_alt = r(min)
						qui drop `ASC_'`min_alt'
				}
				// Dummy ASC
				mata: st_view(ASC = ., ., "`ASC_'*")
				// Coefficients ASC
				matrix `ASC_hat' = `b_all'[1,`e(rank)'-(`n_altern' - 2)..`e(rank)']				
				// Coefficients RRM part
				matrix `b_hat' = `b_all'[1,1..`e(rank)'-`n_altern'+1]
				}
			else{ 
				// Dummy ASC (=0 for conformability)
				mata: ASC = J(1,1,0) 	
				// Coefficients ASC
				matrix `ASC_hat' =J(1,1,0)
				// Coefficients RRM part
				matrix `b_hat' = `b_all'
			}		
			// Mata allocation of estimates
			mata: b_hat= st_matrix("`b_hat'")
			mata: ASC_hat= st_matrix("`ASC_hat'")
			// Predicted Probability Computations
			mata: prediction= pbb_pred(X , ASC, panvar) 
			qui gen	 `varlist' = .	
			mata: empty_view = .
			mata: st_view(empty_view, ., "`varlist'")
			mata: empty_view[.,.] =prediction[.,.]		
			qui replace `varlist' =. if e(sample)  != 1 
			qui replace `varlist' =. if   `touse' !=1 
		
		}
		else if ("`e(rrmfn)'"!="pure") {
			// X's
			loc cmdline =  "`e(cmdline)'" 
			gettoken randregret rhs : cmdline ,parse(" ,")  // cmdline
			gettoken  y x_options  :  rhs ,parse(" ,") 	    // dependent name
			gettoken covars options : x_options ,parse(",") // covars names
			
			// Mata allocation of the regressors.
			mata: st_view(X = ., ., "`covars'") 	
			
			// Parsing parameter vector
			tempname b_all b_hat ASC_hat b_hat aux_star_hat 	
			matrix `b_all' = e(b)
			if "${cons_demanded}" =="NO"{ 
				if ("`e(rrmfn)'" == "classic") {
					matrix `b_hat' = `b_all' 	
				} 
				else if ("`e(rrmfn)'" == "gene") | ("`e(rrmfn)'" == "mu") {
					matrix `b_hat' = `b_all'[1,1..`e(rank)'-1]
				}
			}
			else if "${cons_demanded}" =="YES"{
				qui tab `alternatives' 
				local n_altern = r(r) 
				if ("`e(rrmfn)'"== "classic") {
					matrix `b_hat' 	 = `b_all'[1,1..(`e(rank)'-`n_altern')+1] 		 
					matrix `ASC_hat' = `b_all'[1,`e(rank)'-(`n_altern'-2)..`e(rank)']
				}
				else if ("`e(rrmfn)'" == "gene") | ("`e(rrmfn)'" == "mu") {
					matrix `b_hat' 	 = `b_all'[1,1..(`e(rank)'-`n_altern')] 		 
					matrix `ASC_hat' = `b_all'[1,`e(rank)'-(`n_altern'-1)..`e(rank)'-1]    
				}
			}
			

			// generate ASC if needed
			if "`e(ASC)'"=="YES" {
				// Generate ASC as tempvars
				qui levelsof `alternatives', local(levels_altern)
				// tailored made ASC variables
				tempvar ASC_
				foreach i of local levels_altern {
					tempvar ASC_`i'
					qui gen `ASC_'`i' = (`alternatives' == `i')
				}				
				qui tab `alternatives' 
				local n_altern = r(r) 
				if "`e(basealternative)'"!=""{
					drop `ASC_'`e(basealternative)'
				}
				else{ //drop the alternative with lower number
						qui sum  `alternatives' , meanonly 
						local min_alt = r(min)
						qui drop `ASC_'`min_alt'
				}
				// Mata allocation of the ASC variables
				mata: st_view(ASC = ., ., "`ASC_'*") 		
			}
			else{ 
				// Dummy ASC (=0 for conformability)
				mata: ASC = J(1,1,0) 	
				// Coefficients ASC
				matrix `ASC_hat' =J(1,1,0)
			}				
			// Ancillary parameter 			
			tempname mu_hat gamma_hat
			matrix `mu_hat'    =J(1,1,1)
			matrix `gamma_hat' =J(1,1,1)
			if ("`e(rrmfn)'"=="mu") {
				matrix `mu_hat' =  e(mu)
		}
		else if ("`e(rrmfn)'"=="gene")  {
			matrix `gamma_hat' = e(gamma)
		}
		
		mata: b_hat= st_matrix("`b_hat'")
		mata: ASC_hat= st_matrix("`ASC_hat'")
		mata: gamma_hat= st_matrix("`gamma_hat'")
		mata: mu_hat= st_matrix("`mu_hat'")
		
		// Predicted Probability Computations
		mata: prediction= pbb_pred(X , ASC, panvar) 

		qui gen double	`varlist' = .	
		mata: empty_view = .
		mata: st_view(empty_view, ., "`varlist'")
		mata: empty_view[.,.] =prediction[.,.]
		qui replace `varlist' =. if e(sample)  != 1  
		qui replace `varlist' =. if   `touse' !=1 
		}	
		

		
end

