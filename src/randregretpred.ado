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
			
			// (2) Dealing with ASC (if specified)	
			if "`e(ASC)'"=="YES" {
				tempvar ASC
				qui tab `alternatives' ,gen(`ASC')
				local n_altern = r(r) 
				if "`e(basealternative)'"!=""{
					drop `ASC'`e(basealternative)'
						}
				else{ //drop the alternative with lower number
					qui sum  `alternatives' , meanonly 
					local min_alt = r(min)
					drop `ASC'`min_alt'
				}
				// Dummy ASC
				mata: st_view(ASC = ., ., "`ASC'*")
				// Coefficients ASC
				matrix ASC_hat = e(b)[1,`e(rank)'-(`n_altern'-2)..`e(rank)']				
				// Coefficients RRM part
				matrix b_hat = e(b)[1,1..`e(rank)'-`n_altern'+1]
				}
			else{ 
				// Dummy ASC (=0 for conformability)
				mata: ASC = J(1,1,0) 	
				// Coefficients ASC
				matrix ASC_hat =J(1,1,0)
				// Coefficients RRM part
				matrix b_hat = e(b)
			}		
			*restore
			// Predicted Probability Computations
			mata: prediction= pbb_pred(X , ASC, panvar) 

			qui gen	 `varlist' = .	
			mata: empty_view = .
			mata: st_view(empty_view, ., "`varlist'")
			mata: empty_view[.,.] =prediction[.,.]		
			qui replace `varlist' =. if e(sample)  != 1 
			qui replace `varlist' =. if   `touse' !=1 
			qui drop `prefix'*
			if "`e(ASC)'"=="YES" drop `ASC'*
		
		}
		else if ("`e(rrmfn)'"!="pure") {
			// X's
			loc cmdline =  "`e(cmdline)'" 
			gettoken randregret rhs : cmdline ,parse(" ,")  // cmdline
			gettoken  y x_options  :  rhs ,parse(" ,") 	    // dependent name
			gettoken covars options : x_options ,parse(",") // covars names
			
			// Mata allocation of the regressors.
			mata: st_view(X = ., ., "`covars'") 	

		if "`e(ASC)'"=="YES" {
			tempvar ASC
			qui tab `alternatives' ,gen(`ASC')
			local n_altern = r(r) 
			if "`e(basealternative)'"!=""{
				drop `ASC'`e(basealternative)'
			}
			else{ //drop the alternative with lower number
				qui sum  `alternatives' , meanonly 
				local min_alt = r(min)
				qui drop `ASC'`min_alt'
			}
			
		mata: st_view(ASC = ., ., "`ASC'*") 		
		matrix ASC_hat=  e(b)[1, "ASC:"]	
		}
		else{ 
			// Dummy ASC (=0 for conformability)
			mata: ASC = J(1,1,0) 	
			// Coefficients ASC
			matrix ASC_hat =J(1,1,0)
		}		
		// Regret function parameters		
		matrix b_hat = e(b)[1, "RRM:"]
		// Ancillary parameter 
		matrix mu_hat    =J(1,1,1)
		matrix gamma_hat =J(1,1,1)
		if ("`e(rrmfn)'"=="mu") {
			matrix mu_hat =  e(mu)
		}
		else if ("`e(rrmfn)'"=="gene")  {
			matrix gamma_hat = e(gamma)
		}		
		// Predicted Probability Computations
		mata: prediction= pbb_pred(X , ASC, panvar) 

		if "`e(ASC)'"=="YES" qui drop `ASC'*

		qui gen double	`varlist' = .	
		mata: empty_view = .
		mata: st_view(empty_view, ., "`varlist'")
		mata: empty_view[.,.] =prediction[.,.]
		qui replace `varlist' =. if e(sample)  != 1  
		qui replace `varlist' =. if   `touse' !=1 
		}	
		
end