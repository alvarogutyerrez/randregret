*! randregret 1.0.0 26Oct2020
*! author aagv

program  randregret_pure, eclass 
        version 12
	
        syntax varlist(min=1) [if] , 	///
		GRoup(varname) 			 		///
		SIGNbeta(string)				///
		PREFIX(string)
		 
 
		marksample touse
		markout `touse' `group'	

		// check possible errors //	
	
		if ("`group'" != "") {
		capture confirm numeric var `group'
			if _rc != 0 {
				di in r "The group() variable must be numeric"
				exit 498
			}
		}

		*Cheking signbeta option
		if ("`signbeta'" != "pos") & ("`signbeta'" != "neg") {
			di as r "option signbeta() misspecified. Only 'pos' and 'neg' options are allowed"
			}
		else {
			mata: mata_signbeta = st_local("signbeta") 
		}
		
		foreach k in `varlist' {
			qui gen `prefix'`k' = . if `touse'			
			local pure_list_var "`pure_list_var' `prefix'`k'"
		}
		
	mata: st_view(X = ., ., "`varlist'")  
	mata: st_view(panvar = ., ., "`group'") 
	mata: PuRe  = r_pure_per_altern(X , panvar ) 
	mata: PuRe = PuRe *-1

	mata: new_data=.
	mata: st_view(new_data, ., "`pure_list_var'")
	local num : list sizeof local(varlist)
	mata: new_data[.,.] = PuRe[.,.]
	

end

mata:
real matrix r_pure(real matrix x_n )
{
	
	external mata_signbeta
	signbeta = mata_signbeta


	real scalar i, j ,r_m
	real matrix regret_n 
	real matrix zeros 
	zeros = J(1, cols(x_n), 0)
	regret_n = J(rows(x_n), cols(x_n), 0)
	for(i=1; i <= rows(x_n); ++i) { 
		for(j=1; j <= rows(x_n); ++j) { 
			if (i!=j) { 
				dif = (x_n[j,.] :-  x_n[i, . ])
				dif_zero = dif \ zeros
				if 		(signbeta=="pos"){
					r_m = colmax(dif_zero)	
				}
				else if (signbeta=="neg"){
					r_m = colmin(dif_zero)	
				}
				regret_n[i, . ] = regret_n[i, . ] :+ r_m				
				} 
			}
			
		}
return(regret_n)
}

end

mata:
real matrix r_pure_per_altern(real matrix X, real colvector panvar)
{
	paninfo = panelsetup(panvar, 1) 	
	npanels = panelstats(paninfo)[1] 
	for(n=1; n <= npanels; ++n) { 
		x_n = panelsubmatrix(X, n, paninfo)
		regret_n =r_pure(x_n) 
		if (n==1) PuRe_vars =regret_n			  
		else PuRe_vars=PuRe_vars \ regret_n			
	}	
return(PuRe_vars)	
}

end

exit
