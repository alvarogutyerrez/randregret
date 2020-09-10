*! randregret 1.0.0 26Oct2020
*! author aagv


/***********************************************************
   ___   ____          __    ___   ___  __   ___   ___  ____
  /__/  ____/  /\  /  /  \  /__/  /__  / _  /__/  /__    /
 /  \  /___/  /  \/  /___/ /  \  /__  /__/ /  \  /__    /   

 
 V.1.0:  d0 ml evaluator that run the following RRM models:
		
		->RRM      (Chorus, 2010)
		->muRRM    (S. Van Cranenburgh et.al, 2015)
		->pure-RRM (S. Van Cranenburgh et.al, 2015)
		->G-RRM    (Chorus, 2014)	
		
************************************************************/

program randregret 
	version 12

	if replay() {
	if ("`e(cmd)'" != "randregret") error 301
	Replay `0'
	}
	else Estimate `0'
end


program Estimate, eclass sortpreserve 
	syntax varlist(fv) [if] [in] ,		///
			GRoup(varname) 				///
			RRMfn(string) 				///
			ALTernatives(varname)		///
			[BASEalternative(string)	///
			UPPERMU(integer 5)			///
			POSitive(varlist)			///
			NEGative(varlist)			///
			NOCONStant					///
			INITGAMMA(real 0)     	 	///
			INITMU(real 0)      		///
			NOTLR						///
			SHOWancillary				///
			Robust						///
			CLuster(varname)			///
			noLOg 						/// 
			Level(cilevel)				///
			TRace						///
			GRADient					///
			HESSian						///
			SHOWSTEP					///
			ITERate(passthru)			///
			TOLerance(passthru)			///
			LTOLerance(passthru)		///
			GTOLerance(passthru)		///
			NRTOLerance(passthru)		///
			CONSTraints(passthru)		///
			TECHnique(passthru)			///
			DIFficult					///
			FRom(string)				/// 
	]

	
	
	
	local mlopts `trace' `gradient' `hessian' `showstep' `iterate' `tolerance' ///
	`ltolerance' `gtolerance' `nrtolerance' `constraints' `technique' `difficult' `from'
	

	// globals to mata
	global group_mata =  "`group'"
	global alternatives_mata =  "`alternatives'"
	global cluster_mata =  "`cluster'"
	
	
	// Previous checks.
	if ("`technique'" == "technique(bhhh)") {
	di in red "technique(bhhh) is not allowed."
	exit 498
	}

	capture confirm numeric var `group'
	if _rc != 0 {
		di in r "The group variable must be numeric"
		exit 498
	}

	capture confirm numeric var `alternatives'
	if _rc != 0 {
		di in r "The alternative variable must be numeric"
		exit 498
	}
	if ("`cluster'" != "") {
		capture confirm numeric var `cluster'
		if _rc != 0 {
			di in r "The cluster variable must be numeric"
			exit 498
		}
		else{
			tempvar last
			bys `group': gen `last' = cond(_n==_N,1,0)
			mata: cluster_var_mata = st_data(., st_global("cluster_mata"), "`last'")
		}	
	}


	// check syntax
	gettoken lhs rhs : varlist

	// mark the estimation sample
	marksample touse
	markout `touse' `group' `alternatives' `cluster'

	preserve
	qui keep if `touse'
	// Check that the independent variables vary within groups //
	sort `group'
	foreach var of varlist `rhs' `positive' `negative'  {
		capture by `group': assert `var'==`var'[1]
		if (_rc == 0) {
			di in red "Variable `var' has no within-group variance"
			exit 459		
		}
	}	
	// Check that the dependent variable only takes values 0-1 //
	capture assert `lhs' == 0 | `lhs' == 1
	if (_rc != 0) {
		di in red "The dependent variable must be a 0-1 variable indicating which alternatives are chosen"
		exit 450		
	}
	// Check that each group has only one chosen alternative //
	tempvar chonum
	sort `group'
	qui by `group': egen `chonum' = sum(`lhs')
	capture assert `chonum' == 1
	if (_rc != 0) {
		di in red "At least one group has more than one chosen alternative"
		exit 498		
	}
	restore
	
	
	
	// Type of Regret Function //
	if  ("`rrmfn'" =="mu"	) 		|  	///
		("`rrmfn'" =="classic") 	|	///
		("`rrmfn'" =="gene") 		|	///
		("`rrmfn'" =="pure") 		{
		global regret_fn "`rrmfn'"
		}
	else {	
		di as red "Regret function misspecified: rrmfn()  "
		exit 111
	}	
	
	// ASC Checks //
	if "`noconstant'" ==""{ 
		global cons_demanded = "YES"
		global neq = 2 /*ml display purposes*/
	}
	else {
		global cons_demanded = "NO" 	
		global neq = 1 /*ml display purposes*/
	}
	if "`showancillary'" != "" global showancillary  ""	
	else global showancillary  "neq(${neq})"	

	
	// Checking if ASC are demanded and well specified
	// for Classic RRM, Generalized  RRM and muRRM 
	if "${cons_demanded}" =="YES"{  	
		// basealternative(#) , # must be numeric. 
		capture confirm number `basealternative'
		if _rc == 7 & "`basealternative'"!="" { // basealternative not numeric
			di in red "In option basealternative(#), # must be numeric." 
			exit 198
			}
		// basealternative(#) , # must be within values on alternatives().	
		else {
			qui levelsof `alternatives', miss local(mylevs)
			loc check_base_altern = strpos("`mylevs'", "`basealternative'")
			if `check_base_altern' == 0 {
				di in red "Variable in alternatives() does not contain basealternative(#) provided " 
				exit 198
				}
			}
		// Generate ASC as tempvars
		qui levelsof `alternatives', local(levels_altern)
		// tailored made ASC variables
		tempvar ASC_
		foreach i of local levels_altern {
			tempvar ASC_`i'
			qui gen `ASC_'`i' = (`alternatives' == `i')
		}				
		if "`basealternative'" ==""{ // if basealternative not specified use min 
			qui sum  `alternatives' , meanonly 
			local min_alt = r(min)
			drop `ASC_'`min_alt' // drop ASC with min value from `alternatives'
		} 
		else{
			capture drop `ASC_'`basealternative'
		}
		// Generate local with ASC equation for ML //
		local ASC_vars = "(ASC: `ASC_'*   , nocons)"
	}
	else { // if no constant are specified,  ASC_vars is empty
		if "`basealternative'" !=""{
			di in red "Option basealternative() not compatible with noconstant."
			exit 198
			}
		else {
		    local ASC_vars = ""
		}		
	}

	
	// ----- Different Regret Functions ------ //
if "`rrmfn'" == "pure" {
	
	
	// Not possible to include explanatory variables in rhs 
	if "`rhs'" !=""  {
		di in red "Explanatory variables must be fill in positive() or negative() options when using rrmfn('pure')"
		exit 198
	}
	// Need to include at least a block of variables in pos() or neg()
	if ("`positive'" =="") &  ("`negative'" =="")  {
		di in red "Explanatory variables must be fill in positive() or negative() options when using rrmfn('pure')"
		exit 198
	}
	// show not allowed
	if ("`showancillary'" !="")   {
		di in red "showancillary is not compatible with rrmfn(pure)"
		exit 198
	}
	// noLR: not possible
	if ("`notlr'" != "")  {
		di in red "nolrt is not compatible with rrmfn(pure)"
		exit 198
	}
	
	preserve 
	qui keep if `touse'==1

	tempvar prefix _temp

	capture  randregret_pure `positive' if `touse', gr(`group')  sign(pos) prefix(`prefix') 
	capture  randregret_pure `negative' if `touse', gr(`group')  sign(neg) prefix(`prefix')

	loc old_covars `positive' `negative'
	foreach k in `old_covars' {
		local pure_covars "`pure_covars' `prefix'`k'"
		local _temp_old_vars "`_temp_old_vars' `_temp'`k'" 
	}
		
	// old variables now in a temp name
	rename ( `old_covars' )  ( `_temp_old_vars' )	
	// pure_vars into the name of old vars. Just for the sake of presentation
	rename ( `pure_covars' )  ( `old_covars' )

	// cluster errors
	if "`cluster'" !="" loc vce_cluster  "vce(cluster `cluster')"
	
	if "${cons_demanded}" =="YES"{
		foreach k of varlist `ASC_'* {
			// Reverse sign of ASC
			qui replace `k' = -`k' 
			}		
		qui clogit `lhs' `old_covars' `ASC_'*   if `touse' ,group(`group') ///
											from(`from') `robust' `vce_cluster'
	}
	else{
		qui clogit `lhs' `old_covars' 		  if `touse' ,group(`group') ///
											from(`from') `robust' `vce_cluster'
	}
		
	restore		
	
	if "${cons_demanded}" =="YES"{
		tempname b_pure
		matrix `b_pure' = e(b)
		// replace tempvar names of ASC for meaningfull names
		local names : colnames `b_pure'
		qui tokenize `names'
		foreach i of var `ASC_'*{
			*position ASC_i
			loc pos_col_i = colnumb(`b_pure', "`i'")
			*id ASC_i
			loc id_ASC_i =  substr("`i'",-1,.) 	
			*replacement
			local `pos_col_i' ASC_`id_ASC_i'
			local newnames "`*'"
			*assign corrected names of ASC.
			matrix colnames `b_pure' = `newnames'
		}
		*Adjust ereturn names
		ereturn repost b=`b_pure' ,  rename esample(`touse')
	}
	else{
		// recovering esample
		tempname b_pure
		matrix `b_pure' = e(b)
		ereturn repost b=`b_pure' ,  rename esample(`touse')
	}
	// number of distinct choice situations 
	mata: st_view(panvar = ., ., st_global("group_mata"))
	mata: paninfo = panelsetup(panvar, 1) 	
	mata: npanels = panelstats(paninfo)[1]
	mata: st_numscalar("__n_cases", npanels)
	ereturn scalar n_cases = __n_cases

	// wald statistic over attribute-specific variables
	qui test `old_covars' 
	ereturn scalar chi2 = r(chi2)
	ereturn scalar p =  r(p)
	ereturn scalar df_m = r(df)
	
	
	ereturn local title "PRRM: Pure Random Regret Minimization Model"	
	ereturn local  rrmfn  "`rrmfn'"
	ereturn local positive "`positive'"
	ereturn local negative "`negative'"	
	ereturn local basealternative  "`basealternative'"
	ereturn local ASC  "${cons_demanded}"
	ereturn local cmd "randregret"
	
	
	
	Header
	randregret 
	di `"The Pure-RRM uses a transformation of the original regressors using options"' 		 
	di `"positive() and negative() as detailed in {browse "https://www.sciencedirect.com/science/article/pii/S0965856415000166":S. van Cranenburgh et. al (2015)}"'
	di "Afterward, randregret invokes clogit using these transormed regresors"
	di as text "{hline 78}"
	exit		

	}

		
	
else if "`rrmfn'" == "classic" {	
	// Previous checks //	
	if ("`positive'" !="") |  ("`negative'" !="")  {
		di in red "Options positive() and negative() only posible using rrmfn('pure')"
		exit 198
	}
	if ("`notlr'" != "")   {
		di in red "Option notlr not allowed; only possible with rrmfn(mu) or rrmfn(gene)"
		exit 198
	}
	if ("`showancillary'" !="")   {
		di in red "Option show not allowed; only possible with rrmfn(mu) or rrmfn(gene)"
		exit 198
	}
	if (`initmu' !=0)  {
		di in red "Options initmu() only posible using rrmfn('mu')"
		exit 198
	}
	if (`initgamma' !=0)  {
		di in red "Options initgamma() only posible using rrmfn('gene')"
		exit 198
	}
		
	
	local RR_log = "(RRM: `lhs' = `rhs', nocons)"
	local LL  `RR_log' `ASC_vars'
	local title "RRM: Classic Random Regret Minimization Model"		
	local fitting_message = "Fitting Classic RRM Model " 
	}


else if "`rrmfn'" =="mu" {
	// Previous checks
	if ("`positive'" !="") |  ("`negative'" !="")  {
		di in red "Options positive() and negative() only posible using rrmfn('pure')"
		exit 198
	}
	if (`initgamma' !=0)  {
		di in red "Options initgamma() only posible using rrmfn('gene')"
		exit 198
	}
		
	// Ancilliary param mu upperbond 
	global uppermu = `uppermu'
	// ML equations
	local RR_mu = "(mu_star:   )"		
	local RR_log = "(RRM: `lhs' = `rhs', nocons) " 
	local LL  `RR_log' `ASC_vars' `RR_mu'
	local title "muRRM: Mu-Random Regret Minimization Model"
	
	// Param display
	local diparm diparm(mu_star, /*
		*/ function(invlogit(@)*`uppermu') /*
		*/ derivative( invlogit(@)*(1- invlogit(@))*`uppermu' ) /*
		*/  label("mu"))
	
	local fitting_message = "Fitting muRRM Model" 
	}				

else if  "`rrmfn'" == "gene" {
	// Previous checks
	if ("`positive'" !="") |  ("`negative'" !="")  {
		di in red "Options positive() and negative() only posible using rrmfn('pure')"
		exit 198
	}

	if (`initmu' !=0)  {
		di in red "Options initmu() only posible using rrmfn('mu')"
		exit 198
	}
		
	
	
	// ML equations		
	local RR_log = "(RRM: `lhs' = `rhs', nocons)"
	local RR_gamma = "(gamma_star:   )"
	local LL  `RR_log' `ASC_vars' `RR_gamma'
	local title "GRRM: Generalized Random Regret Minimization Model"
	
	// Param display
	local diparm diparm(gamma_star,  invlogit   label("gamma" ) ) 
	local fitting_message = "Fitting Generalized RRM Model " 
	}		

	
	//  --------------------------------------------------  //	
	//  -----------------  ML fitting  -------------------  //
	//  --------------------------------------------------  //	
	
	// initial values	
	if "`from'"!=""{
		if ("`rrmfn'" == "classic") | ///
		   ("`rrmfn'" == "gene") 	| ///
		   ("`rrmfn'" == "mu") 	   {
			local init "init(`from', copy skip )"
		}
	}
	else { // Fitting Classic RRM for initial values
		if ("`rrmfn'" == "gene") | ("`rrmfn'" == "mu") {
			global regret_fn_user = "`rrmfn'"
			global regret_fn "classic"
			di as text "{hline 78}"	
			di "Fitting Classic RRM for Initial Values"
			di as text "{hline 78}"	

			// Fitting classic RRM model for initial values
			ml model d0 randregret_LL() `RR_log' `ASC_vars' if `touse', ///
				`vce' missing first `log' search(on) 					/// 
				maximize 	

			tempname classic_lr 			
			estimate store `classic_lr' 
			
			if ("`rrmfn'" == "gene") matrix RRM_beta=e(b),`initgamma'
			else if ("`rrmfn'" == "mu") matrix RRM_beta=e(b),`initmu'
			

			loc b_classic RRM_beta		
			local init "init(`b_classic', copy skip )"
			// back to the original demanded model in rrmfn()
			global regret_fn = "${regret_fn_user}"	
		}
	}
		
		
	// LR test for Generalized and muRRM
	global notlr = "`notlr'"
	if ("`notlr'" == "") & (("`rrmfn'" == "gene") | ("`rrmfn'" == "mu")) {
		tempname mu_lr_1 					///
				 mu_unrestricted			///
				 gene_lr_0 					///
				 gene_lr_1 					///
				 gene_unrestricted 
		
		if ("`rrmfn'" == "gene") {
			// First run restricted model with gamma= 0
			// Chorus (2013) proved the equivalence between the loglikelihood of
			// RUM model and GRRM model with gamma= 0 therefore the restricted
			// model is a RUM model.
			di as text "{hline 78}"	
			di "Fitting Conditional Logit as a Restricted Model (gamma=0) for LR test"
			di as text "{hline 78}"	
			
			if "${cons_demanded}" =="YES" qui clogit `lhs'  `rhs' `ASC_'* if `touse' ,group(`group') 
			else qui clogit `lhs'  `rhs'  if `touse'  ,group(`group')
			estimate store `gene_lr_0' 
		}
	}

	// Fitting demanded model
	di as text "{hline 78}"	
	di "`fitting_message'"
	di as text "{hline 78}"	
	ml model d0 randregret_LL()						///
				`LL'								///
				if `touse', 						///
				`vce' 								///
				missing								///
				first 								///	
				`log'								///
				`mlopts' 							///
				`diparm'  							///
				search(on)	 						///
				repeat(10)	 						///				
				`init'								/// 
				title(`title') 						///
				maximize

	// Utilities for robust/cluster standard errors.
	tempname b_all ASC_hat b_hat mu_star_hat gamma_star_hat D
	matrix `b_all' = e(b)
	
	capture mat `ASC_hat' = e(b)[1,"ASC:"]
	mata: ASC_hat = st_matrix("`ASC_hat'")
	
	
	// replace tempvar names of ASC for meaningfull names
	if "${cons_demanded}" =="YES"{ 
		local names : colnames `b_all'
		qui tokenize `names'
		foreach i of var `ASC_'*{
			*position ASC_i
			loc pos_col_i = colnumb(`b_all', "`i'")
			*id ASC_i
			loc id_ASC_i =  substr("`i'",-1,.) 	
			*replacement
			local `pos_col_i' ASC_`id_ASC_i'
			local newnames "`*'"
			*assign corrected names of ASC.
			matrix colnames `b_all' = `newnames'
		}
	}
	
	
	// Utilities for cluster/robust
	capture mat `b_hat' = e(b)[1,"RRM:"]
	mata: b_hat = st_matrix("`b_hat'")
	
	capture mat `mu_star_hat'= e(b)[1,"mu_star:"]
	mata: mu_star_hat = st_matrix("`mu_star_hat'")
	
	capture mat `gamma_star_hat'= e(b)[1,"gamma_star:"]
	mata: gamma_star_hat = st_matrix("`gamma_star_hat'")
	
	matrix `D' = e(V)	
	mata: D = st_matrix("`D'")
	mata: st_view(Y= ., ., "`lhs'", "`touse'" )
	mata: st_view(X = ., ., "`rhs'", "`touse'") 			
	capture mata: st_view(ASC = ., ., "`ASC_'*", "`touse'") 		
	mata: st_view(panvar=., ., "`group'", "`touse'")


	if ("`cluster'" != "") | ("`robust'" !=""){
		tempname  V_r	
		mata: st_matrix("`V_r'", std_errs(D,Y,X,ASC,panvar))
		ereturn repost b=`b_all' V=`V_r',  rename
		ereturn local vcetype Robust
		if ("`cluster'" != "") ereturn scalar n_clusters = __n_clusters
	}
	else{
		ereturn repost b=`b_all' ,  rename
	}
	
	tempname fitted_model
	qui estimate store `fitted_model'
	// saving coef and sd error of ancilliary params in original scale. 
	if  "`rrmfn'" == "mu" {
		// Display Params back into original scale
		qui _diparm mu_star, function(invlogit(@)*`uppermu') /*
			*/derivative(invlogit(@)*(1- invlogit(@))*`uppermu')  label("mu")
		ereturn scalar mu = r(est)
		ereturn scalar mu_se = r(se)
		
		// LR test if demanded
		if ("`notlr'" == ""){
			*H_0: mu=1
			qui lrtest  `fitted_model' `classic_lr'  ,force
			global chi_mu_1 = r(chi2) 
			ereturn scalar  chi_mu_1 = $chi_mu_1
		}	
	}
	else if  "`rrmfn'" == "gene" {
		qui _diparm gamma_star, invlogit label("gamma")
		ereturn scalar gamma = r(est)
		ereturn scalar gamma_se = r(se)
		// LR test if demanded
		if ("`notlr'" == ""){
			*H_0: gamma=1
			qui lrtest  `fitted_model' `classic_lr'  ,force
			global chi_gamma_1 = r(chi2) 
			ereturn scalar lr_gamma_1 = $chi_gamma_1
			*H_0: gamma=0
			qui lrtest  `fitted_model' `gene_lr_0'  ,force
			global chi_gamma_0 = r(chi2) 
			ereturn scalar  lr_gamma_0 = $chi_gamma_0	
		}	
	}
	// correction wald test on display
	qui test `rhs' 
	ereturn scalar chi2 = r(chi2)
	ereturn scalar p =  r(p)
	ereturn scalar df_m = r(df)
	// Display Results
	ereturn local cmd randregret
	// # Cases
	ereturn scalar n_cases = __n_cases
	Header
	Replay , level(`level')

	/*randregretpred utilities */
	ereturn local basealternative  "`basealternative'"
	ereturn local ASC  "${cons_demanded}"
	ereturn local negative "`negative'"	
	ereturn local rrmfn  "`rrmfn'"
	ereturn local cmdline  "randregret `0'"
end

program Header

	loc df_m = e(df_m)
	di ""
	di in gr e(title)
	di ""
	di in gr `"Case ID variable: `=abbrev("${group_mata}",24)'"' _col(48) /// Case ID
	 "Number of cases" _col(67) "= " in ye %10.0g e(n_cases) 

	di in gr `"Alternative variable: `=abbrev("${alternatives_mata}",24)'"' _col(48) /// 
	 "Number of obs" _col(67) "= " in ye %10.0g e(N) 	 
	 
	di in gr _col(48) "Wald chi2(`df_m')" _col(67) "= "  in ye  %10.2f e(chi2) 
 
	di in gr `"Log likelihood = "' in ye %10.0g e(ll)   _col(48) /// 
	in gr _col(48)"Prob > chi2"  _col(67) "= "  in ye %10.4f e(p) 
 
    if (e(n_clusters)!=.){
		di in gr _col(34) "(Std. Err. adjusted for " in ye %5.0f e(n_clusters) in gr`" clusters in `=abbrev("${cluster_mata}",24)')"'  
	}

end 
			
 program Replay
	syntax [, Level(cilevel) ]
	ml display , level(`level')  $showancillary noheader 
	if  ("$regret_fn" == "gene") & ("$notlr" == "") {
		_lrtest_gene ,  msg("LR test of gamma=0") chi_test(${chi_gamma_0})  chibar01upperthan
		_lrtest_gene ,  msg("LR test of gamma=1") chi_test(${chi_gamma_1})  chibar01lowerthan
		di as text "{hline 78}"	
		}
	else if ("$regret_fn" == "mu") & ("$notlr" == "") {
		_lrtest_gene ,  msg("LR test of mu=1") chi_test(${chi_mu_1})  chi2
		di as text "{hline 78}"	
	}

 end

program define _lrtest_gene, sclass
	version 12

	syntax , msg(string) chi_test(real) 	///
			 [chibar01upperthan 			///
			 chibar01lowerthan 				///
			 chi2]                       	///
			 
	tempname pval
	if "`chibar01upperthan'" != "" {
		scalar `pval' =  chiprob(1, `chi_test')*0.5
		if `chi_test'<=1e-7 scalar `pval'= 1
		
		if ((`chi_test' > 0.005) & (`chi_test'<1e4)) | (`chi_test'==0) {
			local fmt "%8.2f"
		}
		else  local fmt "%8.2e"

		local chi : di `fmt' `chi_test'
		local chi = trim("`chi'")

		di "{txt}`msg': " _c
		di "{help j_chibar##|_new:chibar2(01) = }{res}`chi'" _c ///
			_col(56) "{txt}Prob >= chibar2 = {res}" %5.3f `pval' _n
			
	}
	else if "`chibar01lowerthan'" != "" {
	
		scalar `pval' =   chiprob(1, `chi_test')*0.5
		if `chi_test'<=1e-7 scalar `pval'= 1
		
		if ((`chi_test' > 0.005) & (`chi_test'<1e4)) | (`chi_test'==0) {
			local fmt "%8.2f"
		}
		else    local fmt "%8.2e"

		local chi : di `fmt' `chi_test'
		local chi = trim("`chi'")

		di "{txt}`msg': " _c
		di "{help j_chibar##|_new:chibar2(01) = }{res}`chi'" _c ///
			_col(56) "{txt}Prob >= chibar2 = {res}" %5.3f `pval' _n
	}
	else  if "`chi2'" != "" {
		
		scalar `pval' =  chiprob(1, `chi_test')
		if `chi_test'==0 scalar `pval'= 1
		if ((`chi_test' > 0.005) & (`chi_test'<1e4)) | (`chi_test'==0) {
			local fmt "%8.2f"
		}
		else    local fmt "%8.2e"

		local chi : di `fmt' `chi_test'
		local chi = trim("`chi'")

		di "{txt}`msg': " _c
		di "chi2(1) ={res}`chi'" _c ///
			_col(56) "{txt}Prob >= chibar2 = {res}" %5.3f `pval' _n
		}
		
end

// include mata functions from randregret_LL.mata
findfile "randregret_LL.mata"
capture do "`r(fn)'"

exit
