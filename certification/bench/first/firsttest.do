cscript "stub of test script"
display "to be filled in"


clear all
set more off
which randregret


cd data 
use VoT




keep obs id cs tt1 tc1 tt2 tc2 tt3 tc3 choice
list obs id cs tt1 tc1 tt2 tc2 tt3 tc3 choice in 1/4,sepby(obs)
rename (choice) (choice_w)
reshape long tt tc , i(obs) j(altern)
generate choice = 0
replace choice = 1 if choice_w==altern
label define alt_label 1 "First" 2 "Second" 3 "Third"
label values altern alt_label
list obs altern choice id cs tt tc in 1/12, sepby(obs)





//---- Check: Optimization Technique ----//

* tech(hbbb) not allowed
rcof "noisily randregret choice  tt tc ,   altern(altern) gr(obs) rrm(classic)   tech(bhhh) nocons" == 498


//---- Check: -group()- & -alternative()-  ----//

*variable group() needs to be numeric 
gen non_numeric = "a"
rcof "noisily randregret choice  tt tc ,   altern(altern) gr(non_numeric) rrm(classic)   nocons" == 498

*variable alternative() needs to be numeric 
rcof "noisily randregret choice  tt tc ,   altern(non_numeric) gr(obs) rrm(classic)   nocons" == 498

*variable cluster() needs to be numeric 
rcof "noisily randregret choice  tt tc ,   altern(altern) gr(obs) rrm(classic) cl(non_numeric)  nocons" == 498



// ---- Checks on regressors and response variable ----//


// Check that the independent variables vary within groups //
gen x_const = 5
rcof "noisily randregret choice  tt tc x_const ,   altern(altern) gr(obs) rrm(classic)   nocons" == 459				

*The dependent variable must be a 0-1 variable
gen y_constant = 7
rcof "noisily randregret y_constant  tt tc ,   altern(altern) gr(obs) rrm(classic)   nocons" == 450		

* At least one group has more than one chosen alternative
gen y_only_1= 1
rcof "noisily randregret y_only_1  tt tc ,   altern(altern) gr(obs) rrm(classic)   nocons" == 498		



// ---- Check on Alternative Specific Constants (ASC) ----//

*Option basealternative() not compatible with noconstant.
rcof "noisily randregret choice  tt tc ,  base(non_numeric)  altern(altern) gr(obs) rrm(classic)   tech(bfgs) nocons" == 198

*In option basealternative(#), # must be numeric.
rcof "noisily randregret choice  tt tc ,  base(non_numeric)  altern(altern) gr(obs) rrm(classic)   tech(bfgs) " == 198

*Variable in alternatives() does not contain basealternative(#) provided 
rcof "noisily randregret choice  tt tc ,  base(5)  altern(altern) gr(obs) rrm(classic)   tech(bfgs) " == 198


// ---- Checks regret function ----  //

*Non existing regret function
rcof "noisily randregret choice  tt tc ,   altern(altern) gr(obs) rrm(NOT_FOUND)  " == 111



///---Pure: rrmfn(pure)---///
*Explanatory variables must be fill in positive() or negative() options 
rcof "noisily randregret choice  tt tc , altern(altern) gr(obs) rrm(pure)  " == 198

*Explanatory variables are required
rcof "noisily randregret choice  ,  altern(altern) gr(obs) rrm(pure)  " == 100

*notlr not allowed
rcof "noisily randregret choice, pos(tt)  altern(altern) gr(obs) rrm(pure) notlr " == 198

*notlr not allowed
rcof "noisily randregret choice, pos(tt) notlr altern(altern) gr(obs) rrm(pure) show " == 198



///---Classic: rrmfn(classic)---///
*Options positive() and negative() only posible using rrmfn('pure')
rcof "noisily randregret choice   , pos(tt ) neg(tc)  altern(altern) gr(obs) rrm(classic)  " == 198

*Option notlr not allowed
rcof "noisily randregret choice tt tc  ,  altern(altern) gr(obs) rrm(classic) notlr " == 198

*Option show not allowed
rcof "noisily randregret choice tt tc  ,  altern(altern) gr(obs) rrm(classic) show " == 198

*Init value for gamma_star where there is no gamma_star
rcof "noisily randregret choice  tc tt ,   altern(altern) gr(obs) rrm(classic) iter(0) initgamma(-1) " == 198

*Init value for mu_star where there is no mu_star
rcof "noisily randregret choice  tc tt ,   altern(altern) gr(obs) rrm(classic) iter(0) initmu(-1) " == 198

*Init values regressors
matrix init = J(1,4,1)
rcof "noisily randregret choice  tc tt , from(init)  altern(altern) gr(obs) rrm(classic) iter(0)  " == 0




///---Generalized: rrmfn(gamma)---///
*Options positive() and negative() only posible using rrmfn('pure')
rcof "noisily randregret choice   , pos(tt ) neg(tc)  altern(altern) gr(obs) rrm(gene) " == 198

*lr test suppres
rcof "noisily randregret choice  tc tt ,   altern(altern) gr(obs) rrm(gene) notlr iter(0) " == 0

*lr test included
rcof "noisily randregret choice  tc tt ,   altern(altern) gr(obs) rrm(gene) iter(0)  " == 0

*show gamma_star
rcof "noisily randregret choice  tc tt ,   altern(altern) gr(obs) rrm(gene) iter(0) show " == 0

*Init value for gamma_star = -1
rcof "noisily randregret choice  tc tt ,   altern(altern) gr(obs) rrm(gene) iter(0) initgamma(-1) show" == 0

*Init value for mu_star where it should be initgamma
rcof "noisily randregret choice  tc tt ,   altern(altern) gr(obs) rrm(gene) iter(0) initmu(-1) show" == 198

*Init value for all parameters 
matrix init = J(1,5,1)
rcof "noisily randregret choice  tc tt , from(init)  altern(altern) gr(obs) rrm(gene) iter(0)  show" == 0


///---Mu: rrmfn(mu)---///
*Options positive() and negative() only posible using rrmfn('pure')
rcof "noisily randregret choice   , pos(tt ) neg(tc)  altern(altern) gr(obs) rrm(mu)  " == 198

*lr test suppres
rcof "noisily randregret choice  tc tt ,   altern(altern) gr(obs) rrm(mu) notlr iter(0) " == 0

*lr test included
rcof "noisily randregret choice  tc tt ,   altern(altern) gr(obs) rrm(mu) iter(0)  " == 0

*show mu_star
rcof "noisily randregret choice  tc tt ,   altern(altern) gr(obs) rrm(mu) iter(0) show " == 0

*Init value for mu_star = -1
rcof "noisily randregret choice  tc tt ,   altern(altern) gr(obs) rrm(mu) iter(0) initmu(-1) show" == 0

*Init value for gamma_star where it should be initmu
rcof "noisily randregret choice  tc tt ,   altern(altern) gr(obs) rrm(mu) iter(0) initgamma(-1) show" == 198

*Init value for all parameters 
matrix init = J(1,5,1)
rcof "noisily randregret choice  tc tt , from(init)  altern(altern) gr(obs) rrm(mu) iter(0)  show" == 0


///----------------------------------------------------////
///        Tests from version 1.1.0 onwards            ////
///----------------------------------------------------////

//  Check that dependent variable is not part of the regressors 
rcof "noisily randregret choice choice  tc tt , from(init)  altern(altern) gr(obs) rrm(gene) iter(0)  show" == 498
rcof "noisily randregret choice choice  tc tt , from(init)  altern(altern) gr(obs) rrm(classic) iter(0)  show" == 498
rcof "noisily randregret choice choice  tc tt , from(init)  altern(altern) gr(obs) rrm(mu) iter(0)  show" == 498


// The variable choice is specified to have both positive and negative coefficients.
rcof "noisily randregret choice ,  positive( tc) negative( tc)  altern(altern) gr(obs) rrm(pure) iter(0)  " == 498
//  The dependent variable has been used as explanatory variable inside positive() or negative()."
rcof "noisily randregret choice ,  positive( choice ) negative( tc)  altern(altern) gr(obs) rrm(pure) iter(0)  " == 498
//  The dependent variable has been used as explanatory variable inside positive() or negative()."
rcof "noisily randregret choice ,  positive( tc ) negative( choice)  altern(altern) gr(obs) rrm(pure) iter(0)  " == 498

// Omitted regressors when clustering/robust standard errors. 
gen tt_collinear = tt

// Classic / Generalized / Mu (with and without constants.)
// robust 
rcof "noisily randregret choice  tc tt tt_collinear, r altern(altern) gr(obs) rrm(classic) iter(0)  " == 0
rcof "noisily randregret choice  tc tt tt_collinear, r altern(altern) gr(obs) rrm(gene)    iter(0)  show" == 0
rcof "noisily randregret choice  tc tt tt_collinear, r altern(altern) gr(obs) rrm(mu)      iter(0)  show" == 0

rcof "noisily randregret choice  tc tt tt_collinear, r altern(altern) gr(obs) rrm(classic) iter(0)   nocons    " == 0
rcof "noisily randregret choice  tc tt tt_collinear, r altern(altern) gr(obs) rrm(gene)    iter(0)  show nocons" == 0
rcof "noisily randregret choice  tc tt tt_collinear, r altern(altern) gr(obs) rrm(mu)      iter(0)  show nocons" == 0


// cluster
rcof "noisily randregret choice  tc tt tt_collinear, cl(obs) altern(altern) gr(obs) rrm(classic)    " == 0
rcof "noisily randregret choice  tc tt tt_collinear ,cl(obs) altern(altern) gr(obs) rrm(gene)   show" == 0
rcof "noisily randregret choice  tc tt tt_collinear, cl(obs) altern(altern) gr(obs) rrm(mu)     show" == 0

rcof "noisily randregret choice  tc tt tt_collinear, cl(obs) altern(altern) gr(obs) rrm(classic)     nocons" == 0
rcof "noisily randregret choice  tc tt tt_collinear ,cl(obs) altern(altern) gr(obs) rrm(gene)   show nocons " == 0
rcof "noisily randregret choice  tc tt tt_collinear, cl(obs) altern(altern) gr(obs) rrm(mu)     show nocons " == 0


// Pure-RRM.
*Preventing the user from typos. 
//  The variable tc is specified twice in positive() option.
rcof "noisily randregret choice ,  positive( tc tc )  altern(altern) gr(obs) rrm(pure) iter(0)  " == 498
//  The variable tc is specified twice in negative() option.
rcof "noisily randregret choice ,  negative( tc tc )  altern(altern) gr(obs) rrm(pure) iter(0)  " == 498

//  Ready to omitte collinear regressors
rcof "noisily randregret choice ,  negative( tt tt_collinear  )  altern(altern) gr(obs) rrm(pure) iter(0)  " == 0
//  Ready to omitte collinear regressors
rcof "noisily randregret choice ,  positive( tt tt_collinear  )  altern(altern) gr(obs) rrm(pure) iter(0)  " == 0
















