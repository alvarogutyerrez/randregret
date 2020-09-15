clear all
// Loading data
scalar server = "https://surfdrive.surf.nl/files/"
scalar doi = "index.php/s/pyTXgpXx0JPTHhv/download"
import delimited   "`=server + doi'" ,clear

// Data procesing
rename (choice)  (choice_wide)
gen id_ = _n
reshape long fsg fso tt , i(id_) j(alternative)
gen choice_long =0
replace choice_long =1 if  choice_wide==alternative  
*Follow along Sanders' numerical manipulation
replace fsg  = fsg/1000
replace fso  = fso/1000
replace tt  = tt/100

/*Classic RRM*/ 
randregret choice_long  fsg fso tt ,gr(id_) altern(alternative) rrm(classic) nocons 

/*muRRM*/ 
randregret choice_long  fsg fso tt ,gr(id_) altern(alternative) rrm(mu) show nocons tech(bfgs)

/*Generalized RRM*/ 
randregret choice_long  fsg fso tt ,gr(id_) altern(alternative) rrm(gene) show nocons

/*Pure RRM*/ 
randregret choice_long    , pos(fsg fso) neg(tt)  gr(id_) altern(alternative) rrm(pure)  nocons 
