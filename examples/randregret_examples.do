/*====================================*/
/*====      DATA LOANDING     ========*/
/*====================================*/

/*----  PREAMBLE  ----*/
clear all
set more off
capture log close 

cd "C:\"


// Data loanding
scalar server = "https://data.4tu.nl/ndownloader/"
scalar doi = "files/24015353"
import delimited   "`=server + doi'" ,clear
keep obs id cs  tt1 tc1 tt2 tc2 tt3 tc3 choice 
list obs id cs  tt1 tc1 tt2 tc2 tt3 tc3 choice in 1/4,sepby(obs)

// Data processing
rename (choice)  (choice_w)
reshape long tt tc  , i(obs) j(altern)
generate choice = 0
replace  choice = 1 if  choice_w==altern  
label define  alt_label 1 "First" 2 "Second" 3 "Third" 
label values  altern alt_label
list obs altern choice id cs tt tc   in 1/12, sepby(obs)


// Classic RRM 
randregret choice  tc tt , gr(obs) alt(altern) rrmfn(classic) nocons   


// Classic RRM + cluster(id)
randregret choice  tc tt, gr(obs) alt(altern) rrmfn(classic) ///
nocons cluster(id) nolog 


// Generalized RRM + cluster(id)
randregret choice  tc tt , gr(obs) alt(altern) rrmfn(gene) ///
nocons  cluster(id) show  


// muRRM + cluster(id)
randregret choice tc tt, gr(obs) alt(altern) rrm(mu) ///
nocons  show  cluster(id) 


// Pure RRM + cluster(id)
randregret choice  , neg(tc tt)  gr(obs) alt(altern) rrmfn(pure) ///
nocons  cluster(id)     


// replication PRRM using randregret_pure
randregret_pure tc tt , sign(neg)  gr(obs)  prefix(p_)
list obs altern choice tt p_tt tc p_tc  in 1/3, sepby(obs)
clogit choice p_tc  p_tt, gr(obs) vce(cluster id)
randregret choice tc tt, gr(obs) alt(altern) base(1) rrmfn(classic) nolog   


// Classic RRM + ASC
randregret choice tc tt, cl(id) gr(obs) alt(altern) nocons rrmfn(classic) nolog   

// Classic RRM + randregret_pred for predictions.
qui randregret choice  tc tt , gr(obs) alt(altern) rrmfn(classic) nocons nolog  
randregretpred prob ,gr(obs) alt(altern) proba
randregretpred xb ,gr(obs) alt(altern) xb
list obs altern choice id cs tt tc prob xb  in 1/12, sepby(obs)



