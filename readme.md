## *```randregret```*: A command for fitting random regret minimization models using Stata 

```
*! randregret 1.1.0 24Dec2020
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
```


Here we describe the ```randregret``` command (forthcoming at [The Stata Journal](https://www.stata-journal.com/)), which implements a variety of Random Regret Minimization (RRM) models. The command allows the user to apply the classic RRM model introduced in Chorus (2010), the Generalized RRM model introduced in Chorus (2014), and also the muRRM and Pure RRM models, both introduced in van Cranenburgh (2015). We illustrate the usage of the ```randregret``` command using stated choice data on route preferences using data from van Cranenburgh (2018). The command offers robust and cluster standard error correction using analytical expressions of the scores functions. It also offers likelihood ratio tests that can be used to assess the relevance of a given model specification. Finally, users can obtain the predicted probabilities from each model using the ```randregretpred``` command.

```keywords```: randregret, randregret_pure, randregretpred, discrete choice models,  semi-compensatory behavior, random utility maximization, random regret minimization.


### Install ```randregret``` 

``` 
*Describe randregret
net describe randregret, from("https://raw.githubusercontent.com/alvarogutyerrez/randregret/master/src/")


*Install randregret
cap ado uninstall randregret
net install randregret, from("https://raw.githubusercontent.com/alvarogutyerrez/randregret/master/src/")
```


### Examples 

```
/*----  PREAMBLE  ----*/
clear all
set more off
capture log close 

// randregret instalation 
cap ado uninstall randregret
net install randregret, from("https://raw.githubusercontent.com/alvarogutyerrez/randregret/master/src/")

// Data download
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

// Different Regret Models:
// Classic RRM+ cluster(id)
randregret choice  tc tt , gr(obs) alt(altern) rrmfn(classic)  cluster(id)	nocons

// muRRM + cluster(id)
randregret choice  tc tt , gr(obs) alt(altern) rrmfn(mu) cluster(id) show  	nocons

// Generalized RRM + cluster(id)
randregret choice  tc tt , gr(obs) alt(altern) rrmfn(gene) cluster(id) show nocons 

// Pure RRM + cluster(id)
randregret choice  , neg(tc tt) gr(obs) alt(altern) rrmfn(pure) cluster(id) nocons     
```

### Documentation 

The full documentation/working paper of the command can be found [here](https://www.dropbox.com/s/3ke1x3syjtmrdtg/randregret_WP.pdf?dl=0)!

### Conferences

> *   [Stata Meeting 2020, London, UK](https://events.timberlake.co.uk/event/2020-stata-conference). Slides available [here](http://repec.org/usug2020/Gutierrez_u20.pdf). 

> *   [Stata Meeting 2020, Bern, Switzerland](https://ritme.com/CH-de/ritme/unsere-neuigkeiten/2020-swiss-stata-conference/). Slides are coming soon. 


### Further information 

Professor [Sander van Cranenburgh](https://www.tudelft.nl/en/tpm/about-the-faculty/departments/engineering-systems-and-services/people/assistant-professors/drir-s-sander-van-cranenburgh/) has a [website](https://www.advancedrrmmodels.com/) where you can find additional information about RRM models together with routines in other languages to fit them (R, Python, and Matlab included). A dofile that replicated professor Sander's results can be found [here](https://github.com/alvarogutyerrez/randregret/blob/master/examples/randregret_replication_Sanders_routines.do).  







### References 

> * Chorus. C. 2010.  A New Model of Random Regret Minimization.  European Journal of Transport and Infrastructure Research 10: pp. 181-196.

> *    Chorus. C. 2013.  A Generalized Random Regret Minimization model.  Transportation Research Part B: Methodological 68: pp. 224-238.

> *    Van Cranenburgh S., C.A. Guevara and C.G. Chorus 2015.  New insights on random regret minimization models.  Transportation Research Part A: Policy and Practice 74: pp. 91-109.

> *    Van Cranenburgh S. and Alwosheel A. 2019.  An artificial neural network based approach to investigate travellers’ decision rules.  Transportation Research Part C: Emerging Technologies 98: pp. 152-166.

> *    Van Cranenburgh, S. 2018.  Small value-of-time experiment, Netherlands. 4TU.Centre for Research Data. Dataset.











