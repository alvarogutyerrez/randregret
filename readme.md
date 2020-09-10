# randregret: A command for fitting random regret minimization models using Stata 

```
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
```

Here we describe the ```randregret``` command, which implements a variety of Random Regret Minimization (RRM) models. The command allows the user to apply the classic RRM model introduced in Chorus (2010, _European Journal of Transport and Infrastructure Research_ 10: 181-196), the Generalized RRM model introduced in Chorus (2014, Transportation Research Part B: Methodological 68: 224-238), and also the $\mu$RRM and Pure RRM models, both introduced in van Cranenburgh (2015, _Transportation Research Part A: Policy and Practice_ 74: 91-109). We illustrate the usage of the ```randregret``` command using stated choice data on route preferences. The command offers robust and cluster standard error correction using analytical expressions of the scores functions. It also offers likelihood ratio tests that can be used to assess the relevance of a given model specification. Finally, users can obtain the predicted probabilities from each model using the ```randregretpred``` command.

```keywords```: randregret, randregret_pure, randregretpred, discrete choice models,  semi-compensatory behavior, random utility maximization, random regret minimization.

