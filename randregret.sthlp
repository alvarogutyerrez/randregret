{smcl}
{* *! version 1.0.0  15may2020}{...}
{title:Title}

{p2colset 1 15 17 2}{...}
{p2col:{bf: randregret} {hline 2}}Classic-, Generalized-, Mu-, and Pure Random 
			Regret Minimization (RRM) Models{p_end}
{p2colreset}{...}

{phang}
Classic Random Regret Minimization (RRM) Model

{p 8 15 2}{cmd:randregret} {depvar} [{indepvars}] {ifin}
[{cmd:, {cmdab:rrm:fn(}classic{cmd:)}} {it:{help randregret##classicoptions:classic_rrm_options}}]

{phang}
Generalized Random Regret Minimization (GRRM) Model

{p 8 16 2}{cmd:randregret} {depvar} [{indepvars}] {ifin}
[{cmd:, {cmdab:rrm:fn(}gene{cmd:)}} {it:{help randregret##geneoptions:generalized_rrm_options}}]

{phang}
Mu Random Regret Minimization (muRRM) Model

{p 8 16 2}{cmd:randregret} {depvar} [{indepvars}] {ifin}
[{cmd:, {cmdab:rrm:fn(}mu{cmd:)}} {it:{help randregret##muoptions:mu_rrm_options}}]

{phang}
Pure Random Regret Minimization (PRRM) Model

{p 8 16 2}{cmd:randregret} {depvar}  {ifin}
[{cmd:, {cmdab:rrm:fn(}pure{cmd:)}} {it:{help randregret##pureoptions:pure_rrm_options}}]

{phang}
{it:depvar} equal to 1 identifies the chosen alternatives,
whereas a 0 indicates the alternatives that were not chosen.
There can be only one chosen alternative for each case.

{marker classicoptions}{...}
{synoptset 23 tabbed}{...}
{synopthdr :classic_rrm_options}
{synoptline}
{syntab:Model}
{p2coldent :* {opth rrm:fn(classic)}}uses the classic systematic regret 
	function model to fit the model.{p_end}
{p2coldent :* {opth gr:oup(varname)}}use {it:varname} to identify cases{p_end}
{p2coldent :* {opth alt:ernatives(varname)}}use {it:varname} to identify the
        alternatives available for each case{p_end}
{synopt:{opt base:alternative(#)}}sets base Alternative Specific 
		Constants (ASC){p_end}
{synopt:{opt nocons:ant}}suppress the alternative specific constants{p_end}

{syntab:SE/Robust}
{synopt :{opth cl:uster(varname)}}implement robust cluster estimator using {it:clustvar} {p_end}
{synopt :{opt r:obust}}implement robust to heteroskedasticity estimator {p_end}
	
{syntab:Reporting}
{synopt :{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}

{syntab:Maximization}
{synopt:{it:maximize_options}}
{opt dif:ficult},
{opt tech:nique(algorithm_spec)}, 
{opt iter:ate(#)}, {opt tr:ace}, {opt grad:ient}, 
{opt showstep}, {opt hess:ian}, {opt tol:erance(#)}, 
{opt ltol:erance(#)} {opt gtol:erance(#)}, {opt nrtol:erance(#)}, 
{opt from(init_specs)}; see {helpb maximize}.
{cmd:technique(bhhh)} is not allowed.
{synoptline}
{p2colreset}{...}

{marker geneoptions}{...}
{synoptset 23 tabbed}{...}
{synopthdr :generalized_rrm_options}
{synoptline}
{syntab:Model}
{p2coldent :* {opth rrm:fn(gene)}}uses the generalized systematic regret 
	function model to fit the model.{p_end}
{p2coldent :* {opth gr:oup(varname)}}use {it:varname} to identify cases{p_end}
{p2coldent :* {opth alt:ernatives(varname)}}use {it:varname} to identify the
        alternatives available for each case{p_end}
{synopt:{opt base:alternative(#)}}sets base Alternative Specific 
		Constants (ASC){p_end}
{synopt:{opt nocons:ant}}suppress the alternative specific constants{p_end}
{synopt:{opt notlr:}}suppress the computations of the Likelihood Ratio (LR)
			test over gamma{p_end}

{syntab:SE/Robust}
{synopt :{opth cl:uster(varname)}}implement robust cluster estimator using {it:clustvar} {p_end}
{synopt :{opt r:obust}}implement robust to heteroskedasticity estimator {p_end}
			
{syntab:Reporting}
{synopt :{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt :{opt show}}show the value of the estimated ancillary 
gamma* parameters.{p_end}

{syntab:Maximization}
{synopt:{it:maximize_options}}{opt dif:ficult},
{opt tech:nique(algorithm_spec)},
{opt iter:ate(#)}, {opt tr:ace}, {opt grad:ient}, 
{opt showstep}, {opt hess:ian},{opt tol:erance(#)}, 
{opt ltol:erance(#)} {opt gtol:erance(#)}, {opt nrtol:erance(#)}, 
{opt from(init_specs)}; see {helpb maximize}.
{cmd:technique(bhhh)} is not allowed.

{synopt:{opt initgamma(#):}}declares the initial value for ancillary parameter 
	gamma*; default value 0.{p_end}
{synoptline}
{p2colreset}{...}

{marker muoptions}{...}
{synoptset 23 tabbed}{...}
{synopthdr :mu_rrm_options}
{synoptline}
{syntab:Model}
{p2coldent :* {opth rrm:fn(mu)}}uses the mu systematic regret 
	function model to fit the model.{p_end}
{p2coldent :* {opth gr:oup(varname)}}use {it:varname} to identify cases{p_end}
{p2coldent :* {opth alt:ernatives(varname)}}use {it:varname} to identify the
        alternatives available for each case{p_end}
{synopt:{opt base:alternative(#)}}sets base Alternative Specific 
		Constants (ASC){p_end}
{synopt:{opt nocons:ant}}suppress the alternative specific constants{p_end}
{synopt:{opt notlr:}}suppress the computations of the Likelihood Ratio (LR)
			test over mu{p_end}
			
{syntab:SE/Robust}
{synopt :{opth cl:uster(varname)}}implement robust cluster estimator using {it:clustvar} {p_end}
{synopt :{opt r:obust}}implement robust to heteroskedasticity estimator {p_end}
			
{syntab:Reporting}
{synopt :{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt :{opt show}}show the value of the estimated ancillary 
mu* parameters.{p_end}

{syntab:Maximization}
{synopt:{it:maximize_options}}{opt dif:ficult},
{opt tech:nique(algorithm_spec)},
{opt iter:ate(#)}, {opt tr:ace}, {opt grad:ient}, 
{opt showstep}, {opt hess:ian},{opt tol:erance(#)}, 
{opt ltol:erance(#)} {opt gtol:erance(#)}, {opt nrtol:erance(#)}, 
{opt from(init_specs)}; see {helpb maximize}.
{cmd:technique(bhhh)} is not allowed.

{synopt:{opt uppermu(#):}}declares the upper bound for the searching space for 
ancillary parameter	mu*; default value 5.{p_end}
{synopt:{opt initmu(#):}}declares the initial value for ancillary parameter 
	mu*; default value 0.{p_end}
{synoptline}
{p2colreset}{...}

{marker pureoptions}{...}
{synoptset 23 tabbed}{...}
{synopthdr :pure_rrm_options}
{synoptline}
{syntab:Model}
{p2coldent :* {opth rrm:fn(pure)}}uses the mu systematic regret 
	function model to fit the model.{p_end}
{p2coldent :* {opth pos:itive(varname)}}use {it:varname} to create transformed
	attributes asuming they have positive sign. {p_end}
{p2coldent :* {opth neg:ative(varname)}}use {it:varname} to create transformed
	attributes asuming they have negative sign. {p_end}	
{p2coldent :* {opth gr:oup(varname)}}use {it:varname} to identify cases{p_end}
{p2coldent :* {opth alt:ernatives(varname)}}use {it:varname} to identify the
        alternatives available for each case{p_end}
{synopt:{opt base:alternative(#)}}sets base Alternative Specific 
		Constants (ASC){p_end}
{synopt:{opt nocons:ant}}suppress the alternative specific constants{p_end}

{syntab:SE/Robust}
{synopt :{opth cl:uster(varname)}}implement robust cluster estimator using {it:clustvar} {p_end}
{synopt :{opt r:obust}}implement robust to heteroskedasticity estimator {p_end}

{syntab:Reporting}
{synopt :{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}

{syntab:Maximization}
{synopt:{it:maximize_options}}{opt dif:ficult},
{opt tech:nique(algorithm_spec)},
{opt iter:ate(#)}, {opt tr:ace}, {opt grad:ient}, 
{opt showstep}, {opt hess:ian},{opt tol:erance(#)}, 
{opt ltol:erance(#)} {opt gtol:erance(#)}, {opt nrtol:erance(#)}, 
{opt from(init_specs)}; see {helpb maximize}.
{cmd:technique(bhhh)} is not allowed.

{synoptline}
{p2colreset}{...}

{marker description}{...}
{title:Description}

{pstd}
{cmd:randregret} fits Random Regret Minimization models allowing for four 
systematic regret functions. In particular, {cmd:randregret} with 
{opt rrmfn(classic)} utilizes the systematic regret function described in 
{help randregret##chorus2010:Chorus. C. (2010)}; 
with {opt rrmfn(gene)} it fits the Generalized RRM model (GRRM) using the 
systematic regret function described in 
{help randregret##chorus2013:Chorus. C. (2013)} 
which allows to have an additional flexibility in the concavity of the function; 
with {opt rrmfn(mu)} it fits the Mu-RRM model, which can estimate the 
scale parameter (mu) of the systematic regret function described in 
{help randregret##cranenburgh2015: Van Cranenburgh S. et al. (2015)}. 
Finally, with  {opt rrmfn(pure)} the model perform the Pure RRM model, which 
invokes {cmd: prrm} to generate transformed regressors and the uses {cmd: clogit}
to fit the model. For details see 
{help randregret##cranenburgh2015: Van Cranenburgh S. et al. (2015)}.

{marker examples}{...}
{title:Examples}

{pstd}
The example comes from data referred to in
{help randregret##cranenburgh2019:Van Cranenburgh.S et at. (2019)}. The data 
was collected as part of a Stated Choice (SC) experiment.
The Choice Situation in this experiment consisted of three unlabeled route
alternatives, each consisting of two generic attributes:
Travel Cost (tc) and Travel Time (tt).
In this experiment each respondent answered a total of 10 choice situation.
The dataset can be downloaded from
{help randregretpred##cranenburgh2018:Van Cranenburgh.S (2018)}.

{pstd}
Open dataset{p_end}
{phang2}{cmd:. scalar server = "https://data.4tu.nl/ndownloader/"}{p_end}
{phang2}{cmd:. scalar doi = "files/24015353"}{p_end}
{phang2}{cmd:. import delimited   "`=server + doi'" ,clear}{p_end}
{phang2}{cmd:. keep obs id cs  tt1 tc1 tt2 tc2 tt3 tc3 choice }{p_end}
{phang2}{cmd:. list obs id cs  tt1 tc1 tt2 tc2 tt3 tc3 choice in 1/4,sepby(obs)}{p_end}
{pstd}
Data preparation{p_end}
{phang2}{cmd:. rename (choice)  (choice_w)}{p_end}
{phang2}{cmd:. reshape long tt tc  , i(obs) j(altern)}{p_end}
{phang2}{cmd:. generate choice = 0}{p_end}
{phang2}{cmd:. replace  choice = 1 if  choice_w==altern  }{p_end}
{phang2}{cmd:. label define  alt_label 1 "First" 2 "Second" 3 "Third" }{p_end}
{phang2}{cmd:. label values  altern alt_label}{p_end}
{phang2}{cmd:. list obs altern choice id cs tt tc   in 1/12, sepby(obs)}{p_end}
{pstd}
Classical RRM ({help randregret##chorus2010:Chorus. C. (2010)}){p_end}
{phang2}{cmd:. randregret choice tc tt, gr(obs) alt(altern) rrmfn(classic) nocons}{p_end}
{pstd}
Generalized RRM ({help randregret##chorus2010:Chorus. C. (2013)}) {p_end}
{phang2}{cmd:. randregret choice tc tt, gr(obs) alt(altern) rrmfn(gene) nocons}{p_end}
{pstd}
muRRM ({help randregret##cranenburgh2015:Van Cranenburgh.S et at. (2015)}){p_end}
{phang2}{cmd:. randregret choice tc tt, gr(obs) alt(altern) rrmfn(mu) nocons }{p_end}
{pstd}
Pure-RRM ({help randregret##cranenburgh2015:Van Cranenburgh.S et at. (2015)}){p_end}
{phang2}{cmd:. randregret choice, neg(tc tt) gr(obs) alt(altern) rrmfn(pure) nocons}{p_end}
{pstd}
Using robust cluster estimation across individuals (id) on Classic RRM {p_end}
{phang2}{cmd:. randregret choice choice tc tt, gr(obs) alt(altern) rrmfn(classic) nocons cluster(id)}{p_end}
{pstd}
Using robust estimation on Classic RRM {p_end}
{phang2}{cmd:. randregret choice choice tc tt, gr(obs) alt(altern) rrmfn(classic) nocons robust}{p_end} 
{pstd}

{marker author}{...}
{title:Authors}

{pstd}
Álvaro A. Gutiérrez Vargas{break}
Faculty of Economics and Business{break}
KU Leuven{break} 
Research Centre for Operations Research and Statistics (ORSTAT){break}
Leuven, Belgium{break}
alvaro.gutierrezvargas@kuleuven.be

{pstd}
Michel Meulders{break}
Faculty of Economics and Business{break}
KU Leuven{break} 
Leuven, Belgium{break}
michel.meulders@kuleuven.be

{pstd}
Martina Vandebroek{break}
Faculty of Economics and Business{break}
KU Leuven{break} 
Leuven, Belgium{break}
martina.vandebroek@kuleuven.be>

{title:References}

{marker chorus2010}{...}
{phang}Chorus. C. 2010.
{browse "https://ojs-lib2.tudelft.nl/ejtir/article/view/2881":A New Model of Random Regret Minimization}.
{it:European Journal of Transport and Infrastructure Research} 10: pp. 181-196.

{marker chorus2013}{...}
{phang}Chorus. C. 2013.
{browse "https://www.sciencedirect.com/science/article/pii/S0191261514001167":A Generalized Random Regret Minimization model}.
{it:Transportation Research Part B: Methodological} 68: pp. 224-238.

{marker cranenburgh2015}{...}
{phang} Van Cranenburgh S., C.A. Guevara and C.G. Chorus 2015.
{browse "https://www.sciencedirect.com/science/article/pii/S0965856415000166":New insights on random regret minimization models}.
{it:Transportation Research Part A: Policy and Practice} 74: pp. 91-109.

{marker cranenburgh2019}{...}
{phang} Van Cranenburgh S. and Alwosheel A. 2019.
{browse "https://www.sciencedirect.com/science/article/pii/S0968090X18305230?via%3Dihub":An artificial neural network based approach to investigate travellers’ decision rules}.
{it:Transportation Research Part C: Emerging Technologies} 98: pp. 152-166.

{marker cranenburgh2018}{...}
{phang} Van Cranenburgh, S. 2018.
{browse "https://doi.org/10.4121/uuid:1ccca375-68ca-4cb6-8fc0-926712f50404":Small value-of-time experiment, Netherlands.} {it:4TU.Centre for Research Data. Dataset}.


