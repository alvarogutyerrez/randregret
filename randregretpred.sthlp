{smcl}
{* *! version 1.0.0  15may2020}{...}
{title:Title}

{p2colset 7 25 17 2}{...}
{p2col:{bf: randregretpred} {hline 2}}Obtain predictions of systematic regret
and probabilities after
{help randregret:randregret}.{p_end}
{p2colreset}{...}

{p 7 15 2}{cmd:randregretpred} {newvar} {ifin}
{cmd:, {cmdab:alt:ernatives(}varname{cmd:)}}
{cmdab:gr:oup(}varname{cmd:)}  
[{cmd:xb} 
{cmd:proba}]


{synoptset 23 tabbed}{...}
{synopthdr :Options}
{synoptline}
{syntab:Model}
{p2coldent :* {opth gr:oup(varname)}}use {it:varname} to identify cases{p_end}
{p2coldent :* {opth alt:ernatives(varname)}}use {it:varname} to identify the
        alternatives available for each case{p_end}
{synopt:{opt proba:}}calculate probability of a positive outcome; 
		the default{p_end}
{synopt:{opt xb:}}calculate linear prediction of the systematic 
		regret{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:randregretpred} allows to the user to generate predictions of the 
systematic regret and the predicted probabilities for the Random Regret
Minimization models estimated using {help randregret:randregret}.
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
Predicted probabilities from Classical RRM{p_end}
{phang2}{cmd:. randregretpred pred_proba, gr(obs) alt(altern) proba }{p_end}
{pstd}
Predicted systematic regret from Classical RRM{p_end}
{phang2}{cmd:. randregretpred pred_regret, gr(obs) alt(altern) xb}{p_end}
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

{marker cranenburgh2019}{...}
{phang} Van Cranenburgh S. and Alwosheel A. 2019.
{browse "https://www.sciencedirect.com/science/article/pii/S0968090X18305230?via%3Dihub":An artificial neural network based approach to investigate travellers’ decision rules}.
{it:Transportation Research Part C: Emerging Technologies} 98: pp. 152-166.

{marker cranenburgh2018}{...}
{phang} Van Cranenburgh, S. 2018.
{browse "https://doi.org/10.4121/uuid:1ccca375-68ca-4cb6-8fc0-926712f50404":Small value-of-time experiment, Netherlands.} {it:4TU.Centre for Research Data. Dataset}.

