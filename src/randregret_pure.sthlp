{smcl}
{* *! version 1.0.0  15may2020}{...}
{title:Title}

{p2colset 7 28 17 2}{...}
{p2col:{bf: randregret_pure} {hline 2}}Create transformed alternative-specific attributes 
for Pure Random Regret Minization Model.{p_end}
{p2colreset}{...}

{p 7 15 2}{cmd:randregret_pure} {newvar} {ifin}
{cmdab:gr:oup(}varname{cmd:)}  
{cmdab:sign:beta(}varname{cmd:)}  
{cmdab:prefix(}stubname{cmd:)}  


{synoptset 23 tabbed}{...}
{synopthdr :Options}
{synoptline}
{syntab:Model}
{p2coldent :* {opth gr:oup(varname)}}use {it:varname} to identify cases{p_end}
{p2coldent :* {opth sign:beta(string)}}define the expected sign of the attribute.
{cmd:signbeta(pos)} must be used for attributes with positive expected sign, while
{cmd:signbeta(neg)} for attributes with negative expected sign{p_end}
{p2coldent :* {opth prefix(stubname)}}create indicator variables for {it:stubname} {p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:randregret_pure} creates the transformed regressors requiered by the Pure Random 
Regret Minimization (PRRM) model implemented by {help randregret:randregret}. After the 
creation the transformed regressors {help randregret:randregret} invokes
{help clogit:clogit} to implement the {help randregret##pureoptions:PRRM model}.

{marker examples}{...}
{title:Examples}

{pstd}
The example comes from data referred to in
{help randregret##cranenburgh2019:Van Cranenburgh. S. et at. (2019)}. The data 
was collected as part of a Stated Choice (SC) experiment
{help randregret##cranenburgh2018:Van Cranenburgh. S . (2018)}
The Choice Situation in this experiment consisted of three unlabeled route
alternatives, each consisting of two generic attributes:
Travel Cost (tc) and Travel Time (tt).
In the example we will generate the transformed regressors for attributes tc 
and tt, as explained in
{help randregret##cranenburgh2015:Van Cranenburgh. S. et at. (2015)}
under the assumption that both attributes have negative signs, 
given that cheaper and faster routes are preferred to costlier and slower ones. 


{pstd}
Open dataset{p_end}
{phang2}{cmd:. scalar server = "https://data.4tu.nl/ndownloader/"}{p_end}
{phang2}{cmd:. scalar doi = "files/24015353"}{p_end}
{phang2}{cmd:. import delimited   "`=server + doi'" ,clear}{p_end}
{phang2}{cmd:. keep obs id cs  tt1 tc1 tt2 tc2 tt3 tc3 choice }{p_end}
{pstd}
Data preparation{p_end}
{phang2}{cmd:. rename (choice)  (choice_w)}{p_end}
{phang2}{cmd:. reshape long tt tc  , i(obs) j(altern)}{p_end}
{phang2}{cmd:. generate choice = 0}{p_end}
{phang2}{cmd:. replace  choice = 1 if  choice_w==altern  }{p_end}
{phang2}{cmd:. label define  alt_label 1 "First" 2 "Second" 3 "Third" }{p_end}
{phang2}{cmd:. label values  altern alt_label}{p_end}
{pstd}
Generating transformed regressors for Total Cost and Total Time.{p_end}
{phang2}{cmd:. randregret_pure choice tc tt , sign("neg")  group(obs)  prefix(p_)}{p_end}
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

