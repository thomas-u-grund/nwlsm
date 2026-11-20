{smcl}
{* *! version 1.0.2  20nov2026 author: Thomas Grund}{...}
{title:Title}

{p2colset 9 22 29 2}{...}
{p2col :nwlsmselect  {hline 2}}Choose the number of latent clusters of an {help nwlsm} model by BIC{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:nwlsmselect} {it:netname} [{cmd:,} {opt clu:sters(numlist)} {opt eff:obs(string)} {opt ref:it} {opt nol:ist} {opt graph} {opt gname(name)} {opt gopts(string)} {it:nwlsm_options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt clu:sters(numlist)}}numbers of clusters to compare: 0 (no clusters) and 2 to 6; default {cmd:0 2 3 4}{p_end}
{synopt :{opt eff:obs(string)}}effective number of observations of the BIC: {cmd:ties} (default), {cmd:dyads} or {cmd:actors}{p_end}
{synopt :{opt ref:it}}fit the model with the smallest BIC again, so that it is the active model afterwards{p_end}
{synopt :{opt nol:ist}}do not list the table{p_end}
{synopt :{opt graph}}plot the BIC (and the BIC of the ties, of the positions and, if present, of the random effects) against the number of clusters; a dashed line marks the smallest BIC{p_end}
{synopt :{opt gname(name)}}name of the graph{p_end}
{synopt :{opt gopts(string)}}further {help twoway} options for the graph, such as {cmd:title()} or {cmd:scheme()}{p_end}
{synopt :{it:nwlsm_options}}options of {help nwlsm} for every fit, such as {opt dim()}, {opt edgecov()}, {opt absdiff()}, {opt sender}, {opt nsample()} and {opt seed()}; not {opt influence()}, {opt bilinear} or {opt family()}{p_end}
{synoptline}

{title:Description}

{pstd}
{cmd:nwlsmselect} fits {help nwlsm} to the binary network {it:netname} once for every number of clusters in {opt clusters()}, computes the BIC of latentnet's {cmd:bic.ergmm()} for each fit ({cmd:estat bic} of {help nwlsm}) and lists the fits ranked by BIC, with the parts of the BIC for the ties given the positions and for the positions,
the posterior mean log likelihood and the difference to the smallest BIC. All fits use the same options and the same {opt seed()} (without it, a seed is drawn for each fit from Stata's generator). Without {opt refit} the active estimation results are those of the last fit.
The BIC is a rough guide: latentnet itself warns that it is not certain to be appropriate for choosing the dimension, the random effects or between clustered and unclustered models.

{title:Stored results}

{synoptset 16 tabbed}{...}
{synopt :{cmd:r(best)}}number of clusters with the smallest BIC{p_end}
{synopt :{cmd:r(table)}}matrix with one row per fit: clusters, BIC, BIC of the ties, BIC of the positions, BIC of the random effects, mean log likelihood, clusters used by the EM fit, difference to the smallest BIC{p_end}

{title:Example}

{phang2}{cmd:. nwwebuse lazega, nwclear}{p_end}
{phang2}{cmd:. nwlsmselect lazega_fr, clusters(0 2 3 4) absdiff(seniority) seed(12) refit graph}{p_end}

{title:Reference}

{p 4 8 2}Krivitsky, P. N. and Handcock, M. S. 2008. Fitting position latent cluster models for social networks with latentnet. Journal of Statistical Software 24(5).{p_end}

{hline}
{pstd}{it:nwlsm version 1.0.2, 20 November 2026}{p_end}{* nwversion}
