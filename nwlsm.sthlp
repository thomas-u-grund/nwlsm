{smcl}
{* *! version 0.2.0  19nov2026 author: Thomas Grund}{...}
{title:Title}

{p2colset 9 17 24 2}{...}
{p2col :nwlsm  {hline 2}}Latent space model for the ties of a binary network, fitted by MCMC{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:nwlsm} {it:netname} [{cmd:,} {it:options}]

{p 8 16 2}
{cmd:nwlsm}  (replay after a fit)

{p 8 16 2}
{cmd:nwlsm version}

{pstd}
{it:netname} is a one-mode network of the nwcommands family, undirected or directed and binary. Self-loops are ignored; ties with a missing value are left out of the likelihood.

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt dim(#)}}number of latent dimensions, 1 to 10; default {cmd:2}{p_end}
{synopt :{opt edgecov(list)}}edge covariates: networks of the same actors, or Stata matrices with one row and column per actor{p_end}
{synopt :{opt absdiff(varlist)}}numeric node variables (one observation per actor, in the order of the actors) entered as the absolute difference |x_i - x_j|{p_end}
{synopt :{opt nsample(#)}}number of posterior draws kept; default {cmd:4000}{p_end}
{synopt :{opt burnin(#)}}sweeps discarded after the tuning runs; default {cmd:10000}{p_end}
{synopt :{opt interval(#)}}sweeps between two kept draws; default {cmd:10}{p_end}
{synopt :{opt pilot(#)}}number of tuning runs for the step sizes; default {cmd:4}{p_end}
{synopt :{opt plen(#)}}sweeps in each tuning run; default {cmd:500}{p_end}
{synopt :{opt seed(#)}}random-number seed (default: drawn from Stata's generator, so {cmd:set seed} makes a fit reproducible){p_end}
{synopt :{opt randoms:tart}}random start positions instead of classical scaling{p_end}
{synopt :{opt savedraws(filename)}}save the draws of the coefficients, the log likelihood and the position variance to a Stata data set{p_end}
{synoptline}

{title:Description}

{pstd}
{cmd:nwlsm} fits the latent space model of Hoff, Raftery and Handcock (2002) in the form of the {cmd:euclidean} term of the R package
latentnet (Krivitsky and Handcock): every actor {it:i} has a position z_i in a {it:d}-dimensional space, and the log odds of a tie from {it:i} to {it:j} are

{p 12 12 2}
logit P(y_ij = 1) = b_0 + b_1 x_ij1 + ... + b_k x_ijk - |z_i - z_j|

{pstd}
where |.| is the Euclidean distance. Actors that are close in the latent space are more likely to be tied; ties are independent given the positions.
For an undirected network every pair counts once. The posterior is sampled by Markov chain Monte Carlo in a compiled (C) plugin: a random-walk
Metropolis-Hastings update of one position at a time, of each coefficient in turn, and a Gibbs update of the variance of the positions.

{title:Model and priors}

{pstd}
The defaults follow latentnet. The coefficients have independent normal priors with mean 0 and variance 9; for a covariate the variance is 9
divided by the mean square of the covariate over the dyads. The positions are z_i ~ N(0, s2 I), and s2 has a scaled inverse chi-squared prior with
scale 0.125 n^(2/d) and sqrt(n) degrees of freedom, where {it:n} is the number of actors. The prior on the positions pulls toward a common spread, so
the intercept and the spread of the positions are only weakly identified in small or sparse networks; posterior means of the intercept can then differ
from the values used to simulate a network, exactly as in latentnet.

{pstd}
The step sizes of the random-walk proposals are tuned in {opt pilot()} runs of {opt plen()} sweeps toward an acceptance rate of 0.234 and are then held fixed
during burn-in and sampling. One sweep updates every position once, then every coefficient, then the variance. The start positions come from classical
multidimensional scaling of the shortest-path distances between the actors (disconnected pairs get the largest distance plus one) and the start coefficients
from a logit regression given those distances, unless {opt randomstart} is given. The sample of the posterior is kept: every {opt interval()}th sweep after the burn-in.

{title:Results}

{pstd}
The table shows the posterior mean, the posterior standard deviation, the Monte Carlo standard error of the mean (batch means) and the 95 percent
posterior interval (2.5th and 97.5th percentile) of every coefficient. {cmd:e(b)} holds the posterior means and {cmd:e(V)} the sample covariance of the draws.

{pstd}
Positions are identified only up to rotation, reflection and translation. {cmd:e(D)} is the matrix of posterior mean distances between the actors, which does not
depend on the rotation. {cmd:e(Z)} holds the mean positions: every draw is centred and rotated onto a reference by Procrustes analysis (rotation and reflection, no scaling),
the reference being the draw with the highest log posterior, refined twice by the mean of the aligned draws. The minimum Kullback-Leibler positions of latentnet are
not implemented. {cmd:e(D)} and {cmd:e(Z)} are stored when the number of actors fits into a Stata matrix.

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Stored results}{p_end}
{p2colreset}{...}
{synopt :{cmd:e(N)}}number of dyads (pairs for an undirected network){p_end}
{synopt :{cmd:e(N_nodes)}}number of actors{p_end}
{synopt :{cmd:e(dim)}, {cmd:e(nsample)}, {cmd:e(burnin)}, {cmd:e(interval)}, {cmd:e(seed)}}settings of the run{p_end}
{synopt :{cmd:e(ll_mean)}, {cmd:e(ll_sd)}}posterior mean and standard deviation of the log likelihood{p_end}
{synopt :{cmd:e(sigma2_mean)}, {cmd:e(sigma2_sd)}}posterior mean and standard deviation of s2{p_end}
{synopt :{cmd:e(prior_scale)}, {cmd:e(prior_df)}}scale and degrees of freedom of the prior of s2{p_end}
{synopt :{cmd:e(b)}, {cmd:e(V)}}posterior means and covariance of the coefficients{p_end}
{synopt :{cmd:e(summary)}}mean, sd, Monte Carlo s.e. and 95 percent interval of each coefficient{p_end}
{synopt :{cmd:e(tuning)}}acceptance rates (sampling phase) and final step sizes: latent positions, then each coefficient{p_end}
{synopt :{cmd:e(D)}, {cmd:e(Z)}}posterior mean distances; Procrustes-aligned mean positions{p_end}
{p2colreset}{...}

{pstd}
The full posterior sample stays in Mata after the fit: {cmd:mata: LsmSample("beta")} (draws by coefficient), {cmd:"llk"}, {cmd:"s2"}, {cmd:"Z"}
(one row per draw: the positions by actor, then dimension) and {cmd:"Zalign"} (the aligned mean positions).

{title:Examples}

{phang2}{cmd:. nwlsm friends, dim(2)}{p_end}
{phang2}{cmd:. nwlsm friends, dim(2) absdiff(age) nsample(2000) burnin(5000) seed(12)}{p_end}
{phang2}{cmd:. matrix list e(D)}{p_end}

{title:Reference}

{p 4 8 2}Hoff, P. D., Raftery, A. E. and Handcock, M. S. 2002. Latent space approaches to social network analysis. Journal of the American Statistical Association 97: 1090-1098.{p_end}
{p 4 8 2}Krivitsky, P. N. and Handcock, M. S. 2008. Fitting position latent cluster models for social networks with latentnet. Journal of Statistical Software 24(5).{p_end}

{hline}
{pstd}{it:nwlsm version 0.2.0, 19 November 2026}{p_end}{* nwversion}
