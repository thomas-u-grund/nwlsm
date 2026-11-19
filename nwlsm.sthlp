{smcl}
{* *! version 0.3.0  19nov2026 author: Thomas Grund}{...}
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
{cmd:estat} {c -(}{cmd:positions}|{cmd:mcmcdiag}|{cmd:predict}{c )-} [{cmd:,} {it:options}]

{p 8 16 2}
{cmd:nwlsm version}

{pstd}
{it:netname} is a one-mode network of the nwcommands family, undirected or directed and binary. Self-loops are ignored; ties with a missing value are left out of the likelihood.

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt dim(#)}}number of latent dimensions, 1 to 10; default {cmd:2}{p_end}
{synopt :{opt edgecov(list)}}edge covariates: networks of the same actors, or Stata matrices with one row and column per actor{p_end}
{synopt :{opt absdiff(varlist)}}numeric node variables (one observation per actor, in the order of the actors) entered as the absolute difference |x_i - x_j|{p_end}
{synopt :{opt sender}}random sender effect for every actor (directed networks){p_end}
{synopt :{opt receiver}}random receiver effect for every actor (directed networks){p_end}
{synopt :{opt soc:iality}}random sociality effect for every actor, added to both ends of a tie (undirected and directed networks){p_end}
{synopt :{opt sendervar(#)}, {opt senderdf(#)}}scale (default 1) and degrees of freedom (default 3) of the prior of the sender variance; likewise {opt receivervar()}, {opt receiverdf()}, {opt socialityvar()}, {opt socialitydf()}{p_end}
{synopt :{opt bil:inear}}bilinear (inner-product) latent term instead of the distance{p_end}
{synopt :{opt nsample(#)}}number of posterior draws kept; default {cmd:4000}{p_end}
{synopt :{opt burnin(#)}}sweeps discarded after the tuning runs; default {cmd:10000}{p_end}
{synopt :{opt interval(#)}}sweeps between two kept draws; default {cmd:10}{p_end}
{synopt :{opt pilot(#)}}number of tuning runs for the step sizes; default {cmd:4}{p_end}
{synopt :{opt plen(#)}}sweeps in each tuning run; default {cmd:500}{p_end}
{synopt :{opt mklstarts(#)}}number of starting points of the minimum Kullback-Leibler fit; default {cmd:5}{p_end}
{synopt :{opt seed(#)}}random-number seed (default: drawn from Stata's generator, so {cmd:set seed} makes a fit reproducible){p_end}
{synopt :{opt randoms:tart}}random start positions instead of classical scaling{p_end}
{synopt :{opt savedraws(filename)}}save the draws of the coefficients, the log likelihood and the variances to a Stata data set{p_end}
{synoptline}

{title:Description}

{pstd}
{cmd:nwlsm} fits the latent space model of Hoff, Raftery and Handcock (2002) in the form of the {cmd:euclidean} and {cmd:bilinear} terms of the R package
latentnet (Krivitsky and Handcock): every actor {it:i} has a position z_i in a {it:d}-dimensional space, and the log odds of a tie from {it:i} to {it:j} are

{p 12 12 2}
logit P(y_ij = 1) = b_0 + b_1 x_ij1 + ... + b_k x_ijk + s_i + r_j + o_i + o_j - |z_i - z_j|

{pstd}
where |.| is the Euclidean distance and s_i, r_j and o_i are the optional random sender, receiver and sociality effects. With {opt bilinear} the last term is
replaced by the inner product z_i'z_j, which is added to the log odds. Actors that are close in the latent space (or point in the same direction in the bilinear model)
are more likely to be tied; ties are independent given the positions. For an undirected network every pair counts once. The posterior is sampled by Markov chain Monte
Carlo in a compiled (C) plugin: a random-walk Metropolis-Hastings update of one position at a time, of one random effect at a time and of each coefficient in turn,
and Gibbs updates of the variances of the positions and of the random effects.

{title:Model and priors}

{pstd}
The defaults follow latentnet. The coefficients have independent normal priors with mean 0 and variance 9; for a covariate the variance is 9
divided by the mean square of the covariate over the dyads. The positions are z_i ~ N(0, s2 I), and s2 has a scaled inverse chi-squared prior with
scale 0.125 n^(2/d) and sqrt(n) degrees of freedom, where {it:n} is the number of actors. A random effect has a normal distribution with mean 0 and its own
variance; that variance has a scaled inverse chi-squared prior with scale 1 and 3 degrees of freedom unless {opt sendervar()} and the options like it say otherwise.
The prior on the positions pulls toward a common spread, so the intercept and the spread of the positions are only weakly identified in small or sparse networks;
posterior means of the intercept can then differ from the values used to simulate a network, exactly as in latentnet. A sender, a receiver and a sociality effect
together with the intercept are not separately identified by the likelihood; the proper priors of the random effects identify them.

{pstd}
The step sizes of the random-walk proposals are tuned in {opt pilot()} runs of {opt plen()} sweeps toward an acceptance rate of 0.234 and are then held fixed
during burn-in and sampling. One sweep updates every position once, then every random effect, then every coefficient, then the variances. The start positions come from classical
multidimensional scaling of the shortest-path distances between the actors (disconnected pairs get the largest distance plus one) and the start coefficients
from a logit regression given those distances, unless {opt randomstart} is given. The random effects start at 0. The sample of the posterior is kept: every {opt interval()}th sweep after the burn-in.

{title:Results}

{pstd}
The table shows the posterior mean, the posterior standard deviation, the Monte Carlo standard error of the mean (batch means) and the 95 percent
posterior interval (2.5th and 97.5th percentile) of every coefficient, followed by the same summary for the variance of the positions and of the random effects.
{cmd:e(b)} holds the posterior means and {cmd:e(V)} the sample covariance of the draws.

{pstd}
Positions are identified only up to rotation, reflection and translation (in the bilinear model only up to rotation and reflection). {cmd:e(Z)} holds the minimum
Kullback-Leibler (MKL) positions of latentnet (Handcock, Raftery and Tantrum 2007): the positions, coefficients and random effects that maximise the likelihood of the
posterior mean tie probabilities, so that the fitted tie probabilities are as close as possible, in the Kullback-Leibler sense, to the posterior mean probabilities.
The optimisation starts from the {opt mklstarts()} draws whose tie probabilities are closest to the posterior mean probabilities and keeps the best result; the
MKL fit can have several local optima, and {cmd:mklstarts(1)} starts from the closest draw only, as latentnet does. The fit carries a very weak penalty on the positions and
the random effects (one thousandth of the posterior precision) that keeps the optimiser from running away along directions in which the likelihood is flat; latentnet has no such
penalty. The effects of the MKL fit are centred and their means moved into the intercept {cmd:e(mkl_b)}; {cmd:e(mkl_effects)} holds them. {cmd:e(Zmean)} is the mean of the draws, each aligned
to the MKL positions by Procrustes analysis (rotation and reflection, no scaling); it is shrunken toward the centre. {cmd:e(D)} is the matrix of posterior mean distances between the actors (for
{opt bilinear}: the posterior mean of minus the inner products), which does not depend on the rotation, and {cmd:e(P)} the matrix of posterior mean tie probabilities.
These matrices are stored when the number of actors fits into a Stata matrix.

{synoptset 28 tabbed}{...}
{p2col 5 28 32 2: Stored results}{p_end}
{p2colreset}{...}
{synopt :{cmd:e(N)}}number of dyads (pairs for an undirected network){p_end}
{synopt :{cmd:e(N_nodes)}}number of actors{p_end}
{synopt :{cmd:e(dim)}, {cmd:e(nsample)}, {cmd:e(burnin)}, {cmd:e(interval)}, {cmd:e(seed)}, {cmd:e(mklstarts)}}settings of the run{p_end}
{synopt :{cmd:e(ll_mean)}, {cmd:e(ll_sd)}}posterior mean and standard deviation of the log likelihood{p_end}
{synopt :{cmd:e(sigma2_mean)}, {cmd:e(sigma2_sd)}}posterior mean and standard deviation of s2{p_end}
{synopt :{cmd:e(prior_scale)}, {cmd:e(prior_df)}}scale and degrees of freedom of the prior of s2{p_end}
{synopt :{cmd:e(mkl_start)}, {cmd:e(mkl_iter)}, {cmd:e(mkl_obj)}}draw that started the MKL fit, its iterations, the maximised log likelihood of the posterior mean tie probabilities{p_end}
{synopt :{cmd:e(b)}, {cmd:e(V)}}posterior means and covariance of the coefficients{p_end}
{synopt :{cmd:e(summary)}}mean, sd, Monte Carlo s.e. and 95 percent interval of each coefficient{p_end}
{synopt :{cmd:e(vsummary)}}the same for sigma2 and the variances of the random effects{p_end}
{synopt :{cmd:e(tuning)}}acceptance rates (sampling phase) and final step sizes: latent positions, each coefficient, sender, receiver, sociality{p_end}
{synopt :{cmd:e(Z)}, {cmd:e(Zmean)}}MKL positions; mean of the aligned draws{p_end}
{synopt :{cmd:e(D)}, {cmd:e(P)}}posterior mean distances; posterior mean tie probabilities{p_end}
{synopt :{cmd:e(sender)}, {cmd:e(receiver)}, {cmd:e(sociality)}}posterior mean and sd of each effect (one row per actor){p_end}
{synopt :{cmd:e(mkl_b)}, {cmd:e(mkl_effects)}}coefficients and effects of the MKL fit{p_end}
{synopt :{cmd:e(effects)}, {cmd:e(latent)}}the random effects in the model; {cmd:euclidean} or {cmd:bilinear}{p_end}
{p2colreset}{...}

{pstd}
The full posterior sample stays in Mata after the fit: {cmd:mata: LsmSample("beta")} (draws by coefficient), {cmd:"llk"}, {cmd:"s2"}, {cmd:"ev"}
(variances of the sender, receiver and sociality effects), {cmd:"E"} (the effects by draw), {cmd:"Z"}
(one row per draw: the positions by actor, then dimension), {cmd:"Zalign"} (the aligned mean positions) and {cmd:"mkl"} (the MKL positions).

{title:Post-estimation}

{phang}
{cmd:estat positions} [{cmd:, mean generate(}{it:stub}{cmd:) nolist}] lists the MKL positions with the posterior mean effects ({opt mean}: the mean of the aligned draws instead) and returns them in {cmd:r(Z)}.
With {opt generate()} the coordinates go into the variables {it:stub}1, {it:stub}2, ... of a data set with one observation per actor, in the order of the actors.

{phang}
{cmd:estat mcmcdiag} [{cmd:, nolist}] shows, for the coefficients, the log likelihood and the variance parameters, the posterior mean and standard deviation, the effective sample size
(initial monotone positive sequence estimator of Geyer) and the Geweke z-score (first 10 percent of the draws against the last 50 percent) and returns them in {cmd:r(diag)}. It needs the draws of the last fit, which stay in memory.

{phang}
{cmd:estat predict} [{cmd:, generate(}{it:netname}{cmd:) replace}] returns the posterior mean tie probabilities in {cmd:r(P)} and, with {opt generate()}, stores them as a valued network of the family.

{pstd}
{help nwlsmplot} plots the positions with the ties as light lines, the actors coloured by a variable and sized by an effect or a variable.

{title:Examples}

{phang2}{cmd:. nwlsm friends, dim(2)}{p_end}
{phang2}{cmd:. nwlsm friends, dim(2) absdiff(age) nsample(2000) burnin(5000) seed(12)}{p_end}
{phang2}{cmd:. nwlsm advice, dim(2) sender receiver}{p_end}
{phang2}{cmd:. nwlsm friends, dim(2) sociality}{p_end}
{phang2}{cmd:. estat mcmcdiag}{p_end}
{phang2}{cmd:. estat positions, generate(pos)}{p_end}
{phang2}{cmd:. estat predict, generate(pfriends)}{p_end}
{phang2}{cmd:. nwlsmplot, color(gender) size(sociality)}{p_end}

{title:Reference}

{p 4 8 2}Hoff, P. D., Raftery, A. E. and Handcock, M. S. 2002. Latent space approaches to social network analysis. Journal of the American Statistical Association 97: 1090-1098.{p_end}
{p 4 8 2}Handcock, M. S., Raftery, A. E. and Tantrum, J. M. 2007. Model-based clustering for social networks. Journal of the Royal Statistical Society A 170: 301-354.{p_end}
{p 4 8 2}Krivitsky, P. N. and Handcock, M. S. 2008. Fitting position latent cluster models for social networks with latentnet. Journal of Statistical Software 24(5).{p_end}
{p 4 8 2}Krivitsky, P. N., Handcock, M. S., Raftery, A. E. and Hoff, P. D. 2009. Representing degree distributions, clustering, and homophily in social networks with latent cluster random effects models. Social Networks 31: 204-213.{p_end}

{hline}
{pstd}{it:nwlsm version 0.3.0, 19 November 2026}{p_end}{* nwversion}
