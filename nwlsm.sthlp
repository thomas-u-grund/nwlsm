{smcl}
{* *! version 0.6.0  19nov2026 author: Thomas Grund}{...}
{title:Title}

{p2colset 9 17 24 2}{...}
{p2col :nwlsm  {hline 2}}Latent space model for the ties of a binary or valued network, and for social influence, fitted by MCMC{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:nwlsm} {it:netname} [{cmd:,} {it:options}]

{p 8 16 2}
{cmd:nwlsm}  (replay after a fit)

{p 8 16 2}
{cmd:estat} {c -(}{cmd:positions}|{cmd:mcmcdiag}|{cmd:clusters}|{cmd:predict}|{cmd:gof}|{cmd:ppcheck}|{cmd:influence}{c )-} [{cmd:,} {it:options}]

{p 8 16 2}
{cmd:nwlsm version}

{pstd}
{it:netname} is a one-mode network of the nwcommands family, undirected or directed, binary or (with {opt family()}) valued. Self-loops are ignored; ties with a missing value are left out of the likelihood.

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt fam:ily(name)}}distribution of the tie values: {cmd:bernoulli} (default; binary ties), {cmd:poisson}, {cmd:binomial} or {cmd:normal} (alias {cmd:gaussian}){p_end}
{synopt :{opt tri:als(spec)}}binomial family: the number of trials, as a number, a network or a Stata matrix{p_end}
{synopt :{opt dispv:ar(#)}, {opt dispd:f(#)}}normal family: scale (default 1) and degrees of freedom (default 3) of the prior of the dispersion (the variance of the values){p_end}
{synopt :{opt dim(#)}}number of latent dimensions, 1 to 10; default {cmd:2}{p_end}
{synopt :{opt edgecov(list)}}edge covariates: networks of the same actors, or Stata matrices with one row and column per actor{p_end}
{synopt :{opt absdiff(varlist)}}numeric node variables (one observation per actor, in the order of the actors) entered as the absolute difference |x_i - x_j|{p_end}
{synopt :{opt sender}}random sender effect for every actor (directed networks){p_end}
{synopt :{opt receiver}}random receiver effect for every actor (directed networks){p_end}
{synopt :{opt soc:iality}}random sociality effect for every actor, added to both ends of a tie (undirected and directed networks){p_end}
{synopt :{opt sendervar(#)}, {opt senderdf(#)}}scale (default 1) and degrees of freedom (default 3) of the prior of the sender variance; likewise {opt receivervar()}, {opt receiverdf()}, {opt socialityvar()}, {opt socialitydf()}{p_end}
{synopt :{opt bil:inear}}bilinear (inner-product) latent term instead of the distance{p_end}
{synopt :{opt clus:ters(#)}}number of latent clusters, 2 to 6 (latentnet's {cmd:G}); default none{p_end}
{synopt :{opt inf:luence(varlist)}}latent space model for social influence: the outcome, the lagged outcome and any other regressors (variables with one observation per actor){p_end}
{synopt :{opt nei:ghbors(#)}}influence model: number of nearest actors in the latent space whose lagged outcomes act on an actor; default {cmd:0}, all other actors{p_end}
{synopt :{opt wei:ght(name)}}influence model: {cmd:invdist} (weights proportional to 1/distance; default) or {cmd:expdist} (proportional to exp(-distance)){p_end}
{synopt :{opt obv:ar(#)}}influence model: prior variance of the outcome coefficients; default {cmd:100}{p_end}
{synopt :{opt osig:ma(a b)}}influence model: inverse-gamma prior of the error variance; default {cmd:0.01 0.01}{p_end}
{synopt :{opt nsample(#)}}number of posterior draws kept (in all chains together); default {cmd:4000}{p_end}
{synopt :{opt burnin(#)}}sweeps discarded after the tuning runs; default {cmd:10000}{p_end}
{synopt :{opt interval(#)}}sweeps between two kept draws; default {cmd:10}{p_end}
{synopt :{opt pilot(#)}}number of tuning runs for the step sizes; default {cmd:4}{p_end}
{synopt :{opt plen(#)}}sweeps in each tuning run; default {cmd:500}{p_end}
{synopt :{opt mklstarts(#)}}number of starting points of the minimum Kullback-Leibler fit; default {cmd:1}{p_end}
{synopt :{opt chai:ns(#)}}number of independent chains, run in parallel and pooled; default {cmd:1}{p_end}
{synopt :{opt thr:eads(#)}}number of threads; default: one per core of the computer; the result does not depend on it{p_end}
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

{title:Latent space model for social influence}

{pstd}
With {opt influence(}{it:outcome lag} [{it:covariates}]{opt )} {cmd:nwlsm} fits the latent space model for social influence of Sweet and Adhikari (2020): the latent space model of the ties above and, jointly with it,
a linear model for an actor-level outcome in which social influence acts through the latent positions that are estimated from the network. The variables hold one observation per actor, in the order
of the actors of {it:netname}: first the outcome y_t, measured after the network was formed, second the same outcome y_(t-1) measured before (the lagged outcome), then any other regressors. The outcome of actor {it:k} is

{p 12 12 2}
y_k = b_0 + b_1 y_(t-1),k + b_2 x_k2 + ... + b_influence m_k + e_k,   e_k ~ N(0, s2_y),   m_k = sum_j w_kj y_(t-1),j

{pstd}
where m_k is the weighted mean of the lagged outcomes of the other actors, and the weights depend on the distances between the latent positions: w_kj is proportional to f(|z_k - z_j|) over
the {opt neighbors()} actors that are nearest to {it:k} in the latent space (all other actors by default) and sum to 1 for each {it:k}; f(d) = 1/d with {cmd:weight(invdist)} (the weights of the paper) or exp(-d) with {cmd:weight(expdist)}.
Actors who are close in the latent space, in practice those with many ties to {it:k} or to the same others, are therefore most influential for {it:k}; ties to {it:k} are not needed to influence {it:k}.
The coefficient b_influence is the effect of social influence (the "influence" row of the outcome equation). The lagged outcome is always a regressor.
The weights change with the positions, so the outcome also informs the positions; the network and the outcome parameters are estimated in one MCMC run.

{pstd}
Sampler: the positions are updated one actor at a time by random-walk Metropolis-Hastings, whose acceptance ratio contains the change of the log likelihood of the ties and of all outcomes
(the weights of the actor who moves and its weight in the others' means); coefficients of the outcome equation and error variance are drawn by Gibbs steps given the weights (the conditional distributions are a normal and an inverse gamma
distribution). Priors: outcome coefficients N(0, {opt obvar()}), error variance inverse gamma with shape and rate {opt osigma()}. The network part (coefficients, effects, clusters, family, priors) is the model described above. The influence
model needs the Euclidean distance; it cannot be combined with {opt bilinear}.

{pstd}
Differences from Sweet and Adhikari (2020). The intercept of the tie model is estimated (the paper fixes it at 2 to resolve a weak trade-off with the variance of the positions; the prior of the positions here is that of latentnet, whose scale grows with the number of actors).
The variance of the positions has the latentnet prior instead of an inverse gamma. The outcome coefficients and the error variance are updated by Gibbs steps (the same posterior as random-walk steps). The default prior of the error variance is inverse gamma(0.01, 0.01), weaker than the paper's.
Covariates of the ties are allowed, the paper uses none. One network and one outcome wave are modelled; the hierarchical version for several networks with random coefficients and the paper's indicator that makes actors without close neighbours immune to influence are not implemented.
The positions reported ({cmd:e(Z)}) are the minimum Kullback-Leibler positions of the tie model, as for the other models. Influence needs variation in the lagged outcome among the actors who are weighted: with all other actors as neighbours the weighted means are nearly constant
and the influence effect is weakly identified; a moderate {opt neighbors()} (as in the application of the paper, about the number of ties an actor can name) gives more information, and a {opt neighbors()} much larger than the true number of influential actors biases the effect.

{pstd}
Results: {cmd:e(b)} and {cmd:e(V)} hold the posterior means and covariance of the tie coefficients (equation {cmd:network}) and the outcome coefficients (equation {cmd:outcome}, the last is the influence effect); {cmd:e(summary)}, {cmd:e(rhat)} and
{cmd:estat mcmcdiag} include them, the error variance is the row {cmd:sigma2_y} of {cmd:e(vsummary)}. {cmd:e(W)} is the posterior mean weight matrix (row {it:k}: the weights of the actors on {it:k}) and {cmd:e(pred)} holds, for each actor, the outcome, the posterior mean predicted outcome, the
exposure m_k and the part of the prediction due to influence. {cmd:e(prob_positive)} is the posterior probability that the influence effect is positive.

{title:Valued networks}

{pstd}
With {opt family()} the values of the ties are modelled instead of their presence, as latentnet does with its {cmd:response} and {cmd:family} arguments. The linear predictor
eta_ij = b_0 + b_1 x_ij1 + ... + s_i + r_j + o_i + o_j - |z_i - z_j| (or with the inner product) is the same in all families; the families differ in how the value y_ij depends on it:

{phang2}{cmd:poisson}: y_ij ~ Poisson(exp(eta_ij)); non-negative integer values (counts).{p_end}
{phang2}{cmd:binomial}: y_ij ~ Binomial(m_ij, logistic(eta_ij)) with the number of trials m_ij given by {opt trials()}: a positive integer (the same for all dyads), a network of the same actors or a Stata matrix with one row and column per actor. The values are integers from 0 to m_ij.{p_end}
{phang2}{cmd:normal}: y_ij ~ N(eta_ij, tau2) with the dispersion tau2 (the variance of the values) having a scaled inverse chi-squared prior with scale {opt dispvar()} and {opt dispdf()} degrees of freedom; any real values. latentnet requires the prior of the dispersion to be given; here the defaults are 1 and 3.{p_end}

{pstd}
The network must carry its values: a network created with {cmd:nwset ..., valued} (or from a weighted import) is valued; a valued network without {opt family()} is refused unless all its ties equal 1. The values of the network are
used as they are: a dyad without a tie has the value 0 (a missing dyad is left out of the likelihood). Sender, receiver and sociality effects, covariates, clusters, the bilinear term, chains and threads work as in the binary model.
The priors of the coefficients, the positions and the random effects are those of the binary model. Starting values: the positions from classical scaling of a dissimilarity made from the values, the coefficients from a generalised linear model of the family given those distances. The variance of the normal family is sampled by a Gibbs step and kept as the last row of {cmd:e(vsummary)} (and as the variable {cmd:dispersion} of {opt savedraws()}).
The minimum Kullback-Leibler fit uses the family's log likelihood of the posterior mean expected values; for the normal family that objective is scaled by the posterior mean of the dispersion. With a valued family {cmd:e(P)} is the matrix of the posterior mean expected values E(y_ij) (the tie probabilities for the binary model; m_ij times the probability for the binomial family).
{cmd:estat gof} for a valued network works on the network of the dyads with a positive value, and tabulates the distribution of the values.

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
The optimisation starts from the draw whose tie probabilities are closest to the posterior mean probabilities, as in latentnet. The MKL fit can have several local optima;
{opt mklstarts(#)} starts from the {it:#} closest draws instead and keeps the best result. The fit carries a very weak penalty on the positions and
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
{synopt :{cmd:e(dim)}, {cmd:e(nsample)}, {cmd:e(burnin)}, {cmd:e(interval)}, {cmd:e(seed)}, {cmd:e(mklstarts)}}settings of the run ({cmd:e(nsample)}: draws in all chains){p_end}
{synopt :{cmd:e(chains)}, {cmd:e(nsample_chain)}, {cmd:e(threads)}, {cmd:e(clusters)}}number of chains, draws per chain, threads requested (0: one per core), number of latent clusters{p_end}
{synopt :{cmd:e(rhat)}}Gelman-Rubin statistics (several chains){p_end}
{synopt :{cmd:e(Q)}, {cmd:e(class)}, {cmd:e(cmeans)}, {cmd:e(lambda)}}latent clusters: membership probabilities, classification, cluster centres, mixing proportions{p_end}
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
{synopt :{cmd:e(influence)}, {cmd:e(outcome)}, {cmd:e(lagvar)}, {cmd:e(ocovariates)}, {cmd:e(weight)}}influence model: "yes"; the outcome and lagged outcome variables; the names of the outcome coefficients; the weight function{p_end}
{synopt :{cmd:e(neighbors)}, {cmd:e(obvar)}, {cmd:e(osigma_a)}, {cmd:e(osigma_b)}, {cmd:e(prob_positive)}, {cmd:e(ll_outcome)}}number of nearest actors, priors, posterior probability of a positive influence effect, posterior mean log likelihood of the outcomes{p_end}
{synopt :{cmd:e(W)}, {cmd:e(pred)}}posterior mean weight matrix; outcome, predicted outcome, exposure and influence part by actor{p_end}
{synopt :{cmd:e(effects)}, {cmd:e(latent)}}the random effects in the model; {cmd:euclidean} or {cmd:bilinear}{p_end}
{synopt :{cmd:e(family)}, {cmd:e(trials)}, {cmd:e(dispvar)}, {cmd:e(dispdf)}}the family; the number of trials if one number was given; the prior of the dispersion (normal family){p_end}
{p2colreset}{...}

{pstd}
The full posterior sample stays in Mata after the fit: {cmd:mata: LsmSample("beta")} (draws by coefficient), {cmd:"llk"}, {cmd:"s2"} (one column for each cluster with {opt clusters()}), {cmd:"class"} (the corrected memberships by draw), {cmd:"lambda"}, {cmd:"chains"}, {cmd:"ev"}
(variances of the sender, receiver and sociality effects), {cmd:"E"} (the effects by draw), {cmd:"Z"}
(one row per draw: the positions by actor, then dimension), {cmd:"Zalign"} (the aligned mean positions) and {cmd:"mkl"} (the MKL positions).

{title:Chains, threads and speed}

{pstd}
With {opt chains(#)} the command runs {it:#} independent chains and pools their draws. {opt nsample()} is the number of draws in all chains together (divided between the chains, rounded up); every chain
has its own burn-in and tuning runs, so that more chains shorten the run only as far as the burn-in allows. The first chain starts from the classical-scaling positions; the others start
from over-dispersed versions of them, which is what the Gelman-Rubin statistic needs. Chain {it:k} uses the seed {opt seed()} + 1000003 {it:k}. Positions of different chains are rotated differently, so
only rotation-invariant quantities are compared between chains: {cmd:estat mcmcdiag} shows the Gelman-Rubin potential scale reduction factor of the coefficients, the log likelihood and the variance
parameters, and {cmd:e(rhat)} holds it; values close to 1 (below 1.1, say) indicate that the chains agree. The minimum Kullback-Leibler fit uses the pooled draws.

{pstd}
{opt threads()} sets how many threads the plugin uses (default: one per core). Chains are distributed over the threads first. Threads that are left over (more threads than chains) work
inside a chain on the passes over all dyads (the updates of the coefficients) and on the minimum Kullback-Leibler fit, in networks of about 100 actors or more. A single chain cannot be
made much faster than that, because its updates of the positions are sequential. The results for a given seed and {opt chains()} are exactly the same for every value of {opt threads()}.

{title:Latent clusters}

{pstd}
{opt clusters(G)} replaces the single normal distribution of the positions by a mixture of {it:G} spherical normal distributions as in latentnet: z_i ~ N(mu_g, s2_g I) for actor {it:i} in cluster {it:g},
mu_g ~ N(0, omega2 I), s2_g with a scaled inverse chi-squared prior, and mixing proportions with a Dirichlet prior. The defaults follow latentnet: scale 0.125 (n/G)^(2/d) and sqrt(n/G) degrees of freedom for the cluster variances,
omega2 = 2 x 0.125 n^(2/d), and Dirichlet parameter sqrt(n/G). The memberships are sampled together with the other parameters. The labels of the clusters are arbitrary in every draw; they are corrected by
permuting the labels of each draw to agree as well as possible with the posterior membership probabilities (first those of the draw that started the minimum Kullback-Leibler fit, then three rounds with the probabilities found so far).
The result: {cmd:e(Q)} (posterior probability that an actor belongs to a cluster), {cmd:e(class)} (the most probable cluster of each actor), {cmd:e(cmeans)} (posterior mean cluster centres, aligned to the
minimum Kullback-Leibler positions), {cmd:e(lambda)} (mixing proportions: mean and sd) and one variance for each cluster in {cmd:e(vsummary)}. The classification and the centres are label-corrected and aligned
posterior means of the draws; they are not the mixture fitted to the positions by latentnet's minimum Kullback-Leibler procedure.

{title:Post-estimation}

{phang}
{cmd:estat positions} [{cmd:, mean generate(}{it:stub}{cmd:) nolist}] lists the MKL positions with the posterior mean effects ({opt mean}: the mean of the aligned draws instead) and returns them in {cmd:r(Z)}.
With {opt generate()} the coordinates go into the variables {it:stub}1, {it:stub}2, ... of a data set with one observation per actor, in the order of the actors.

{phang}
{cmd:estat mcmcdiag} [{cmd:, nolist graph name(}{it:name}{cmd:)}] shows, for the coefficients, the log likelihood and the variance parameters, the posterior mean and standard deviation, the effective sample size
(initial monotone positive sequence estimator of Geyer), the Geweke z-score (first 10 percent of the draws against the last 50 percent) and, with several chains, the Gelman-Rubin statistic, and returns them in {cmd:r(diag)}. It needs the draws of the last fit, which stay in memory. With {opt graph} the trace plots of all these quantities are drawn (one line per chain, with the Gelman-Rubin statistic in the title) and combined into one graph.

{phang}
{cmd:estat clusters} [{cmd:, generate(}{it:stub}{cmd:) nolist}] (after {opt clusters()}) lists the posterior membership probabilities and the classification and returns them in {cmd:r(Q)} and {cmd:r(class)}. With {opt generate()} the variables
{it:stub}_class and {it:stub}_p1, {it:stub}_p2, ... are added to a data set with one observation per actor.

{phang}
{cmd:estat predict} [{cmd:, generate(}{it:netname}{cmd:) replace}] returns the posterior mean tie probabilities in {cmd:r(P)} and, with {opt generate()}, stores them as a valued network of the family.

{phang}
{cmd:estat gof} [{cmd:, nsim(}{it:#}{cmd:) seed(}{it:#}{cmd:) gof(}{it:terms}{cmd:) plot maxdeg(}{it:#}{cmd:) maxdist(}{it:#}{cmd:) name(}{it:name}{cmd:) threads(}{it:#}{cmd:)}] is the goodness of fit of latentnet's {cmd:gof()}: one network is simulated from the model at each of {opt nsim()}
(default 100) posterior draws, taken evenly spaced from the stored draws, and statistics of the observed network are compared with their distribution over the simulated networks. {opt gof()} lists the statistics:
{cmd:degree} (undirected), {cmd:idegree} and {cmd:odegree} (directed), {cmd:espartners} (edgewise shared partners; for directed networks the outgoing two-paths i-k-j), {cmd:dspartners} (dyadwise shared partners), {cmd:distance} (geodesic distances, the last cell is unreachable pairs),
{cmd:model} (the total value of the ties and the sum of each covariate over the ties) and {cmd:values} (valued networks: the number of dyads with each value; for the normal family in ten bins of equal observed frequency). The default is latentnet's:
{cmd:degree distance} (undirected) or {cmd:idegree odegree distance} (directed), plus {cmd:values} for the Poisson and binomial families, and {cmd:model values} for the normal family (whose network structure is not defined).
For each cell the table shows the observed count, the minimum, mean and maximum over the simulated networks and the Monte Carlo p-value 2 min(P(simulated <= observed), P(simulated >= observed)), at most 1, as in latentnet.
With {opt plot} the shares of the observed and the simulated networks are drawn in box plots (the observed network as a red line). The tables are returned in {cmd:r(gof_degree)}, {cmd:r(gof_distance)} and so on. Missing dyads of the observed network count as no tie.
Deviation from latentnet: its {cmd:simulate()} draws an undirected network by drawing every ordered pair and keeping a tie if either draw is a tie, so that its simulated undirected networks are too dense and its undirected {cmd:gof()} is biased;
{cmd:estat gof} draws one value for every pair. For directed networks the simulation agrees with latentnet's.

{phang}
{cmd:estat ppcheck} [{cmd:, groups(}{it:#}{cmd:) graph name(}{it:name}{cmd:) nolist}] is a posterior predictive check of the tie probabilities: the dyads are sorted by their posterior mean predicted probability (expected value) and cut into {opt groups()} (default 10) groups of about equal size; the table compares the mean prediction with the observed mean
in every group. Below it: for binary networks the Brier score, the area under the ROC curve and the mean log loss; for valued networks the root mean squared error, the correlation and the mean absolute error of the posterior mean expected values with the observed values. {opt graph} draws the calibration plot. Returned in {cmd:r(groups)}, {cmd:r(score)}, {cmd:r(auc_or_corr)} and {cmd:r(N)}.

{phang}
{cmd:estat influence} [{cmd:, top(}{it:#}{cmd:) predict(}{it:var}{cmd:) exposure(}{it:var}{cmd:) exerted(}{it:var}{cmd:) graph name(}{it:name}{cmd:) nolist}] (after {opt influence()}) shows the outcome equation, the posterior probability of a positive influence effect,
the R-squared of the posterior mean predictions with and without the influence term, the mean effective number of influencing actors (the inverse of the sum of squared weights) and the {opt top()} actors who exert most influence (column sums of the posterior mean weights) with the exposure they receive.
{opt predict()}, {opt exposure()} and {opt exerted()} store the posterior mean prediction, the exposure m_k and the influence exerted in new variables (one observation per actor); {opt graph} plots observed against predicted outcomes. Results in {cmd:r(influence)} (posterior mean), {cmd:r(influence_sd)},
{cmd:r(influence_lower)}, {cmd:r(influence_upper)}, {cmd:r(prob_positive)}, {cmd:r(r2)}, {cmd:r(r2_noinfluence)}, {cmd:r(eff_neighbors)} and {cmd:r(exerted)}.

{pstd}
{help nwlsmplot} plots the positions with the ties as light lines, the actors coloured by a variable and sized by an effect, the influence exerted ({cmd:size(influence)}) or a variable.

{title:Examples}

{phang2}{cmd:. nwlsm friends, dim(2)}{p_end}
{phang2}{cmd:. nwlsm friends, dim(2) absdiff(age) nsample(2000) burnin(5000) seed(12)}{p_end}
{phang2}{cmd:. nwlsm advice, dim(2) sender receiver}{p_end}
{phang2}{cmd:. nwlsm friends, dim(2) sociality}{p_end}
{phang2}{cmd:. estat mcmcdiag}{p_end}
{phang2}{cmd:. estat positions, generate(pos)}{p_end}
{phang2}{cmd:. estat predict, generate(pfriends)}{p_end}
{phang2}{cmd:. nwlsm friends, dim(2) chains(4) nsample(4000)}{p_end}
{phang2}{cmd:. nwlsm friends, dim(2) clusters(3)}{p_end}
{phang2}{cmd:. estat clusters, generate(grp)}{p_end}
{phang2}{cmd:. nwlsm counts, family(poisson) dim(2)}{p_end}
{phang2}{cmd:. nwlsm votes, family(binomial) trials(10) sender receiver}{p_end}
{phang2}{cmd:. nwlsm ratings, family(normal) dispvar(1) dispdf(3)}{p_end}
{phang2}{cmd:. estat gof, nsim(200) plot}{p_end}
{phang2}{cmd:. estat ppcheck, graph}{p_end}
{phang2}{cmd:. estat mcmcdiag, graph}{p_end}
{phang2}{cmd:. nwlsm advice, dim(2) influence(belief2 belief1 tenure) neighbors(10)}{p_end}
{phang2}{cmd:. estat influence, top(5) predict(pbelief)}{p_end}
{phang2}{cmd:. nwlsmplot, size(influence)}{p_end}
{phang2}{cmd:. nwlsmplot, clusters}{p_end}
{phang2}{cmd:. nwlsmplot, color(gender) size(sociality)}{p_end}

{title:Reference}

{p 4 8 2}Hoff, P. D., Raftery, A. E. and Handcock, M. S. 2002. Latent space approaches to social network analysis. Journal of the American Statistical Association 97: 1090-1098.{p_end}
{p 4 8 2}Handcock, M. S., Raftery, A. E. and Tantrum, J. M. 2007. Model-based clustering for social networks. Journal of the Royal Statistical Society A 170: 301-354.{p_end}
{p 4 8 2}Sweet, T. M. and Adhikari, S. 2020. A latent space network model for social influence. Psychometrika 85: 251-274.{p_end}
{p 4 8 2}Krivitsky, P. N. and Handcock, M. S. 2008. Fitting position latent cluster models for social networks with latentnet. Journal of Statistical Software 24(5).{p_end}
{p 4 8 2}Krivitsky, P. N., Handcock, M. S., Raftery, A. E. and Hoff, P. D. 2009. Representing degree distributions, clustering, and homophily in social networks with latent cluster random effects models. Social Networks 31: 204-213.{p_end}

{hline}
{pstd}{it:nwlsm version 0.6.0, 19 November 2026}{p_end}{* nwversion}
