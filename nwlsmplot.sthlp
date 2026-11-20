{smcl}
{* *! version 1.0.3  20nov2026 author: Thomas Grund}{...}
{title:Title}

{p2colset 9 20 27 2}{...}
{p2col :nwlsmplot  {hline 2}}Plot the latent positions of an {help nwlsm} fit{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:nwlsmplot} [{cmd:,} {it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt dims(# #)}}the two latent dimensions on the axes; default {cmd:1 2}{p_end}
{synopt :{opt mean}}plot the mean of the aligned draws {cmd:e(Zmean)} instead of the minimum Kullback-Leibler positions {cmd:e(Z)}{p_end}
{synopt :{opt color(varname)}}colour the actors by the values of a numeric or string variable (at most 15 values; the groups of a string variable are in alphabetical order); the legend shows its value labels, if it has them, and its variable label as the legend title{p_end}
{synopt :{opt colors(list)}}colours of the groups of {opt color()} or {opt clusters}: a list of Stata colour names, or {cmd:gray}, where the groups differ by gray value and by marker (filled and hollow circles, diamonds, triangles, squares) so that they stay apart in a black-and-white print{p_end}
{synopt :{opt clus:ters}}colour the actors by their latent cluster (after {cmd:nwlsm, clusters()}); not with {opt color()}{p_end}
{synopt :{opt size(spec)}}size the actors by {cmd:sender}, {cmd:receiver} or {cmd:sociality} (posterior mean effects of the fit), by {cmd:influence} (the influence each actor exerts, after the influence model) or by a numeric variable{p_end}
{synopt :{opt label(varname)}}label the actors with a string variable, or with a numeric variable (its value labels, or the numbers where it has none){p_end}
{synopt :{opt labs:ize(string)}}size of the labels; default {cmd:vsmall}{p_end}
{synopt :{opt labp:osition(string)}}position of the labels as a clock position (0 to 12), or {cmd:auto}; default {cmd:auto}{p_end}
{synopt :{opt noedges}}do not draw the ties{p_end}
{synopt :{opt arr:ows}}draw the ties of a directed network as arrows from sender to receiver (the default for directed networks){p_end}
{synopt :{opt noarr:ows}}draw the ties of a directed network as lines without arrow heads{p_end}
{synopt :{opt arrowsize(string)}}size of the arrow heads ({cmd:vsmall}, {cmd:small}, ...); default {cmd:small}, {cmd:vsmall} when the network has more than 400 ties{p_end}
{synopt :{opt edgecolor(string)}}colour of the tie lines; default {cmd:gs11%60}, a lighter {cmd:gs12%35} when the network has more than 400 ties{p_end}
{synopt :{opt edgewidth(string)}}width of the tie lines; default {cmd:vthin}{p_end}
{synopt :{opt msize(string)}}marker size when {opt size()} is not given; default {cmd:medium}{p_end}
{synopt :{opt clusterm:eans}}mark the posterior mean cluster means {cmd:e(cmeans)} with a plus sign (fit with {cmd:clusters()}){p_end}
{synopt :{opt circ:les}}draw a circle around each cluster mean with radius equal to the square root of the posterior mean variance of the cluster (fit with {cmd:clusters()}){p_end}
{synopt :{opt cert:ainty}}with {opt clusters}: opacity of the markers by the certainty of the cluster membership, the largest posterior membership probability {cmd:e(Q)} (four steps){p_end}
{synopt :{opt draws(#)}}draw a faint cloud of {it:#} (1 to 1000) evenly spaced posterior draws of the positions under the actors{p_end}
{synopt :{opt noax:es}}draw no axes{p_end}
{synopt :{opt noeq:ual}}do not force equal scales: the plot is square and each axis spans its own range (by default one unit has the same length on both axes){p_end}
{synopt :{opt exp:ort(filename [, options])}}export the graph with {help graph export}, replacing an existing file; the extension of {it:filename} sets the format{p_end}
{synopt :{opt name(name)}}name of the graph{p_end}
{synopt :{it:twoway_options}}other options of {help twoway}, such as {cmd:title()}, {cmd:scheme()} or {cmd:saving()}; they come last and override the defaults. Titles, subtitles, axis titles and the legend title are black by default, also a {cmd:title()} or {cmd:subtitle()} you give (unless it has its own {cmd:color()}){p_end}
{synoptline}

{title:Description}

{pstd}
{cmd:nwlsmplot} draws the positions of the actors after {cmd:nwlsm} in the plane of two latent dimensions. The ties of the network are drawn as light lines (an undirected
tie once, a directed tie as an arrow from sender to receiver, or with {opt noarrows} as a line). {opt color()}, {opt size()} and {opt label()} take variables with one observation per
actor, in the order of the actors of the network, as the other plot commands of the family do. The positions are identified only up to rotation and reflection, so that the
orientation of the axes carries no meaning. Because only distances carry meaning, the two axes have the same scale (one unit is equally long on both), so that distances and the variance circles are drawn undistorted; the plot region takes the shape of the range of the points, circles, draws and cluster means. {opt noequal} switches this off; an {opt aspectratio()} given as a twoway option replaces the equal-scale aspect ratio.

{pstd}
{opt draws()} shows the uncertainty of the positions: the draws of the sampler kept by the last {cmd:nwlsm} fit are centred and rotated or reflected to the minimum Kullback-Leibler
positions, as in the computation of {cmd:e(Zmean)}, and {it:#} evenly spaced draws are drawn as small dots, coloured as the actors. The draws stay in memory until the next fit or {cmd:clear all}; if
they are gone, or belong to another fit than the one in {cmd:e()}, {cmd:nwlsmplot} stops with an error. {opt clustermeans} and {opt circles} are the cluster means and the spherical-normal variance circles
of the latent cluster model, as {cmd:plot.means} and {cmd:plot.vars} of latentnet's {cmd:plot.ergmm()}; the circle is the cluster in the plane of the two chosen dimensions.
The cluster labels of the fit are those of {cmd:e(class)}, {cmd:e(cmeans)} and {cmd:e(lambda)}. {opt certainty} gives markers with a clearly assigned membership full colour and uncertain ones a pale colour.
{p_end}

{pstd}
Not offered, compared with {cmd:plot.ergmm()}: pie charts of the membership probabilities (see {opt certainty}), density plots, plots of one latent dimension (jittered), 3-D plots, {cmd:Z.ref} and {cmd:zoom.on}.
{p_end}

{title:Examples}

{phang2}{cmd:. nwlsm friends, dim(2) sociality}{p_end}
{phang2}{cmd:. nwlsmplot}{p_end}
{phang2}{cmd:. nwlsmplot, clusters}{p_end}
{phang2}{cmd:. nwlsmplot, color(gender) size(sociality) label(name) labsize(small) labposition(6) title("Friendship space")}{p_end}
{phang2}{cmd:. nwlsm friends, dim(2) clusters(3)}{p_end}
{phang2}{cmd:. nwlsmplot, clusters clustermeans circles certainty draws(100) noedges}{p_end}
{phang2}{cmd:. nwlsmplot, arrows export(positions.png, width(1200))}{p_end}

{title:Also see}

{psee}{help nwlsm}, {help nwplot}{p_end}

{hline}
{pstd}{it:nwlsm version 1.0.3, 20 November 2026}{p_end}{* nwversion}
