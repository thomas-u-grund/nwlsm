{smcl}
{* *! version 0.6.0  19nov2026 author: Thomas Grund}{...}
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
{synopt :{opt color(varname)}}colour the actors by the values of a variable (at most 15 values){p_end}
{synopt :{opt clus:ters}}colour the actors by their latent cluster (after {cmd:nwlsm, clusters()}); not with {opt color()}{p_end}
{synopt :{opt size(spec)}}size the actors by {cmd:sender}, {cmd:receiver} or {cmd:sociality} (posterior mean effects of the fit), by {cmd:influence} (the influence each actor exerts, after the influence model) or by a numeric variable{p_end}
{synopt :{opt label(varname)}}label the actors with a string variable{p_end}
{synopt :{opt noedges}}do not draw the ties{p_end}
{synopt :{opt edgecolor(string)}}colour of the tie lines; default {cmd:gs11%60}{p_end}
{synopt :{opt edgewidth(string)}}width of the tie lines; default {cmd:vthin}{p_end}
{synopt :{opt msize(string)}}marker size when {opt size()} is not given; default {cmd:medium}{p_end}
{synopt :{opt name(name)}}name of the graph{p_end}
{synopt :{it:twoway_options}}other options of {help twoway}, such as {cmd:title()} or {cmd:scheme()}{p_end}
{synoptline}

{title:Description}

{pstd}
{cmd:nwlsmplot} draws the positions of the actors after {cmd:nwlsm} in the plane of two latent dimensions. The ties of the network are drawn as light lines (an undirected
tie once, a directed tie as a line from sender to receiver without arrow head). {opt color()}, {opt size()} and {opt label()} take variables with one observation per
actor, in the order of the actors of the network, as the other plot commands of the family do. The positions are identified only up to rotation and reflection, so that the
orientation of the axes carries no meaning.

{title:Examples}

{phang2}{cmd:. nwlsm friends, dim(2) sociality}{p_end}
{phang2}{cmd:. nwlsmplot}{p_end}
{phang2}{cmd:. nwlsmplot, clusters}{p_end}
{phang2}{cmd:. nwlsmplot, color(gender) size(sociality) label(name) title("Friendship space")}{p_end}

{title:Also see}

{psee}{help nwlsm}, {help nwplot}{p_end}

{hline}
{pstd}{it:nwlsm version 0.6.0, 19 November 2026}{p_end}{* nwversion}
