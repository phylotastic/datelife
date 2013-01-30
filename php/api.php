<html>
<?php include("pagestart.html"); ?>
<h2>API</h2>
<br />
<p>This page describes the current API for DateLife. Note that it is still alpha: it may break or change.
<br />The basic interface is passing arguments through the url (something familiar from other sites: looking at <a href="http://www.ncbi.nlm.nih.gov/nuccore/?term=myrmecocystus">http://www.ncbi.nlm.nih.gov/nuccore/?term=myrmecocystus</a>, you can tell that <a href="http://www.ncbi.nlm.nih.gov/nuccore/?term=rhinostomus">http://www.ncbi.nlm.nih.gov/nuccore/?term=rhinostomus</a> would give info on a different taxon). The arguments correspond to info on the <a href="http://www.datelife.org">home page</a>:</p>

<p>format=<br />
<table border="0">
<tr>
<td>html</td>
<td>Web page with tables and a figure</td>
</tr>
<tr>
<td>newickmed</td>
<td>Newick median tree (median of all the chronograms with sufficient overlap)</td>
</tr>
<tr>
<td>bestguess</td>
<td>The median estimate for the root age: a single floating point number</td>
</tr>
<tr>
<td>bestguessuncert</td>
<td>Comma delimited: Median,Min,Max estimates for the root age</td>
</tr>
</table>
</p>
<?php include("pageend.html"); ?>
