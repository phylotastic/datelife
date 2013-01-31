<html>
<?php include("pagestart.html"); ?>
<h2>API</h2>
<br />
<p>This page describes the current API for DateLife. Note that it is still alpha: it may break or change.</p>
<p>The basic interface is passing arguments through the URL. For example, <br /><a href="http://datelife.org/cgi-bin/R/result?input=Rhinoceros_unicornis%2CEquus_caballus%2CMus_musculus&format=newickmed&partial=liberal&useembargoed=yes&uncertainty=100">http://datelife.org/cgi-bin/R/result?input=Rhinoceros_unicornis%2CEquus_caballus%2CMus_musculus&format=newickmed&partial=liberal&useembargoed=yes&uncertainty=100</a><br /> will return a newick tree linking <i>Rhinoceros unicornis</i>, <i>Equus caballus</i>, and <i>Mus musculus</i>. The arguments correspond to info on the <a href="http://www.datelife.org">home page</a>:</p>

<p><b>format=</b>
<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>html</b>: The standard web page view of results
<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>newickmed</b>: The Newick tree that is the median of all the returned tree ages
<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>bestguess</b>: A single floating point number with the median age of the MRCA across studies
<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>bestguessuncert</b>: Comma-delimited: median age, min age, max age
<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>png</b>: A PNG image of the median tree.
</p>
<p><b>partial=</b>
<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>liberal</b>: Return info from trees that have any taxon overlap with your query (even though this may result in an age underestimate)
<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>conservative</b>: Return info from trees that have complete taxon overlap with your query (even though this may be few, or zero, trees)
</p>
<p><b>useembargoed=</b>
<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>yes</b>: Yes, use all data in the database (even info from as yet unpublished trees by DateLife contributors)
<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>no</b>: No, use only peer-reviewed, published results (almost all of the dataset, except for a few big trees)
</p>
<p><b>uncertainty=</b>
<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b><i>0-100</i></b>: You can enter any number from 0 to 100 here. If there is a single point estimate, how uncertain should the age estimate be interpreted as (60=age of MRCA +/- 60% of this age)
</p>
<p><b>plot.width=</b>
<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b><i>integer</i></b>: The width of the plot, in pixels.
</p>
<p><b>plot.height=</b>
<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b><i>integer</i></b>: The height of the plot, in pixels.
</p>
<?php include("pageend.html"); ?>
