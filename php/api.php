<html>
<?php include("pagestart.html"); ?>
<h2>API</h2>
<br />
<p>This page describes the current API for DateLife. Note that it is still alpha: it may break or change.</p>
<p>The basic interface is passing arguments through the URL. For example, <br /><a href="http://datelife.org/cgi-bin/R/result?input=Rhinoceros_unicornis%2CEquus_caballus%2CMus_musculus&format=newickmed&partial=liberal&useembargoed=yes&uncertainty=100">http://datelife.org/cgi-bin/R/result?input=Rhinoceros_unicornis%2CEquus_caballus%2CMus_musculus&format=newickmed&partial=liberal&useembargoed=yes&uncertainty=100</a><br /> will return a newick tree linking <i>Rhinoceros unicornis</i>, <i>Equus caballus</i>, and <i>Mus musculus</i>. The arguments correspond to info on the <a href="http://www.datelife.org">home page</a>:</p>

<p><b>format=</b>
<br /><b>html</b>: The standard web page view of results
<br /><b>newickmed</b>: The Newick tree that is the median of all the returned tree ages
<br /><b>bestguess</b>: A single floating point number with the median age of the MRCA across studies
<br /><b>bestguessuncert</b>: Comma-delimited: median age, min age, max age
</p>
<?php include("pageend.html"); ?>
