<html>
 <head>
  <title>DateLife</title>
 </head>
 <body>
<p><a href="http://www.evoio.org/wiki/Phylotastic"><img src="datelife1.png"></a></p>
<p>Use the syntax <code>Genus_species,Genus_species</code>: comma-delimited, no spaces, underscores in names</p>
 <p><form action = "http://datelife.org/cgi-bin/R/result" method="get">
 <p>Taxa: <input type="text" name="taxa" value="Rhinoceros_unicornis,Equus_caballus" size="80"></p>

 <p>Return format: <select name="format">
 <option value="html">HTML</option>
 <option value="newick">Newick</option>
 <option value="bestguess">The single best guess</option>
 <option value="nexml">NeXML (in the future)</option>
 <option value="rda">R data file (in the future)</option>
 </select> (note that to return a tree you'll need to provide at least three taxon names)</p>
 
 <p>Partial match: <select name="partial">
 <option value="liberal">Return info from trees that have some taxon overlap with your query</option>
 <option value="conservative">Return info only from trees with complete taxon overlap with your query</option>
 </select></p>
 
 <p>Use embargoed data: <select name="useembargoed">
 <option value="yes">Yes, use all data in the database, even that otherwise hidden until publication</option>
 <option value="no">No, use public peer-reviewed results, only</option>
 </select></p>
 
 
 <p><input type="submit" value="Send"></p>
</form></p>
<br /><hr>
<p>Note source code is available at <a href="https://bitbucket.org/bomeara/datelife">https://bitbucket.org/bomeara/datelife</a></p>
<p>We use the <a href="http://r-forge.r-project.org/projects/phyloorchard/">PhyloOrchard</a> R package to store trees; you can install it by doing<br />
<code>install.packages("PhyloOrchard", repos="http://R-Forge.R-project.org")</code></p>
<br /><hr>
<p>Work on code:</p>
<p><script id="feed-1338987350323483" type="text/javascript" src="http://rss.bloople.net/?url=https%3A%2F%2Fbitbucket.org%2Fbomeara%2Fdatelife%2Frss&showtitle=false&type=js&id=1338987350323483"></script></p> </body>
</html>
