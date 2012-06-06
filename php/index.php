<html>
 <head>
  <title>DateLife</title>
 </head>
 <body>
<p><a href="http://www.evoio.org/wiki/Phylotastic"><img src="http://www.evoio.org/wg/evoio/images/f/f1/Phylotastic_logo.png"></a></p>
<p>Demo after a day of coding: enter a string of mammal names. It will give you the age of their most recent common ancestor from the Bininda-Emonds et al. mammal tree.</p>
<p>Use the syntax <code>Genus_species,Genus_species</code>: comma-delimited, no spaces, underscores in names</p>
 <p><form action = "http://datelife.org/cgi-bin/R/test" method="get">
 <p>Taxa: <input type="text" name="taxa" value="Homo_sapiens,Pan_paniscus,Mus_musculus" size="80"> <input type="submit" value="Send"></p>
</form></p>
<p>Note source code is available at <a href="https://bitbucket.org/bomeara/datelife">https://bitbucket.org/bomeara/datelife</a></p>
<p>We use the <a href="http://r-forge.r-project.org/projects/phyloorchard/">PhyloOrchard</a> R package to store trees; you can install it by doing<br />
<code>install.packages("PhyloOrchard", repos="http://R-Forge.R-project.org")</code></p>
<br /><hr>
<p>Work on code:</p>
<p><script id="feed-1338987350323483" type="text/javascript" src="http://rss.bloople.net/?url=https%3A%2F%2Fbitbucket.org%2Fbomeara%2Fdatelife%2Frss&showtitle=false&type=js&id=1338987350323483"></script></p> </body>
</html>
