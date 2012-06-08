<html>
<?php include("pagestart.html"); ?>

<p>Use the syntax <code>Genus_species,Genus_species</code>: comma-delimited, no spaces, underscores in names</p>
 <p><form action = "http://datelife.org/cgi-bin/R/result" method="get">
 <p>Taxa: <input type="text" name="taxa" value="Rhinoceros_unicornis,Equus_caballus" size="80"></p>

 <p>Return format: <select name="format">
 <option value="html">HTML</option>
 <option value="newick">Newick (in the future)</option>
 <option value="bestguess">Best guess</option>
 <option value="bestguessuncert">Best guess,Min guess,Max guess</option>
 <option value="nexml">NeXML (in the future)</option>
 <option value="rda">R data file (in the future)</option>
 </select><br />Note that to return a tree you'll need to provide at least three taxon names. The "best guesses" are just the medians of the results by study.</p>

 <p>Partial match: <select name="partial">
 <option value="liberal">Return info from trees that have some taxon overlap with your query</option>
 <option value="conservative">Return info only from trees with complete taxon overlap with your query</option>
 </select></p>

 <p>Use embargoed data: <select name="useembargoed">
 <option value="yes">Yes, use all data in the database, even data that is otherwise being hidden until publication</option>
 <option value="no">No, use public peer-reviewed results, only</option>
 </select></p>

 <p>How to deal with single point estimates (which have no uncertainty): <select name="uncertainty">
 <option value="100">Assume +/- 100% of the age</option> 
 <option value="75">Assume +/- 75% of the age</option> 
 <option value="50">Assume +/- 50% of the age</option> 
 <option value="20">Assume +/- 20% of the age</option> 
 <option value="10">Assume +/- 10% of the age</option> 
 <option value="5">Assume +/- 5% of the age</option> 
 <option value="1">Assume +/- 1% of the age</option> 
</select></p>


 <p><input type="submit" value="Send"></p>
</form></p>
<br /><hr>
<p>Note source code is available at <a href="https://bitbucket.org/bomeara/datelife">https://bitbucket.org/bomeara/datelife</a></p>
<p>We use the <a href="http://r-forge.r-project.org/projects/phyloorchard/">PhyloOrchard</a> R package to store trees; you can install it by doing<br />
<code>install.packages("PhyloOrchard", repos="http://R-Forge.R-project.org")</code></p>

<?php include("pageend.html"); ?>
