<html>
<?php include("pagestart.html"); ?>
<p>This will give the age of the most recent common ancestor of a set of taxa, or stretch an input tree, based on available time-calibrated trees. It is still actively in development, so it may break at times. Currently it mostly has mammals in the database, but we are actively adding taxa. See the <a href="http://datelife.org/faq.php">FAQ</a> page for more info on this project, and the <a href="http://phylotastic.org">Phylotastic</a> page for more info on the overall hackathon, of which this was just one outcome.</p>
 <p><form action = "http://datelife.org/cgi-bin/R/result" method="get">
 <p>Use the syntax <code>Genus_species,Genus_species</code>: comma-delimited, no spaces, underscores in names, i.e., <br /><b>Rhinoceros_unicornis,Equus_caballus,Mus_musculus</b><br />for a list of species in the taxa box, <b>OR</b> a Newick tree string, ending with a semicolon: <br /><b>((Rhinoceros_unicornis,Equus_caballus),Mus_musculus);</b><br />The tree can have branch lengths. Note that stretching a tree is much slower than just getting ages for a list of taxa. It uses a method by <a href="http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12051/abstract">Eastman et al. 2013</a></p>
 <p>Input: <input type="text" name="input" value="" size="80"></p>
 
 <p>Return format: <select name="format">
 <option value="html">HTML</option>
 <option value="newickmed">Newick median</option>
 <option value="bestguess">Best guess</option>
 <option value="bestguessuncert">Best guess,Min guess,Max guess</option>
 <option value="png">Image file of median tree (PNG format)</option>
 </select><br />Note that to return a tree you'll need to provide at least three taxon names. The "best guesses" are just the medians of the results by study.</p>

 <p>Partial match: <select name="partial">
 <option value="liberal">Return info from trees that have some taxon overlap with your query</option>
 <option value="conservative">Return info only from trees with complete taxon overlap with your query</option>
 </select></p>

 <p>Use embargoed data: <select name="useembargoed">
 <option value="yes">Yes, use all data in the database (even info from as yet unpublished trees by DateLife contributors)</option>
 <option value="no">No, use public, peer-reviewed results, only</option>
 </select></p>

 <p>How to deal with single point estimates (which have no uncertainty): <select name="uncertainty">
 <option value="100">Assume +/- 100% of the age</option> 
 <option value="75">Assume +/- 75% of the age</option> 
 <option value="50">Assume +/- 50% of the age</option> 
 <option value="20">Assume +/- 20% of the age</option> 
 <option value="10">Assume +/- 10% of the age</option> 
 <option value="5">Assume +/- 5% of the age</option> 
 <option value="1">Assume +/- 1% of the age</option> 
 <option value="0">Assume +/- 0% of the age</option> 
</select></p>

<p>Taxonomic name resolution (converting input names to a standard taxonomy; will slow down analysis but increase chance of matching): <select name="usetnrs">
<option value="no">No</option>
</select><br />TNRS source:<select name="tnrssource">
<option value="NCBI">NCBI (good for wide range of species)</option>
<option value="iPlant_TNRS">iPlant (good for plants, can fix typos)</option>
<option value="MSW3">Mammal Species of the World, Vol. 3</option>
</select><br />Name resolution uses the <a href="http://taxosaurus.org">Taxosaurus</a> service.</p>


 <p><input type="submit" value="Send"></p>
</form></p>
<br /><hr>
<p>Note source code is available at <a href="https://bitbucket.org/bomeara/datelife">https://bitbucket.org/bomeara/datelife</a></p>
<p>We use the <a href="http://r-forge.r-project.org/projects/phyloorchard/">PhyloOrchard</a> R package to store trees; you can install it by doing<br />
<code>install.packages("PhyloOrchard", repos="http://R-Forge.R-project.org")</code></p>

<?php include("pageend.html"); ?>
