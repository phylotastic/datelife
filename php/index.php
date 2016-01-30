<html>
<?php include("pagestart.html"); ?>
<p>This will give the age of the most recent common ancestor of a set of taxa, or stretch an input tree, based on available time-calibrated trees. It is still actively in development, so it may break at times. Currently it mostly has mammals in the database, but we are actively adding taxa. See the <a href="http://datelife.org/faq.php">FAQ</a> page for more info on this project, and the <a href="http://phylotastic.org">Phylotastic</a> page for more info. This was started at a NESCent hackathon, added to at others, and eventually NSF-funded as part of Phylotastic.</p>
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

<p>Taxonomic name resolution (converting input names to a standard taxonomy; will slow down analysis (~0.5s/taxon) but increase chance of matching): <select name="usetnrs">
<option value="yes">yes</option>
<option value="no" selected="selected">No</option>
</select>

<p>Use approximate name matching (makes matching MUCH slower: perhaps 20s/taxon): <select name="approximatematch">
<option value="yes" selected="selected">yes</option>
<option value="no">No</option>
</select>

<p>Prune taxa that don't match (otherwise, keep their original names): <select name="prunenonmatch">
<option value="yes" selected="selected">yes</option>
<option value="no">No</option>
</select>


 <p><input type="submit" value="Send"></p>
</form></p>
<br /><hr>
<p>Note source code is available at <a href="https://github.com/phylotastic/datelife">https://github.com/phylotastic/datelife</a></p>
<br />
<hr />
<p>DateLife is part of the Phylotastic project, funded by NSF. It uses resources from the Open Tree of Life project, rOpenSci, and TreeBase. To help, please <a href="https://tree.opentreeoflife.org/curator">add trees</a> to OpenTree's database! [they can be yours or others']</p>
<p><center><a href="http://www.nsf.gov"><img src="https://www.nsf.gov/images/logos/nsf1.jpg" width=50 height=50></a></center></p>

<?php include("pageend.html"); ?>
