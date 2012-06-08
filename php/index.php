<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<!--

	zenlike1.0 by nodethirtythree design
	http://www.nodethirtythree.com

-->
<html>
<head>
<meta http-equiv="content-type" content="text/html; charset=iso-8859-1" />
<title>DateLife</title>
<meta name="keywords" content="" />
<meta name="description" content="" />
<link rel="stylesheet" type="text/css" href="default.css" />
</head>
<body>

<div id="upbg"></div>

<div id="outer">


	


	<div id="headerpic"></div>

	
	<div id="menu">
		<!-- HINT: Set the class of any menu link below to "active" to make it appear active -->
		<ul>
			<li><a href="#">Home</a></li>
			<li><a href="#">Source code</a></li>
			<li><a href="#">About</a></li>
			<li><a href="#">Contact</a></li>
			<li><a href="#">Phylotastic!</a></li>

		</ul>
	</div>
	<div id="menubottom"></div>

	
	<div id="content">

		<!-- Normal content: Stuff that's not going to be put in the left or right column. -->
		<div id="normalcontent">
			<div class="contentarea">
				<!-- Normal content area start -->

<p>Use the syntax <code>Genus_species,Genus_species</code>: comma-delimited, no spaces, underscores in names</p>
 <p><form action = "http://datelife.org/cgi-bin/R/result" method="get">
 <p>Taxa: <input type="text" name="taxa" value="Rhinoceros_unicornis,Equus_caballus" size="80"></p>

 <p>Return format: <select name="format">
 <option value="html">HTML</option>
 <option value="newick">Newick (in the future)</option>
 <option value="bestguess">The single best guess (in the future)</option>
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
				<!-- Normal content area end -->
			</div>
		</div>
	</div>

	<div id="footer">
			<div class="left">Reuse encouraged</div>
			<div class="right">Design by <a href="http://www.nodethirtythree.com/">NodeThirtyThree Design</a></div>
	</div>
	
</div>

</body>
</html>
