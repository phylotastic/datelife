<html>
<?php include("pagestart.html"); ?>
<h2>FAQ</h2>
<br />
<p><b>Don't you know about <a href="http://www.timetree.org">TimeTree</a>? Why do a new site?</b></p>
<p><a href="http://www.timetree.org"a>TimeTree</a> is a great resource. It has expert-derived dates as well as dates from published trees. It also has a good web interface (for mobile devices, too) and now free <a href="http://itunes.apple.com/us/app/timetree/id372842500?mt=8">iPhone</a> and <a href="http://itunes.apple.com/us/app/timetree-hd/id468820058?mt=8">iPad</a> interfaces. You can even get free book chapters and poster images or buy a physical book or poster. What you may not do is reuse the data in an automated manner (at least when we started our project in June 2012; hopefully this FAQ is now out of date on this point). For example, <a href="http://www.timetree.org"a>TimeTree</a> displays this note (as of June 8, 2012):
<br /><center><img src="images/timetreereuse.png"></center>
<br />and the FAQ states (as of June 8, 2012):
<br /><center><img src="images/timetreedownload.png"></center>
<br />Requests to allow download of the calibrated tree or other widespread reuse have been consistently but politely declined. Reusing existing work makes a great deal of sense, but since this was not possible, we decided to develop our own. It differs in other ways from TimeTree in ways that we hope are improvements, but the large-scale reuse was really the key issue driving this development.</p>
<p><b>What can I do with DateLife?</b></p>
<p>You can interact with it as a website to get estimates for the ages of various groups. You can also run automated queries against it. While still in alpha form, the api may change, but you can see what it is at any point by just looking at the URL. For example, searching for a Newick tree for three particular taxa, using only peer-reviewed studies, sends you to this page:
<br /><code><a href="http://datelife.org/cgi-bin/R/result?taxa=Rhinoceros_unicornis%2CEquus_caballus%2CMus_musculus&format=newickmed&partial=liberal&useembargoed=no&uncertainty=100">http://datelife.org/cgi-bin/R/result?taxa=Rhinoceros_unicornis%2CEquus_caballus%2CMus_musculus&format=newickmed&partial=liberal&useembargoed=no&uncertainty=100</a></code>
<br />you could have a python, perl, or R script call this same URL to get trees back. You can download the source code and run it locally, perhaps hacking it for some other purpose. If you reuse in that way, we'd appreciate being cited (citation being the currency of science). If you have an improvement to make, either send us a patch or sign up as a developer. The main code is R, which is fairly easy to use.</p>
<p><b>Who is "we"?</b></p>
<p>Currently, Brian O'Meara, Luke Harmon, Jonathan Eastman, Peter Midford, Tracy Heath, Joseph Brown, Matt Pennell, Mike Alfaro. Maybe you, too?</p>
<?php include("pageend.html"); ?>
