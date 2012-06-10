<html>
<?php include("pagestart.html"); ?>
<h2>Contact</h2>
<p>The best way to contact us is via the datelife-users Google group.</p>
<br />
<iframe id="forum_embed"
  src="javascript:void(0)"
  scrolling="no"
  frameborder="0"
  width="900"
  height="700">
</iframe>
<script type="text/javascript">
  document.getElementById('forum_embed').src =
     'https://groups.google.com/forum/embed/?place=forum/datelife-users'
     + '&showsearch=true&showpopout=true&showtabs=false'
     + '&parenturl=' + encodeURIComponent(window.location.href);
</script><?php include("pageend.html"); ?>
