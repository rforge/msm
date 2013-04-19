
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="http://<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="http://r-forge.r-project.org/"><img src="http://<?php echo $themeroot; ?>/imagesrf/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like 

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

 -->

<h1>msm: Multi-state models.</h1>

<p> This is the development site for the <b><tt>msm</tt></b> R package for continuous-time multi-state modelling. </p>

<p>For more information about the <b><tt>msm</tt></b> package and multi-state modelling, please read the paper 
<a href="http://www.jstatsoft.org/v38/i08/">Multi-State Models for Panel Data: The msm Package for R</a>
(Jackson, 2011) from Journal of Statistical Software. </p>

<p>Experienced R and C developers with an interest in multi-state models are welcome to volunteer to contribute to this package.   See the <a href="https://r-forge.r-project.org/scm/viewvc.php/pkg/TODO?view=markup&root=msm">TODO list</a> for some ideas.  The R-Forge <strong>project summary page</strong> can be found <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. 

<h2>Downloads</h2>

<p>The current stable version can be found on <a href="http://CRAN.R-project.org/package=msm">CRAN</a>.</p>

<p>The latest development version can usually be found here on <a href="http://r-forge.r-project.org/R/?group_id=1586">R-Forge</a>.  For changes in the development version see the <a href="https://r-forge.r-project.org/scm/viewvc.php/pkg/inst/NEWS?view=markup&root=msm">NEWS</a> or <a href="https://r-forge.r-project.org/scm/viewvc.php/pkg/ChangeLog?view=markup&root=msm">ChangeLog</a> files in the package.</p>

<p><b><tt>msm</tt></b> is an add-on package for the <a href="http://www.r-project.org/">R</a> statistical software.  Please see the <a href="http://cran.r-project.org/doc/manuals/R-admin.html#Installing-packages">R Installation and Administration Manual</a> for information about installing R packages.</p>

<hr>
<address></address>
<!-- hhmts start -->Last modified: Fri Apr 19 17:28:23 BST 2013 <!-- hhmts end -->
</body> </html>

</body>
</html>
