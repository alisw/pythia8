<html>
<head>
<title>Simple Showers</title>
<link rel="stylesheet" type="text/css" href="pythia.css"/>
<link rel="shortcut icon" href="pythia32.gif"/>
</head>
<body>

<script language=javascript type=text/javascript>
function stopRKey(evt) {
var evt = (evt) ? evt : ((event) ? event : null);
var node = (evt.target) ? evt.target :((evt.srcElement) ? evt.srcElement : null);
if ((evt.keyCode == 13) && (node.type=="text"))
{return false;}
}

document.onkeypress = stopRKey;
</script>
<?php
if($_POST['saved'] == 1) {
if($_POST['filepath'] != "files/") {
echo "<font color='red'>SETTINGS SAVED TO FILE</font><br/><br/>"; }
else {
echo "<font color='red'>NO FILE SELECTED YET.. PLEASE DO SO </font><a href='SaveSettings.php'>HERE</a><br/><br/>"; }
}
?>

<form method='post' action='SimpleShowers.php'>
 
<h2>Simple Showers</h2> 
 
<a name="section0"></a> 
<h3>Overview of recent changes</h3> 
 
PYTHIA comes with a complete parton-shower machinery, but also 
allows external shower programs to be linked in to it, see the 
<?php $filepath = $_GET["filepath"];
echo "<a href='ImplementNewShowers.php?filepath=".$filepath."' target='page'>";?>Implement New Showers</a> 
page. Notably the VINCIA and DIRE codes have been structured 
to make use of this functionality. Currently these codes are 
distributed separately, but the intention is to integrate them 
into the PYTHIA distribution in the future. 
 
<p/> 
Originally the <code>TimeShower</code> and <code>SpaceShower</code> 
classes implemented the default PYTHIA showers, but also acted as 
base classes from which the external showers derived. This has some 
disadvantages, so the two aspects are now split. The 
<code>TimeShower</code> and <code>SpaceShower</code> classes remain 
as simple base classes from which the actual showers are derived. 
The physics code has been moved to the new derived 
<code>SimpleTimeShower</code> and <code>SimpleSpaceShower</code>. 
An external shower that does not use any of the existing shower 
algorithms will therefore work as before, which would be the normal 
case, but alternatively a shower could of course derive from the new 
classes and then reuse relevant code in them. 
 
<p/> 
Settings names have been retained, again for reasons of backwards 
compatibility of user code, e.g. in command files. Thus setting names 
beginning with <code>TimeShower:</code>, <code>SpaceShower:</code>, 
<code>WeakShower:</code> or <code>UncertaintyBands:</code> refer 
uniquely to the current baseline "simple" ones. In the future some 
of them may become common with VINCIA and DIRE, notably the 
uncertainty bands ones, whereas ones specific to those two programs 
will have names that spell it out. 
 
<p/> 
The prepending of "Simple" was a minimalistic choice under the 
circumstances; more fancy names could have been chosen. What it 
refers to is that showers like VINCIA and DIRE aim higher, in 
striving to achieve full NLL accuracy, whereas the Simple ones 
operate in an improved LL approximation. In other respects the 
Simple showers can do more different physics than the other two, 
at least currently. Some examples of the broad approach are 
<ul> 
<li>Matrix elements corrections for the first ("hardest") gluon 
emission in most two-body resonance decays, effectively making the 
FSR in these decays NLO accurate.</li> 
<li>There is no corresponding NLO accuracy for ISR in any processes, 
but several examples where reasonably accurate kinematics spectra 
are available over the full phase space, by input of partial 
higher-order information.</li> 
<li>The default dipole-recoil scheme for FSR can be switched to a 
global-recoil option for the first few emissions, in order to simplify 
matching and merging to higher-order calculations.</li> 
<li>The default global-recoil scheme for ISR can be replaced by a 
dipole-recoil scheme, where the other colour dipole end may be in 
the final state.</li> 
<li>Showers off massive objects, within and beyond the Standard 
Model, including e.g. octet onium states.</li> 
<li>Showers interleaved with multiparton interactions, and set up 
to handle <?php $filepath = $_GET["filepath"];
echo "<a href='ASecondHardProcess.php?filepath=".$filepath."' target='page'>";?>two predefined hard 
interactions</a>.</li> 
<li>QED showers, where photons can be emitted and then branch 
into fermion pairs that shower further.</li> 
<li>Weak radiation of <i>W^+-</i> and <i>Z^0</i> off fermions.</li> 
<li>Radiation also in some hadronic decays.</li> 
<li>Possibility to handle both abelian and nonabelian showers in a 
hidden valley sector, where relevant fully interleaved with normal 
QCD and QED radiation.</li> 
<li>A wide selection of further switches and parameters to vary shower 
assumptions: running of <i>alpha_s</i>, <i>p_Tmin</i> value, 
scale choices, gluon polarization effects, mass effects in 
<i> g &rarr; q qbar</i>, etc.</li> 
</ul> 
 
<a name="section1"></a> 
<h3>Shower components</h3> 
 
The <?php $filepath = $_GET["filepath"];
echo "<a href='MasterSwitches.php?filepath=".$filepath."' target='page'>";?>Master Switches</a> for ISR and FSR 
in general, and a switch for QED radiation in 
<?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDecays.php?filepath=".$filepath."' target='page'>";?>Particle Decays</a> to leptons, 
are intended to be common for all shower programs, where applicable. 
 
<p/> 
The full description of settings in the Simple Shower framework 
is spread across several pages: 
<ul> 
<li>The final-state <?php $filepath = $_GET["filepath"];
echo "<a href='TimelikeShowers.php?filepath=".$filepath."' target='page'>";?>Timelike Showers</a> 
cover all aspects of QCD and QED FSR.</li> 
<li>The initial-state <?php $filepath = $_GET["filepath"];
echo "<a href='SpacelikeShowers.php?filepath=".$filepath."' target='page'>";?>Spacelike Showers</a> 
cover all aspects of QCD and QED ISR.</li> 
<li>While the main switches for weak radiation of <i>W^+-</i> and 
<i>Z^0</i> are found in the two previous FSR and ISR pages, 
there also a few common technical 
<?php $filepath = $_GET["filepath"];
echo "<a href='WeakShowers.php?filepath=".$filepath."' target='page'>";?>Weak Showers</a> settings.</li> 
<li>There is a special framework to produce uncertainty bands from 
<?php $filepath = $_GET["filepath"];
echo "<a href='Variations.php?filepath=".$filepath."' target='page'>";?>Automated Variations</a> of basic 
parameters, such as factorization and renormalization scales, or 
choice of parton distributions.</li> 
<li>The settings for final-state showers in a 
<?php $filepath = $_GET["filepath"];
echo "<a href='HiddenValleyProcesses.php?filepath=".$filepath."' target='page'>";?>Hidden Valleys</a> 
are stored along with the switches for such hard processes.</li> 
<li>There is a wide selection of 
<?php $filepath = $_GET["filepath"];
echo "<a href='MatchingAndMerging.php?filepath=".$filepath."' target='page'>";?>Matching and Merging</a> 
approaches that have been implemented so as to work well with these 
showers.</li> 
<li>Tunes that include ISR and FSR parameters are described on the 
<?php $filepath = $_GET["filepath"];
echo "<a href='Tunes.php?filepath=".$filepath."' target='page'>";?>Tunes</a> page.</li> 
<li>The shower evolution can be interrupted or modified with the 
help of <?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?>User Hooks</a>.</li> 
</ul> 
 
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
