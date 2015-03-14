<div>
<table cellpadding='0' width='930' border='0' cellspacing='0'>
<tr>
<td align='center' height='50'>
<img src='http://graphics.ewha.ac.kr/C2A/figs/models.jpg' border='0' width='938' height='134'>
<p align='center'><font face='verdana'><font face='verdana'>
<h1>C<sup>2</sup>A: Controlled Conservative Advancement for Continuous Collision Detection of Polygonal Models</h1>
<br><br>
<font face='Times New Roman'><font size='3'>
<h2>Min Tang<sup>1</sup>,  Young J. Kim<sup>1</sup> and Dinesh Manocha<sup>2</sup></h2>

<sup>1</sup>Department of Computer Science & Engineering<br>
<br>
Ewha Womans University, Seoul, Korea<br>
<br>
<a href='mailto:tangmin@ewha.ac.kr'>tangmin@ewha.ac.kr</a>  <a href='mailto:kimy@ewha.ac.kr'>kimy@ewha.ac.kr</a>

<sup>2</sup>Department of Computer Science<br>
<br>
University of North Carolina at Chapel Hill<br>
<br>
<a href='mailto:dm@cs.unc.edu'>dm@cs.unc.edu</a>
<br><br>
<h2>Download</h2>

1. Publication in ICRA 2009 in <a href='http://graphics.ewha.ac.kr/C2A/C2A.pdf'>PDF</a>  (698KBytes), <a href='http://graphics.ewha.ac.kr/C2A/bibtex_c2a.txt'>BibTex</a>

<b>To Appear in the IEEE International Conference on Robotics and Automation,  May 12 - 17, Japan, 2009.</b>



2. Source Code in <a href='http://code.google.com/p/c2a-ewha'>C++</a><img src='http://graphics.ewha.ac.kr/C2A/new.gif' />
<br><br>


<p align='left'>
<h2>Abstract</h2>

<blockquote>We present a simple and fast algorithm to perform continuous collision detection between polygonal models undergoing rigid motion for interactive applications. Our approach can handle all triangulated models and makes no assumption about the underlying geometry and topology. The algorithm uses the notion of conservative advancement (CA), originally developed for convex polytopes. We extend this formulation to general models using swept sphere volume hierarchy and present a compact formulation to compute the motion bounds along with a novel controlling scheme. We have implemented the algorithm and highlight its performance on various benchmarks. In practice, our algorithm can perform continuous collision queries in few milli-seconds on models composed of tens of thousands of triangles.<br>
<br>
<h2>Benchmarking Scenarios</h2></blockquote>

<h3>1. Club vs Club (104.8K)</h3>
<blockquote>In the figure below, the red, blue, and green objects denote the club model at initial, final and TOC configurations respectively. The TOC configuration is computed by our continuous collision detection algorithm after linearly interpolating the initial and final configurations. The yellow object shows the configuration of an obstacle. There are more than 200 simulation steps in the benchmark, and in all trials no collision-free situation exists.</blockquote>

<img src='http://graphics.ewha.ac.kr/C2A/figs/Club.jpg' />

<img src='http://graphics.ewha.ac.kr/C2A/figs/timing%20of%20Club.jpg' />

<ul><li>Club vs Club<b>Download Video in WMV <a href='http://graphics.ewha.ac.kr/C2A/Videos/club.wmv'><img src='http://graphics.ewha.ac.kr/C2A/figs/WM.gif' /></a> (0.7MBytes)</b></li></ul>

<h3>2. Hammer (1.7K) vs CAD Piece (2.6K)</h3>
<blockquote>A similar set up like the benchmark 1.</blockquote>

<img src='http://graphics.ewha.ac.kr/C2A/figs/hammer.jpg' />

<img src='http://graphics.ewha.ac.kr/C2A/figs/timing%20of%20hammer.jpg' />

<ul><li>Hammer vs CAD Piece <b>Download Video in WMV <a href='http://graphics.ewha.ac.kr/C2A/Videos/hammer&cad.wmv'><img src='http://graphics.ewha.ac.kr/C2A/figs/WM.gif' /></a>(0.4MBytes)</b></li></ul>

<h3>3.Rigid Body Dynamics for Bunnies</h3>

<blockquote>Using the same benchmarking (<b>Rigid Body Dynamics for Bunnies</b>) setup as FAST: <a href='http://graphics.ewha.ac.kr/FAST/'>http://graphics.ewha.ac.kr/FAST/</a>
In the figure, the TOC configurations are shown as green objects. In the following graph, "<code>*</code>" shows the simulation steps when the TOC should be computed.</blockquote>

<img src='http://graphics.ewha.ac.kr/C2A/figs/dybunny.jpg' />

<img src='http://graphics.ewha.ac.kr/C2A/figs/timing%20of%20dybunny.jpg' />

<ul><li>Rigid Body Dynamics for Bunnies <b>Download Video in WMV  <a href='http://graphics.ewha.ac.kr/C2A/Videos/dybunny.wmv'><img src='http://graphics.ewha.ac.kr/C2A/figs/WM.gif' /></a>(0.5MBytes)</b>
<br></li></ul>

<h2>FAQ</h2>

<p align='left'>A: C2A_Model::AddTri<code>(</code>const PQP_REAL <code>*</code>p1, const PQP_REAL <code>*</code>p2, const PQP_REAL <code>*</code>p3, int id) seems have problems. It does not set the index member in C2A_Tri, which is used later.</p>
<p align='left'>Q: In C2A_Model class, we have the function AddTri(const PQP_REAL <code>*</code>p1, const PQP_REAL <code>*</code>p2, const PQP_REAL <code>*</code>p3, int i, int i1, int i2, int i3). The AddTri function is overloaded. You can see our demo for correctly using.</p>


<p align='left'>A: C2A_QueryContact crashes sometime. The reason is not clear yet. Fortunately I do not need to compute contact point, only contact time is needed.</p>
<p align='left'>Q: If you don't need contact points, you can comment C2A_QueryContact in C2A_Solve function, you will save much timing.</p>


<p align='left'>A: The enum name in enum C2A_Result{OK, TOCFound, CollisionFound, CollisionFree, CollisionNotFound}; seems to be somewhat confusing. Like CollisionNotFound is in fact colliding! </p>
<p align='left'>Q: When two objects are intersected at beginning, C2A can not be used to 	find out collision in the case. The CollisionNotFound is returned, which means collision is happened at beginning, and at same time a message "collision at time 0" will be printf in the screen. As shown in demo, we use the bool value collisionfree in result class C2A_Result to represent in-collision or collision-free. Yes, you are right. The name of labels is somewhat confuse.</p>


<h2>RELATED LINKS</h2>
<h3>Fast</h3>
<blockquote><a href='http://graphics.ewha.ac.kr/FAST/'>http://graphics.ewha.ac.kr/FAST/</a></blockquote>

<h3>CATCH</h3>
<blockquote><a href='http://graphics.ewha.ac.kr/CATCH'>http://graphics.ewha.ac.kr/CATCH</a></blockquote>

<h3>PQP</h3>
<blockquote><a href='http://www.cs.unc.edu/%7Egeom/SSV/index.html'>http://www.cs.unc.edu/%7Egeom/SSV/index.html</a>
<font size='2'>
<br><br><br><br>
Copyright 2009 Computer Graphics Laboratory<br>
Dept of Computer Science & Engineering<br>
Ewha Womans University, Seoul, Korea