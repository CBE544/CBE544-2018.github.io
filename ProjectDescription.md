---
layout: page
mathjax: true
permalink: /Project/
---

## Course Project ##
1. [Introduction](#intro)
2. [Deadlines](#deadlines)
3. [Calculations](#calcs)
4. [Analysis](#analysis)
5. [Final Report](#report)


For the Final Project, you will be studying trends in the surface segregation energetics of transition metals doped in a host perovskite oxide (SrTiO<sub>3), and delineate the facet and strain effects on the segregation behavior. The end goal is to use the segregation trends in gaining mechanistic insights into the process of Exsolution. You will be working in groups of two or three on a unique transition metal-<SrTiO<sub>3> combination, performing calculations individually, but complementing one another. Each group will present their results in class and be critiqued by their peers. Finally, each group will jointly write a final report, due on <font color="red"> May 1 at 5:00 PM (hard deadline)</font>.

Please make use of the [Piazza](https://piazza.com/) page for troubleshooting, discussions and for sharing results.

Turn in your final report by emailing a PDF file to:

```
alevoj@seas.upenn.edu, sanraman@seas.upenn.edu
```

<a name='intro'></a>

## Introduction ##

Supported metal catalysts have been used in applications such as catalytic converters in automotive exhaust systems, but pose the issue of particle coarsening/growth, coking etc. which reduce their active lifetime. To circumvent this problem, resarchers at [Daihatsu and Toyota](https://www.nature.com/articles/nature00893) proposed the idea of exsolving metal nanoparticles from a host perovskite, and reversibly re-dispersing it back into the host lattice under reducing and oxidizing conditions respectively. This synthesis method involving exsolution has emerged as a new promising alternative to overcome some of the limitations posed by conventional supported metal catalysts and is an active area of research.

To illustrate the idea of a regenerative 'intelligent' catalytic system, a simple schematic diagram is shown below:

<center><img src="../Images/Regen-catalyst-image.jpg" alt="regen" style="width: 450px;"/>
<br> Self-regenerating 'intelligent' catalyst (<a href="https://www.sciencedirect.com/science/article/pii/S0920586106002719">Tanaka et. al. (2006)</a>)</center>

Due to the high operating pressures and temperatures required for this reaction, alternative catalysts are still needed for this process. [Medford et. al. (2015)](http://dx.doi.org/10.1016/j.jcat.2014.12.033) have suggested that the linear scaling between the dissociation energy of N<sub>2</sub> and its transition state energy prevents most catalysts from achieving a high rate. Assuming that the bond-breaking of N<sub>2</sub> is rate limiting, then traditional metal catalysts have a transition state that is too high in energy. This is illustrated in the filled contour plot below, where the turnover frequency is plotted as a function of the transition state energy of the first N<sub>2</sub> bond breaking (*E*<sub>N-N</sub>) and the dissociation energy (∆*E*<sub>diss</sub>). A catalyst would need to behave differently from these extended surfaces in order to land in a more active region of the map. 

<center><img src="../Images/N2_volcano.png" alt="N2 volcano" style="width: 400px;"/>
<br>Filled contour plot for the turnover frequencies (Singh et. al. (2016))</center>

We will be exploring 2D MXenes where such configurations might be found. 

Your goals for the project will be to: 
(1) explore this reaction to find unique adsorption configurations and possibly more favorable thermodynamics,
(2) explore if the theormodynamics can be changed with functionalization of the MXenes,
(3) explore relations between the reaction intermediates in the pathway,
(4) compare the chemsitry of the carbide and nitride MXenes.

<a name='deadlines'></a>

## Deadlines ##
1. Show energy diagram for bare MXene on Wed April 12 on a Power Point slide (Each group)
2. Presentation date Fri Apr 21 (Each group 15 min including presentation and Q&A. Note: Total 90 mins)
3. Hand in a final written report: Mon 1 May at 5:00 PM

<a name='calcs'></a>

## Calculations ##

For the Final Project, create a `FinalProj_M2X` folder (M2X is the material you are assigned, please check [Assignment](https://cbe544.github.io/Project_Assignments/)) in your `CBE544` directory. For example, if you are assignemnt with Mo2C, please run the following command to create the directories: 

```bash
cdw
cd CBE544
mkdir FinalProj_Mo2C 
```
Please change Mo2C in the above command to the one you are assgined. 

You may run the exercises in any directory (as long as it is under `$WORK`), but keep all the final files for the project organized.

To describe the full reaction on your catalytic system, you will need to calculate the adsorption energies of all intermediates, in their most stable configuration (N\*, NH\*, NH<sub>2</sub>\*, NH<sub>3</sub>\*, H\*). A mean field approximation can be used in the analysis (*e.g.* ∆*E*<sub>2NH</sub> = 2∆*E*<sub>NH</sub>). You are not required to calculate any of the transition states for this assignment. Instead use the universal BEP relations for N<sub>2</sub> dissociation.

First, download and unarchive the files you need via:

```bash
wget https://github.com/CBE544/CBE544.github.io/raw/master/Final_Project.tar.gz
tar -zxvf Final_Project.tar.gz
```

This will create a directory named `Class`. Within, you will find pre-relaxed .traj files for the project. Your team will need the bare, O-terminated, and H-terminated MXenes. They are labeled as M2X.traj for the bare, M2XO2.traj for the O-terminated, and M2XH2.traj for the H-terminated; e.g. for Mo<sub>2</sub>C, the .traj file is Mo2C.traj for the bare MXene, Mo2CO2.traj for the O-terminated, and Mo2CH2.traj for the H-terminated.

In summary:

1. Structural relaxations on both the bare MXene and the two functionalized MXenes; that is, O-terminated and H-terminated. 
2. Adsorption energies for the intermediates in the adsorbed state (N\*, NH\*, NH<sub>2</sub>\*, NH<sub>3</sub>\*, H\*). Check all possible sites in order to determine optimal adsorption configurations. 
3. Energy diagrams for the overall reaction.
<!--4. Calculation of the reaction rate and also a free energy diagram with some temperature and pressure dependence. [Project Part 3](../ASE/Transition_States)-->

**IMPORTANT:**

When you have finished all your calculations. Confirm that your results are organized in the following way:

```bash
$WORK/CBE544/FinalProj_M2X/bare/
$WORK/CBE544/FinalProj_M2X/bare/cleansur/
$WORK/CBE544/FinalProj_M2X/bare/Adsorption/
$WORK/CBE544/FinalProj_M2X/bare/Adsorption/N/
$WORK/CBE544/FinalProj_M2X/bare/Adsorption/N/config
$WORK/CBE544/FinalProj_M2X/bare/Adsorption/NH/
$WORK/CBE544/FinalProj_M2X/bare/Adsorption/NH/config
...
$WORK/CBE544/FinalProj_M2X/O-term/
$WORK/CBE544/FinalProj_M2X/O-term/cleansur/
$WORK/CBE544/FinalProj_M2X/O-term/Adsorption/
$WORK/CBE544/FinalProj_M2X/O-term/Adsorption/N/
$WORK/CBE544/FinalProj_M2X/O-term/Adsorption/N/config
$WORK/CBE544/FinalProj_M2X/O-term/Adsorption/NH/
$WORK/CBE544/FinalProj_M2X/O-term/Adsorption/NH/config
...
$WORK/CBE544/FinalProj_M2X/H-term/
$WORK/CBE544/FinalProj_M2X/H-term/cleansur/
$WORK/CBE544/FinalProj_M2X/H-term/Adsorption/
$WORK/CBE544/FinalProj_M2X/H-term/Adsorption/N/
$WORK/CBE544/FinalProj_M2X/H-term/Adsorption/N/config
$WORK/CBE544/FinalProj_M2X/H-term/Adsorption/NH/
$WORK/CBE544/FinalProj_M2X/H-term/Adsorption/NH/config
...
```

You should rename `config` to describe the binding configuration, such as `fcc`, `hcc`, `top`, `bridge` sites. You should have one calculation per directory.

<a name='analysis'></a>

## Analysis ##

Your analysis and final report should include the following:

1. Structures of the different facets of  ‘M’-SrTiO3; where ‘M’ is the transition metal assigned to your group.
2. Tabulated surface segregation energy for the different facets in the absence and presence of strain.
3. Surface segregation energy vs Facet plot; in the absence and presence of strain
4. Electronic structure descriptors; Work function and Density of States


<a name='report'></a>

## Final Report ##

The final report should be in the form of a 3-5 pages long mini paper including figures and tables. One report for each group. Please be succinct and organize it in the following way:

* Introduction (brief) - don't write too much
* Calculation details
* Results and discussion
* Conclusion (brief)
