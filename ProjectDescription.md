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


For the Final Project, you will be studying trends in the surface segregation energetics of transition metals doped in a host perovskite oxide (SrTiO<sub>3</sub>), and delineate the facet and strain effects on the segregation behavior. The end goal is to use the segregation trends in gaining mechanistic insights into the process of Exsolution. You will be working in groups of two or three on a unique transition metal-SrTiO<sub>3</sub> combination, performing calculations individually, but complementing one another. Each group will present their results in class and be critiqued by their peers. Finally, each group will jointly write a final report, due on <font color="red">May 1 at 5:00 PM (hard deadline)</font>.

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


Current research has foused on trying to understand the mechanism behind the exsolution process, and also develop broader trends in identifying suitable transition metal-host perovksite combinations. Some of the proposed mechanisms include sub-surface nucleation of metal nanoparticles followed by segregation ([Tae-Sik Oh et.al. (2015)](https://pubs.acs.org/doi/full/10.1021/acs.jpclett.5b02292)) or surface segregation followed by sintering and particle growth ([Katz et.al. (2011)](https://pubs.acs.org/doi/abs/10.1021/ja2082284)). 

In this project, you will be exploring the segregation-first mechanism to understand the facet and strain dependence of surface segregation for different transition metals in SrTiO<sub>3</sub>.

Your goals for the project will be to: 
(1) Explore the facet dependence for the thermodynamics of surface segregation,
(2) Explore if the thermodynamics can be changed by applying a tensile strain to the host lattice, and,
(3) Identify suitable electronic structure descriptors for the segregation process.


<a name='deadlines'></a>

## Deadlines ##
1. Show the facet dependence of surface segregation for a single transition metal dopant in the host perovskite; in the absence of any strain (Single plot on 1 Power Point slide): <font color="red"> Wed Apr 11 </font> (Each group)
2. Presentation date: <font color="red"> Tue Apr 24 </font> (Each group 20 min including presentation and Q&A. Note: Total 80 mins)
3. Final written report: <font color='red'> Tue 1 May at 5:00 PM </font>

<a name='calcs'></a>

## Calculations ##

**IMPORTANT:**

Before starting any calculations for the project, please make note of the naming convention to be used (This is to be strictly followed to have your projects properly graded).

All calculations will be run in their corresponding directories within `$WORK/CBE544/` following the naming convention outlined below:

```bash
$WORK/CBE544/FinalProj_M-STO/no-strain/sub-surf/001-AO/
$WORK/CBE544/FinalProj_M-STO/no-strain/sub-surf/001-BO2/
$WORK/CBE544/FinalProj_M-STO/no-strain/sub-surf/110/
$WORK/CBE544/FinalProj_M-STO/no-strain/sub-surf/111/
$WORK/CBE544/FinalProj_M-STO/no-strain/surf/001-AO/
$WORK/CBE544/FinalProj_M-STO/no-strain/surf/001-BO2/
$WORK/CBE544/FinalProj_M-STO/no-strain/surf/110/
$WORK/CBE544/FinalProj_M-STO/no-strain/surf/111/
...
$WORK/CBE544/FinalProj_M-STO/tensile/clean/001-AO/
$WORK/CBE544/FinalProj_M-STO/tensile/clean/001-BO2/
$WORK/CBE544/FinalProj_M-STO/tensile/clean/110/
$WORK/CBE544/FinalProj_M-STO/tensile/clean/111/
...
```
Create a `FinalProj_M-STO` folder (M is the transition metal you are assigned, please check [Assignment](https://cbe544.github.io/CBE544-2018.github.io/Project_Assignments/)) in your `CBE544` directory. For example, if you are assigned Pd, please run the following command to create the directories: 

```bash
cd work
cd CBE544
mkdir FinalProj_Pd-STO 
```
Please change Pd in the above command to the one you are assgined. 

**Instructions on continuing/extending a timed-out job:**
Follow these steps to continue a job that has timed out:

1. Make a new directory within the directory you were originally running the job in. Call it `extend`
   
   For example, `mkdir extend` within `$WORK/CBE544/FinalProj_M-STO/no-strain/sub-surf/001-AO/`
   
2. Copy the relevant files needed to the `extend` directory, this will include `relax.py`, `spede_esp.sub` and most importantly the trajectory from the previous run `opt.traj`. (Note: You can either copy the last frame or the full trajectory, since by default ASE reads the last frame)

3. Go to the `extend` directory, i.e. `cd extend` and submit the job. Please make sure you are reading the correct starting configuration before submitting the job. This can be done by changing the file name within `relax.py` to `opt.traj` or renaming `opt.traj` to the file name provided in `relax.py`


**Task 1:**

Your first task in this project is to get the facet dependence of surface segregation for the transition metal you are assigned, in the absence of any strain. Since the low Miller index surfaces are the ones typically exposed, you will be considering the (001), (110) and (111) surfaces of SrTiO<sub>3</sub>. Note that the (001) facet has two surface terminations in a perovskite oxide (ABO<sub>3</sub>), i.e. either AO terminated or BO<sub>2</sub> terminated. In this project, you will be considering both the terminations. 

To model segregation, you will be replacing one 'B-site' atom (in this case one Ti atom) with the transiton metal assigned to you, in (i) the sub-surface layer and (ii) the surface layer for each of the different facets. You will then perform a geometry optimization calculation to get the total energy of the relaxed structure in each case. Once you obtain this, the segregation energy is defined as: ∆*E*<sub>seg</sub> = E<sub>surface</sub>-E<sub>sub-surface</sub>.

For Task 1, we have already provided you with pre-relaxed trajectories of SrTiO<sub>3</sub> in the different facets, which you will use to setup and perform your calculations. For the next task (Task 2) which involves understanding the effects of strain on the segregation behavior, you will develop your own script for generating the bulk-terminated strained surfaces.

To start with Task 1, first, download and unarchive the files you need, via:

```bash
wget https://cbe544.github.io/CBE544-2018.github.io/ASE/Project.tar.gz
tar -zxvf Project.tar.gz
```

This will create a directory named `Project`. Inside, you will find all of the pre-relaxed .traj files of the different clean  SrTiO<sub>3</sub> surfaces, which you will need for this task. They are named sto-'facet'.traj for the corresponding facets. You will also find a `relax.py` script and a submission script `spede_esp.sub` for submitting your calculations. Please use the current `relax.py` and `spede_esp.sub` scripts provided with the Project and not something you used previously.

To setup and get started with the calculations, you will need to modify the pre-relaxed .traj files provided to you by replacing one Ti atom of the SrTiO<sub>3</sub> system either in (i) the sub-surface layer or (ii) the surface layer for each of the different facets, with the transition metal assigned to you, using the GUI or by writing a python script with ASE commands. As an example, take a look at the 001-BO2 and 001-AO surfaces shown below:

<center><img src="../Images/BO2-1.png" alt="bo2" style="width: 400px;"/> 
<center><img src="../Images/AO-1.png" alt="ao" style="width: 400px;"/>
<br>Setting up your starting configurations (eg. top: 001-BO2 and bottom: 001-AO)</center></center>


<font color="red"> Wed Apr 11 </font> Once you complete all the optimization calculations, write a simple Python script to generate a plot of ∆*E*<sub>seg</sub> vs Facets. Report this plot on a single Power Point slide. Hint: The [Pyplot module](https://matplotlib.org/api/pyplot_summary.html) may come in handy when writing the plotting script.


**Task 2:**

For Task 2 you will be studying the facet dependence of segregation in the presence of an in-plane tensile strain applied to the host perovskite. The host perovskites themselves are typically grown on a support with a different lattice constant, which results in a lattice mismatch. This usually induces an in-plane strain to the host perovskite. Your task is to determine the effect of this strain on the thermodynamics of segregation.  

First, you will write a python script to generate the different bulk-terminated strained surfaces (i.e. (001)-AO terminated, (001)-BO2 terminated, (110) and (111)) of SrTiO<sub>3</sub>. In Task-2 you will consider a SrTiO<sub>3</sub> lattice with a lattice constant of 3.990 Å, which corresponds to a +1.4% in-plane tensile strain applied to the lattice. To start with this, first download the bulk structure of SrTiO<sub>3</sub> using:
```bash
wget https://cbe544.github.io/CBE544-2018.github.io/ASE/SrTiO3.traj 
```
 Please use this structure (which is a slightly different rendition) and not anything else you used previously in HW-5, to be consistent with the surface structures used in Task 1. Use this, and the [build.surface module](https://wiki.fysik.dtu.dk/ase/ase/build/surface.html) (check out the section on Create specific non-common surfaces) to develop your script. Note that you will have to adjust the lattice constant of the bulk structure to the value given above. Also, please make sure your structure has 4 layers along the z-axis, and is centered with a vacuum space of 10 Å, when writing the script. 
Run this script in the login node using for example `python xyz.py`. This should give you a 1x1 cell with 4 layers that is centered with a 10 Å vacuum space. Use the GUI to repeat this structure in the x-y plane to get a 2x2x4 super-cell slab. Make sure to save the super-cell. Finally, since the surface is always in equilibrium with the bulk, apply constraints to all the atoms in the top/bottom 2 layers of the slab to mimic the bulk (this will vary depending on the termination and you can use the structures in Task-1 to get an idea). To do this you may use a python script or the GUI. If you are using the GUI use the `Constraints` option in the `Tools` tab. If you are writing a script, please make use of the [FixAtoms class](https://wiki.fysik.dtu.dk/ase/ase/constraints.html). Save your final structure with the applied constraints.

Once you have the bulk-terminated strained surfaces of SrTiO<sub>3</sub>, you will first perform a geometry optimization of the plain SrTiO<sub>3</sub> surfaces. Make sure to follow the naming convention described earlier when setting up this calculation. You may make use of the `relax.py` and `spede_esp.sub` scripts you used in Task-1 for this calculation. After completing the geometry optimization of the strained SrTiO<sub>3</sub> surfaces, use the optimized structures to repeat the calculations you did for Task-1, by modifying the struture with the transition metal assigned and computing the segregation energy for the different facets.

Plot the segregation energies vs the different facets as in Task-1, to delinate the effect of strain on the segregation thermodynamics.

**Task 3:**

For Task-3, you will be performing some basic analysis involving the electronic structure, with the aim of rationalizing the trend in the segregation energetics computed in **Task 1**. Specifically, you will be tasked with computing the density of states of the dopant transition metal in the surface and sub-surface states for each facet, to explain the specific trend in the energetics observed for the given transition metal dopant.

Before getting started with this task, please provide the optimized trajectories to me following the naming convention below (pool everything per group before sharing with me) :

```bash
M-STO/no-strain/sub-surf/001-AO/
M-STO/no-strain/surf/001-AO/
...
```
To get started with the analysis, first download the following:

```bash
wget https://cbe544.github.io/CBE544-2018.github.io/ASE/dosplot.py 
```
The density of states will be saved to a pickle file called `dos.pickle` for each surface in the `M-STO` directories. Use the ones that are relevant to the transition metal assigned to your group. Please make sure to organize the density of states analysis following the naming convention below:

```bash
$WORK/CBE544/FinalProj_M-STO/no-strain/sub-surf/001-AO/dos/
$WORK/CBE544/FinalProj_M-STO/no-strain/surf/001-AO/dos/
...
```

The `dosplot.py` is a generic script that reads a `dos.pickle` file and plots the projected density of states of the atom specified. You will have to append this file to read multiple `dos.pickle` files and plot the corresponding density of states of the dopant atoms, in the sub-surface and surface configurations for each facet. You can run this script on the login node directly using `python dosplot.py`

To summarize, your final plot for each facet should have the following:
For eg. consider the (110) surface: Density of states of the dopant transition metal in (110) sub-surface and (110) surface. (**All for the no-strain case**)

Finally, in your presentations please also inlude a structural analysis, i.e. compare the relaxed structure with the inital unrelaxed bulk-terminated structure to identify any potential reconstructions induced by the transition metal dopant.

<a name='analysis'></a>

**Archive Calculations**

Please provide the optimized trajectories of the doped systems with tensile strain in a similar format to what was done for Task-3. (Pool everything together per group). Follow the naming convention below:

```bash
M-STO/tensile/sub-surf/001-AO/
M-STO/tensile/surf/001-AO/
...
```
Once you have the directory, zip it and share it over mail (email: sanraman@seas.upenn.edu)

 **Analysis**

Your analysis and final report should include the following:

1. Atomic structures of the different facets of  ‘M’-SrTiO3; where ‘M’ is the transition metal assigned to your group.
2. Tabulated surface segregation energy for the different facets in the absence strain and in the presence of tensile strain.
3. Surface segregation energy vs Facet plot; in the absence and presence of strain.
4. Potential electronic structure quantities to be considered as descriptors: Work function and/or Density of States. Discussion of how and if these suggested parameters couple to the calculated segregation energies.


<a name='report'></a>

## Final Report ##

The final report should be in the form of a 3-5 pages long mini paper including figures and tables. Provide one report for each group. Please be succinct and organize it in the following way:

* Introduction (brief) - don't write too much
* Calculation details
* Results and discussion including analysis.
* Conclusion (brief)
