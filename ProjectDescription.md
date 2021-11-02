---
layout: page
mathjax: true
permalink: /Project/
---

## Course Project 
1. [Introduction](#intro)
2. [Deadlines](#deadlines)
3. [Calculations](#calcs)
4. [Analysis](#analysis)<br>
	4a. [LiCoO<sub>2</sub> and Al-doped Data](../CompData)<br>
	4b. [LiCoO<sub>2</sub> magnetism and strucutre](../MagStruct)
5. [Final Report](#report)


For the Final Project, you will be studying the relationship between a materials reactivity and its stability. The students will work in groups of four. Each group will probe a specific material class and each student will be probing a different material. Each group of students will present their results in class that will be critiqued by the other groups. Finally, each group  will jointly write a final report on the combined data. The due date for the final written report is <font color="red">12/16 at 5:00 PM (hard deadline)</font>.

Please make use of the [Piazza](https://piazza.com/) page for troubleshooting, discussions and for sharing results.

Turn in your final report by emailing a PDF file to:

```
alevoj@seas.upenn.edu, csl1191@seas.upenn.edu
```
<a name='intro'></a>

## Introduction ##

Goal: Determine whether a material's reactivity and stability are linked by testing for correlations between easily calculated descriptors.

Plan: Use DFT to calculate oxygen adsorption energies and defect formation energies on surfaces of transition metals, semiconducting oxides, and metallic oxides.

### Motivation ###

Numerous theoretical studies are performed each year that predict materials that will succeed at carrying out crucial reactions (oxygen evolution, nitrogen reduction, CO<sub>2</sub> reduction, etc.) with high rates and good selectivity. Often these materials are then studied experimentally only to find out that the materials perform on par with others that have been tried before. One reason for this often disappointing result is that many theoretical studies still do not consider stability when they are predicting material activity. It is not a trivial property to study and especially when reaction conditions involve a solvent it can be incredibly difficult to accurately model the interface and understand the mechanisms behind a material's activity and stability.

One hypothesis in the field posits that there exists an 'activity-stability' conundrum limiting material performance under reaction conditions. Similar to the Sabatier Volcano in catalysis where an optimal binding strength exists for reactants and products so that they can both adsorb to the surface strongly enough to react but weaklky enough to desorb as products, this conundrum posits that there exists an optimal stability for a material under reaction conditions such that the material is unstable enough to interact with reactants and products to do catalysis, but not so unstable that it's reactive surfaces are short-lived. Likewise if the material is too stable then the surfaces are likely unreactive and reactants are unlikely to interact strongly enough to undergo bond breaking and forming.

We aim to determine whether 1) we can correlate a material reactivity descriptor with a material stability descriptor, and 2) whether we can relate the activity-stability descriptors to find a 1-D relationship equivalent to a Sabatier descriptor to determine a material's optimal position on the activity-stability scale.

## Calculations ##

To get started we need to get the files needed for the final project. Move into your CBE544 directory and run these commands:

```bash
cd
cd CBE544
wget https://upenncbe544.github.io/CBE544-2021.github.io/FinalProject.tar.gz
tar -zxvf FinalProject.tar.gz

```

In the FinalProject directory you should see a directory called scripts, which contains the scripts you will need for this project. You will also see three directories: rutiles, perovskites, and metals. Please check [Assignment](https://cbe544.github.io/Project_Assignments/) in order to determine which material group you will be working with. The bulks for all of these materials have already been optimized using the calculator settings that we determined from HW 5 for the perovskite SrTiO<sub>3</sub>.

### Task 1 ###

You will build surfaces from the optimized bulks for your assigned materials and run the relaxation. Please skip ahead to the section heading relevant to your assigned material.

Perovskites

Build three surfaces: (001)-AO terminated and BO<sub>2</sub>-terminated surfaces and a (110) surface. The number of layers should be 4 and the final trajectory for the surfaces should resemble the ones below.

Rutile Oxides

Build two surfaces: (110) and (111) surface. The number of layers should be 4 and the final trajectory for the surfaces should resemble the ones below.

Metals

You will not be required to build any surfaces as these have already been built and relaxed in our group.

### Task 2 ###

Using the relaxed surfaces you will adsorb Hydrogen and Oxygen species on several unique sites and at a high and low coverage limit. Please skip ahead to the section heading relevant to your assigned material.

Perovskites

On the (001)-AO terminated surface, adsorb H on the O and M sites at coverages of 1 adsorbate per surface and 1 adsorbate per site. Relax these systems. See the figures below for a guide to the adsorption sites.

On the (001)-BO<sub>2</sub>-terminated surface, adsorb H on the O and M sites at coverages of 1 adsorbate per surface and 1 adsorbate per site. Relax these systems. See the figures below for a guide to the adsorption sites.

On the (110) surface, adsorb H on the O and M sites at coverages of 1 adsorbate per surface and 1 adsorbate per site. Relax these systems. See the figures below for a guide to the adsorption sites.

Rutile Oxides

On the (110) surface, adsorb H on the O and M sites at coverages of 1 adsorbate per surface and 1 adsorbate per site. Relax these systems. See the figures below for a guide to the adsorption sites.

On the (111) surface, adsorb H on the O and M sites at coverages of 1 adsorbate per surface and 1 adsorbate per site. Relax these systems. See the figures below for a guide to the adsorption sites.

Metals

On the (111) surface, adsorb H on each site at coverages of 1 adsorbate per surface cell and 1 adsorbate per site. 

On the (210) surface, adsorb H on each site at coverages of 1 adsorbate per surface cell and 1 adsorbate per site. 

On the (100) surface, adsorb H on each site at coverages of 1 adsorbate per surface cell and 1 adsorbate per site. 

### Task 3 ###

The final task will be creating surface defects and relaxing the systems. From these results, we will be able to calculate the defect formation energies. Please skip ahead to the section heading relevant to your assigned material.

Perovskites

On the (001)-AO terminated surface, create a surface layer A-site defect, O-site defect, and AO-defect. Also create a full surface O site defect (remove all surface layer O atoms), full surface A site defect. 

On the (001)-BO<sub>2</sub>-terminated surface, create a surface layer B-site defect, O-site defect, and BO<sub>2</sub>-defect. Also create a full surface O site defect (remove all surface layer O atoms) and full surface B site defect.

On the (111) surface, create a 

Rutile Oxides

On the (110) surface, create a metal defect, an O defect, and high coverage limit defects

Metals

On the (111) surfaces calculate surface metal defect.

On the (211) surface, calculate each unique surface site metal defect. 

On the (100) surface, calculate a metal defect, an O defect, and the high coverage limit equivalent defects.

### Organization ###

Organization for the project is very important so that the data is accesbile once the class is over. We will structure is to be something like this for the 104 surface

```bash
~/CBE544/FinalProject/104/M-surf/noads/
~/CBE544/FinalProject/104/M-surf/ads1/
~/CBE544/FinalProject/104/M-surf/ads2/
~/CBE544/FinalProject/104/M-surf/ads3/
~/CBE544/FinalProject/104/M-surf/noads/bader
~/CBE544/FinalProject/104/M-surf/ads1/bader
~/CBE544/FinalProject/104/M-surf/ads2/bader
~/CBE544/FinalProject/104/M-surf/ads3/bader
~/CBE544/FinalProject/104/M-sub/noads/
~/CBE544/FinalProject/104/M-sub/ads1/
~/CBE544/FinalProject/104/M-sub/ads2/
~/CBE544/FinalProject/104/M-sub/ads3/
~/CBE544/FinalProject/104/M-sub/noads/bader
~/CBE544/FinalProject/104/M-sub/ads1/bader
~/CBE544/FinalProject/104/M-sub/ads2/bader
~/CBE544/FinalProject/104/M-sub/ads3/bader
```
This is an outline of all of the DFT calculations you will need to do, how to organize the files, and where to run each calculation. For the 001 surface it should be something like this:

```bash
~/CBE544/FinalProject/001/CoTerm/noads/
~/CBE544/FinalProject/001/CoTerm/ads1/
~/CBE544/FinalProject/001/CoTerm/ads2/
~/CBE544/FinalProject/001/CoTerm/ads3/
~/CBE544/FinalProject/001/CoTerm/noads/bader
~/CBE544/FinalProject/001/CoTerm/ads1/bader
~/CBE544/FinalProject/001/CoTerm/ads2/bader
~/CBE544/FinalProject/001/CoTerm/ads3/bader
~/CBE544/FinalProject/001/Literm/noads/
~/CBE544/FinalProject/001/Literm/ads1/
~/CBE544/FinalProject/001/Literm/ads2/
~/CBE544/FinalProject/001/Literm/ads3/
~/CBE544/FinalProject/001/Literm/noads/bader
~/CBE544/FinalProject/001/Literm/ads1/bader
~/CBE544/FinalProject/001/Literm/ads2/bader
~/CBE544/FinalProject/001/Literm/ads3/bader
```

Not everyone will be running on calculations so you only need to have the directories of the calculations you will be running. 

#### Jobs not reaching force convergence ####

Some jobs have been running on chestnut and have not reached the force convergence. This means we must extend the job and continue the calculation. I have written a script that will automatically extend the calculation you will just need to copy this by doing :

```bash
cp /home/antcurto/for/CBE544/extend.sh ~/CBE544/FinalProject/scripts
```
If you have a job that ended with convergence you can uses this script by typing:

```bash
sh ~/CBE544/FinalProject/extend.sh
```

This will make a new directory, copy in some files and submit a new job to be run from where the previous calculation has left off. Please only use this if you need to  


### Task 1: ### 

Once you have accurately completed HW5 you can continue on to the final project. We will use the 104 surface and 001 surface trajectories provided to you in the FinalProject directory (this is only because there are specific starting magnetic strcutrures that we want otherwise the structures that you made could be used) to place a metal dopant on the surface and subsurface (separately, so two total calculations). The locations are shown here as in a top view of the 104 surface. The simplest way to change an atom to the desired dopant is to use ase-gui, click on the atom to change, Edit (or ctrl+Y), and type in the element you want. Be sure to save this new trejactory because ase-gui does not automatically save any changes you make. 

<center><img src="../Images/dopantlocations.png" alt="window" style="width: 800px;"/><br>
104 Dopant Locations
</center>

Once you have substitued the metal dopant you can use the relax.py script in ~/CBE544/FinalProject/scripts to relax this system. Copy over the script and submit the job. Record the final energies which can be used to determine if the preferred dopant location is surface or subsurface.

### Task 2: ### 
Using these models from Task 1 we can now adsorb EC to the three locations (per system) shown below:

<center><img src="../Images/Adsorptionlocations104.png" alt="window" style="width: 800px;"/><br>
Locations for Adsorption on the 104 surface of LiCoO<sub>2</sub>
</center>

Refer to the [Adsorption page](../ASE/Adsorption) for instructions on how to add the EC adsorbate. We will use a different script than the relax.py for adsorption DFT calculations. We will use the script called opt-ads.py which you can find in `FinalProject/scripts`. Please be sure to use this script for these calculations or you may expereince converngence issues. 

### Task 3: ### 

Once you have converged the systems with and without EC we will do a [bader charge analysis](http://theory.cm.utexas.edu/henkelman/code/bader/).  Inside the directory where your calculations were run make a new directory called bader (by doing `mkdir bader`). Copy into this directory  fin.traj and vasp-ase.sub. Rename your fin.traj to init.traj (`mv fin.traj init.traj`). Copy from the `FinalProject/scripts` directory the badercharge.py script. It will look like this:

```python
#!/usr/bin/env python
from ase import Atoms, Atom
from ase.calculators.vasp import Vasp
from ase.io import read,write
import numpy as np

p=read('init.traj')
calc = Vasp(prec='accurate',
            encut=520,
            xc='PBE',
            lreal='Auto',
            kpts=[4,4,1],
            nsw = 0,
            ibrion = -1,
            ispin = 2,
            amix_mag = 0.800000,
            bmix = 0.000100,
            bmix_mag= 0.000100,
            amix = 0.20000,
            sigma = 0.05000,
            ediff = 2.00e-04,
            ediffg = -2.00e-02,
            algo ='fast',
            ismear = -5,
            nelm = 250,
            ncore = 16,
            lasph= True,
            ldautype = 2,
            lmaxmix = 4,
            lorbit = 11,
            ldau = True,
            ldauprint = 2,
            ldau_luj={'Co':{'L':2, 'U':3.32, 'J':0},
                      'Li':{'L':-1, 'U':0.0, 'J':0.0},
                      'O':{'L':-1, 'U':0.0, 'J':0.0}
                      },
            lvtot = False,
            lwave = False,
            lcharg = True,
	    laechg= True,
	    gamma=True,
)
calc.calculation_required = lambda x, y: True
p.set_calculator(calc)
pe=p.get_potential_energy()
#####
ana =  Vasp(restart=True)
pend = ana.get_atoms()

forces=pend.get_forces().ravel()
max_force=max([abs(x) for x in forces])

pe = pend.get_potential_energy()
#mag = pend.get_magnetic_moments()

#pend.set_initial_magnetic_moments(mag)
#print mag
write('fin.traj',pend)
```
This script does a static calculation (nsw=0) of the final trajectory from your previous relaxation and writes the files needed to do a bader charge anaylsis. Use the vasp-ase.sub script to submit the badercharge.py script (`sbatch vasp-ase.sub` with the final line `python badercharge.py`). Once the job has finished you will need an updated bader_get_charge_vasp script. To do this do:

```bash
cp /home/antcurto/for/CBE544/bader_get_charge_vasp ~/CBE544/FinalProject/scripts
```

Next you can attach the bader charge to each atom by typing

```bash
module load ase/3.9.1
python ~/CBE544/FinalProject/scripts/bader_get_charge_vasp
```
This will write a new trajectory file called bader_charge.traj that has attached the bader charge of each atom as a magnetic moment. To see this use ase-gui -> View -> Show Labels -> Magnetic Moments. Analyze how the bader charges differ from each system. It is important to load ase/3.9.1 so that the magnetic moments are readable.

### Task 4: ###

Repeat both Task 1 and Task 2 for the 001 surface. The only difference is instead of surface and subsurface we will use Li-terminated vs CoO termianted. You can use the trajectories in the FinalProject directory to run the substituion calculations but we will be doing the adsorption calculations slightly different than the 104. 

Since the calculation for 001 take a long time you will be provided with two trajectories with EC absorbed to the surface already. From here we will alter or dopant location to see the effect of the dopant location within proximity to the adsorbate. To get the trajectories (which are clean LiCoO<sub>2</sub>) type this:

```bash
cp /home/antcurto/for/CBE544/001-Cotermwithads.traj ~/CBE544/FinalProject  
cp /home/antcurto/for/CBE544/001-Litermwithads.traj ~/CBE544/FinalProject
```
This will give you the trajectories with EC adsorbed. The location of where to substitute Co with your dopant can be seen below. 

<center><img src="../Images/dopantlocation001Co.png" alt="window" style="width: 800px;"/><br>
001 Co terminated Dopant Locations
</center>
<center><img src="../Images/dopantlocation001Li.png" alt="window" style="width: 800px;"/><br>
001 Li terminated Dopant Locations
</center>


Run a bader charge analysis on this system as well. See if there are any clear trends between the two systems through things such as dopant location, charge, etc. [Compare your system to the LiCoO<sub>2</sub> and the Al-doped system shown on this page](../CompData). Look for trends between these systems, your own system, and even those of your classmates (if possible)

<a name='deadlines'></a>

## Deadlines ##
1. HW5 Due: Wed 10 April (Each student)
2. Short update (few slides) on completed calcualtions: Wed 17 April during class (1 per group)
3. Final Presentation: Wed 1 May during class (1 per group)
4. Final Paper: Wed 8 May by 5 PM (1 per group)

## Analysis ##

Do an detailed analysis of your system trying to identifty trends in adsorption due to dopant location, charge, surface, or anything else. Compare your data to that of plain LiCoO<sub>2</sub> and Al-doped LiCoO<sub>2</sub> which can be found [here](../CompData).

### Requirements ###

At a minimum you should accomplish the following:

1. Complete the [HW5](../ASE/Getting_Started).
2. Setup a LiCoO<sub>2</sub> surface (104) and calculate adsorption energies for EC adsorption at three sites for two different metal dopant locations (surface and sub-surface).
3. Do a Bader Charge Analysis on metal doped system and metal doped system w/ EC absorbed (8 total bader charge calculations) and compare to the provided systems without a dopant and with an Al dopant.
4. Repeat this process on the 001 facet. Instead of doing is for surface and subsurface we will do this for Li-terminated and CoO<sub>2</sub> terminated. 
5. Analysis
    1. Does the dopant have a strong preference for surface vs subsurface?
    2. How does the metal dopant affect adsorption vs plain LiCoO<sub>2</sub>? vs Al-doped? How is the adosrption different for surface vs sub-surface? 
    3. Do you notice and trends from the bader charge analysis that may contribute to the change in adsorption? Look at different site and different facet terminations.
6. Report (3~5 pages maximum)

### Presentation ##

The Final Presentation should include a summary of all of the calculations you were able to complete. Present the strucutres (with magnetic moments shown please), adsorption energies, and bader charges. You will present that data and an analysis of this data compared to the provided LiCoO<sub>2</sub> and the Al-doped LiCoO<sub>2</sub> data. Also take note of ay major strucutral changes induced by the dopants. Use the structure and magnetism images shown [here](../MagStruct), the 104-LiCoO2.traj and 001-Literm.traj, and the data [here](../CompData) provided to you for some geometric anaylsis.

Discuss problems you had with convergence, magnetism, or any other problems that you encounctered. 

Each Group will also be required to ask questions and begin a disucssion about another groups work. Here are those assignments:

Ni group will ask Ti group questions.<br>
Ti group will ask Mg group questions.<br>
Mg group will ask Ni group questions.

<a name='report'></a>
### Final Report ###

The final report should be in the form of a 3-5 pages long mini paper including figures and tables. Provide one report for each group. Please be succinct and organize it in the following way:

* Introduction (brief) - don’t write too much
* Calculation details
* Results and discussion including analysis.
* Conclusion (brief)

You are welcome to share data amongst your peers to discuss broader trends. 


