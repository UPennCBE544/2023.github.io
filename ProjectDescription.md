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


For the Final Project, you will be studying the relationship between a material's reactivity and its stability. The students will work in two groups. Each group will probe a specific material class and each student will be probing a different material. Each group of students will present their results in class that will be critiqued by the other groups. Finally, each group  will jointly write a final report on the combined data. The due date for the final written report is <font color="red">12/16 at 5:00 PM (hard deadline)</font>.

Please make use of the [Piazza](https://piazza.com/) page for troubleshooting, discussions and for sharing results.

Turn in your final report by emailing a PDF file to:

```
alevoj@seas.upenn.edu, csl1191@seas.upenn.edu
```
<a name='intro'></a>

## Introduction ##

Goal: Determine whether a material's reactivity and stability are linked by testing for correlations between easily calculated descriptors.

Plan: Use DFT to calculate oxygen adsorption energies and defect formation energies on surfaces of perovskite oxides and rutile oxides.

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

In the FinalProject directory you should see a directory called scripts, which contains the scripts you will need for this project. You will also see two directories: rutiles and perovskites. Please check [Assignment](https://cbe544.github.io/Project_Assignments/) in order to determine which material you will be working with. The bulks for all of these materials have already been optimized using the calculator settings that we determined from HW 5 for the perovskite SrTiO<sub>3</sub>. Move into the directory corresponding to your material.

### Things to keep in mind ###
For all of the calculations that you are running here on out, you will need to submit them using stampede.sub to the cluster. Do not run any of these jobs on the login node. 

Each calculation must be carried out in its own directory. All of the files necessary to carry out the calculation must also be in the directory at the time you submit the job or else it will not work. Usually this means the directory should have an ASE trajectory (init.traj), a relaxation script (relax.py), and a submit script (stampede.sub).

It is critical that you organize your directories consistently so that we can find the data later. See the Organization section below for additional guidance. 

### Task 1 ###

You will build surfaces from the optimized bulks for your assigned materials and run the relaxation. In order to build the surfaces, you will be using the surf_build.py script in the scripts directory. This script utilizes the build module in ASE to cut a surface from the optimized bulk trajectory. It reorients the surface so that the surface is perpendicular to the z-direction in the output file. Please skip ahead to the section heading relevant to your assigned material.

### Perovskites ###

You should already be in the directory corresponding to your material. In this directory make 2 new directories: 110 and 111. Go ahead and build the full directory tree for your material that matches the relevant one in the Organization section below. Please make it exactly like this tree. !!!Note that the example tree below is for SrAgO3 and you will need to rename things to match your assigned material depending on the transition metal in your system (throughout this tutorial write-up you will need to keep this in mind whenever I mention Ag!!! Move into the 110 directory. Move into the clean directory. Copy the relaxed bulk trajectory (opt.traj) from the bulk/relax directory corresponding to your material into your 110/clean directory. Copy the surf_build.py from the scripts directory into this directory as well. Now you will build a 4 layer (110) slab model using the surf_build.py script. Look at the surf_build.py script. It should look like the text below: 

```
#!/usr/bin/env python
import numpy as np

from ase import io
from espresso import espresso
from ase import build
from ase.geometry import get_layers
import sys
from ase.optimize import BFGS

atoms =  io.read('opt.traj') #Optimized bulk FCC
traj=io.Trajectory('init.traj','w')
cell=atoms.get_cell() #Get cell of opt.traj
symb = atoms.get_chemical_symbols()
s1=build.surface(atoms, (1,1,0), 5) #cut a surface normal to the [1,1,0] direction and make it 5 layers deep in the z-direction.
s1.center(vacuum=10, axis=2) #center the cell and add 10 Angstroms of vacuum to both sides of the slab in the z-direction

traj.write(s1)
```

You will only need to change the indices in the ```s1=build.surface(atoms, (1,1,0), 5)``` line to correspond to the surface cut you are making. It is set to (110) by default.

Build the (110) surface by running the surf_build.py script with the command ```python surf_build.py```. The script will cut the surface using the bulk opt.traj that you copied and save the surface slab as init.traj in the same directory. The actual number of layers we want is 4 and the final trajectory for the surfaces should resemble the ones below. In order to get these structures you need to remove the asymmetric ABO<sub>3</sub> atoms from the top and bottom to end up with the desired terminations at 4 total layers. Pay attention to the axes in the images below as you orient yourselves. Open init.traj in the GUI and remove the atoms as needed to create the (110) surface that looks exactly like the one below and save it as init.traj. 

In the 111 directory repeat the same process, but build the (111) surface and make it look exactly like the one below. You should now have an init.traj for both the (110) and the (111) surfaces in their own separate directories.

<center><img src="../Images/perovskites_surfs.png" alt="window" style="width: 800px;"/><br>
Schematic of Perovskite Surfaces
</center>


### Rutile Oxides ###

You should already be in the directory corresponding to your material. In this directory make 2 new directories: 110 and 100. Go ahead and build the full directory tree for your material that matches the relevant one in the Organization section below. Please make it exactly like this tree. !!!Note that the example tree below is for MoO2 and you will need to rename things to match your assigned material depending on the transition metal in your system (throughout this tutorial write-up you will need to keep this in mind whenever I mention Mo!!! Move into the newly created 110 directory. Make a directory called clean and move into it. Copy the relaxed bulk trajectory (opt.traj) from the bulk/relax directory corresponding to your material into your 110/clean directory. Copy the surf_build.py from the scripts directory into this directory as well. Now you will build a 4 layer (110) slab model using the surf_build.py script. Look at the surf_build.py script. It should look like the text below: 

```
#!/usr/bin/env python
import numpy as np

from ase import io
from espresso import espresso
from ase import build
from ase.geometry import get_layers
import sys
from ase.optimize import BFGS

atoms =  io.read('opt.traj') #Optimized bulk FCC
traj=io.Trajectory('init.traj','w')
cell=atoms.get_cell() #Get cell of opt.traj
symb = atoms.get_chemical_symbols()
s1=build.surface(atoms, (1,1,0), 5) #cut a surface normal to the [1,1,0] direction and make it 5 layers deep in the z-direction.
s1.center(vacuum=10, axis=2) #center the cell and add 10 Angstroms of vacuum to both sides of the slab in the z-direction

traj.write(s1)
```

You will only need to change the indices in the ```s1=build.surface(atoms, (1,1,0), 5)``` line to correspond to the surface cut you are making. It is set to (110) by default.

Build the (110) surface by running the surf_build.py script with the command ```python surf_build.py```. The script will cut the surface using the bulk opt.traj that you copied and save the surface slab as init.traj in the same directory. The actual number of layers we want is 4 and the final trajectory for the surfaces should resemble the ones below. In order to get these structures you need to remove the asymmetric MO<sub>2</sub> atoms from the top and bottom to end up with the desired terminations at 4 total layers. Pay attention to the axes in the images below as you orient yourselves. Open init.traj in the GUI and remove the atoms as needed to create the (110) surface that looks exactly like the one below and save it as init.traj. 

In the 100/clean directory repeat the same process but build the (100) surface and make it look exactly like the one below. You should now have an init.traj for both the (110) and the (100) surfaces in their own separate directories.

<center><img src="../Images/rutile_surfs.png" alt="window" style="width: 800px;"/><br>
Schematic of Rutile Oxides Surfaces
</center>


### Everyone ###
Once you have built both surface facets for you material, constrain the bottom half of the atoms to the bulk lattice positions. To do this select the atoms that you want to constrain, click Tools -> Constraints -> Constrain Selected Atoms. The constrained atoms should now have dashed 'X's on them. These should match the images above. Make sure that you are constraining a stoichiometric number of atoms. For the Rutile Oxides this means you should be constraining an integer multiple of MO2 atoms. For the Perovskites you should be constraining an integer multiple of ABO<sub>3</sub> atoms.

Next, we need to make sure we have the appropriate vacuum set up between slabs. We will use 20 Angtroms of vacuum between slabs. To do this, copy the script titled vacuum.py into the directory with your trajectory and look at the script. It reads in a file called 'init.traj', centers it in the unit cell, and adds 10 Angstroms of vacuum on both sides of the slab (axis=2 refers to the z-axis) and then rewrites the file as 'init.traj'.

Finally, we can relax these surfaces to get the initial structures that we will use for all of the adsorption and defect calculations going forward. Copy the relax.py script into each surface folder. Copy the stampede.sub script into each surface folder. Then, from each directory make sure your stampede.sub script is copying and running the relevant files. Then run the relaxation using the command ```sbatch stampede.sub```. 

### Task 2 ###

Using the relaxed surfaces you will adsorb O and OH species on several unique sites and coverages. Please skip ahead to the section heading relevant to your assigned material for detailed instructions.

### Perovskites ###

Step 1: Full coverage O adsorption

On the (110) surface:

Open the relaxed surface trajectory using the ase-gui. We want to repeat the unit cell once in the x and y directions to create a 2x2 surface. This will allow for adsorbate-adsorbate interactions to be more accurately captured and for us to probe different adsorbate concentrations. To do this, click on View -> Repeat and then in the window that opens set the x and y boxes (the first and second boxes) equal to 2. Then click Set unit cell. You should now see the larger unit cell and surface looking like the image below. Save this trajectory as init.traj. You will be using this surface trajectory frequently to create your init.traj files for different adsorbate and defect calculations.

<center><img src="../Images/perovskite_2x2.png" alt="window" style="width: 800px;"/><br>
Perovskite (110) surface repeated to form a 2x2 slab model.
</center>

Now we want to set up and run relaxations for O and OH adsorptions at full coverage. Make and move to a directory called 110/oads/1ML/Ag. Copy the relaxed init.traj to this directory. Open the file with ag init.traj. We know from previous work that the favored adsorption site on these surfaces is above the B site cation. So let's start by setting up an adsorption calculation for O. In order to add an adsorbate, click on the atom that you want to add the adsorbate directly above. Press Ctrl+A or click Edit -> Add Atom and in the window that comes up type O in the top box. In the box below Position you can enter a number corresponding to the number of Angstroms above the highlighted Atom you want to place the adsorbate. For the first oxygen adsorbate try 2.3 A like in the image below. Look at the side view and make sure the adsorbed oxygen is in a similar position to the ones in my system. Repeat the same for all four adsorption sites until you have a system with four O adsorbates located above the B site atoms in the lattice.

<center><img src="../Images/perov_o_ads_full.png" alt="window" style="width: 800px;"/><br>
Process for adsorbing O at 1ML coverage on the (110) surface
</center>

Now in his directory copy the relax.py script and the stampede.sub script from the scripts folder. Submit the job using sbatch stampede.sub. 

Step 2: Full coverage OH adsorption

Move to the directory for ohads/Ag. Copy the init.traj that you generated in the previous step to this directory using the cp command. Open this file using the GUI. Now add an H to the top of each O adsorbate at a position of 1 Angstrom above the O. The final structure should look like the images below.

<center><img src="../Images/perov_oh_ads_full.png" alt="window" style="width: 800px;"/><br>
Process for adsorbing OH at 1ML coverage on the (110) surface
</center>

Copy relax.py and stampede.sub to this directory and submit the job. Now you are getting the hang of this process. 

Step 3: Single adsorbate coverage for O and OH

In two new directories (oads/Ag/0.25ML ohads/Ag/0.25ML) create init.traj files that have only one adsorbate instead of four. I would copy the init.traj files from the above steps into these directories and delete 3 of the adsorbates so that you get a file that looks like the ones below for O and OH respectively.

<center><img src="../Images/perov_ads_25.png" alt="window" style="width: 800px;"/><br>
Schematic of Perovskite (110) surface with 0.25 ML O and OH Coverages
</center>

Submit these relaxations as before.

Step 4: Defects 0.25 ML

Go to the directory clean/. Make a directory clean/vac/. Cd into vac/. Make two directories clean/vac/1ML/ and clean/vac/0.25ML/. Cd into 0.25ML. Make four directories: Sr, SrO, Ag, and AgO2. Copy the clean relaxed init.traj into each of these directories. For the Sr case, you will remove the topmost Sr atom, save the init.traj, and submit the job. For the SrO you will be removing the same atom plus the O atom coordinated to it. Your final structures will resemble those shown below. Relax all of these. For the Ag and AgO2 you will be removing Ag and AgO2 groups as shown highlighted below.

<center><img src="../Images/perovs_110_25ml.png" alt="window" style="width: 800px;"/><br>
Highlighted atoms to delete in order to create defect surfaces on Perovskite (110) surfaces
</center>

Step 5: Defects 1 ML

Move into the clean/vac/1ML/ directory that you made in Step 4 above. Make four directories exactly the same as in Step 4: Sr, SrO, Ag, and AgO2. Repeat the same steps except instead of removing a single atom/unit, you will remove all of the equivalent groups on the surface. Relax these.

Step 6: Repeat everything above for the (111) surface

On the (111) surface, you will repeat all of the same calculations as the (110) done previously (including repeating the relaxed slab model to create a 2x2 surface!). See images below for examples of each setup.

<center><img src="../Images/perov_o_ads_full_111.png" alt="window" style="width: 800px;"/><br>
Process for adsorbing O at 1ML coverage on the (111) surface
</center>

<center><img src="../Images/perov_oh_ads_full_111.png" alt="window" style="width: 800px;"/><br>
Process for adsorbing OH at 1ML coverage on the (110) surface
</center>

<center><img src="../Images/perov_ads_25_111.png" alt="window" style="width: 800px;"/><br>
Schematic of Perovskite (110) surface with 0.25 ML O and OH Coverages
</center>

<center><img src="../Images/perovs_111_25ml.png" alt="window" style="width: 800px;"/><br>
Highlighted atoms to delete in order to create defect surfaces on Perovskite (110) surfaces
</center>

Don't forget to also run the full 1ML vac calculations.

### Rutile Oxides ###

Step 1: O adsorption

On the (110) surface:

Open the relaxed surface trajectory using the ase-gui. We want to repeat the unit cell once in the y direction to create a 1x2 surface. This will allow for adsorbate-adsorbate interactions to be more accurately captured and for us to probe different adsorbate concentrations. To do this, click on View -> Repeat and then in the window that opens set the y box (the second box) equal to 2. Then click Set unit cell. You should now see the larger unit cell and surface looking like the image below. Save this trajectory as init.traj. You will be using this surface trajectory frequently to create your init.traj files for different adsorbate and defect calculations.

<center><img src="../Images/rutile_nx2.png" alt="window" style="width: 800px;"/><br>
Rutile oxide (110) surface repeated in the y to form a 1x2 slab model and rutile oxide (100) surface repeated in the x and y to form a 2x2 slab model.
</center>

Now we want to set up and run relaxations for O and OH adsorptions. It will be useful to define the different adsorption sites for the (110) surface. See the figure below for naming conventions. The surface atoms are cus metal sites and bridging oxygen atoms. 

<center><img src="../Images/rutile_surf_descr_110.png" alt="window" style="width: 800px;"/><br>
Rutile oxide (110) surface repeated in the y to form a 1x2 slab model and rutile oxide (100) surface repeated in the x and y to form a 2x2 slab model.
</center>

Let's start by setting up an adsorption calculation for O on the cus site. In order to add an adsorbate, click on the atom that you want to add the adsorbate directly above. Press Ctrl+A or click Edit -> Add Atom and in the window that comes up type O in the top box. In the box below Position you can enter a number corresponding to the number of Angstroms above the highlighted Atom you want to place the adsorbate. For the first oxygen adsorbate try 2 A like in the image below. Look at the side view and make sure the adsorbed oxygen is in a similar position to the ones in my system. Repeat the same for both cus adsorption sites.

<center><img src="../Images/rutile_o_ads_cus_110.png" alt="window" style="width: 800px;"/><br>
Process for adsorbing O on the (110) surface cus site
</center>

Now in this directory copy the relax.py script and the stampede.sub script here from the scripts folder. Submit the job using sbatch stampede.sub. 

Step 2: OH adsorption

Move to the directory for OHads-cus. Copy the Oads-cus init.traj that you generated in Step 1 to this directory using the cp command. Open this file using the GUI. Now add an H to the top of each O adsorbate at a position of 1 Angstrom above the O. The final structure should look like the images below.

<center><img src="../Images/rutile_oh_ads_cus_110.png" alt="window" style="width: 800px;"/><br>
Process for adsorbing OH on the (110) surface cus site
</center>

Copy relax.py and stampede.sub to this directory and submit the job. Now you are getting the hang of this process. 

Step 3: Complete adsorption calculations
In the two remaining adsorption directories (OHads-cus-br and OHads-br) copy the relaxed, clean init.traj for the (110) surface and create surfaces that are identical to those in the Figure below.

<center><img src="../Images/rutile_110_ads.png" alt="window" style="width: 800px;"/><br>
Schematic of Rutile oxide (110) surface with remaining adsorptions
</center>

Submit these relaxations as before.

Step 4: Defects

Go to the clean/vac/0.5ML directory. For each defect directory copy in the relaxed, clean surface. Delete atoms corresponding to the defect (as highlighted in the Figure below) and then relax the trajectory.

<center><img src="../Images/rutile_110_vacs.png" alt="window" style="width: 800px;"/><br>
Highlighted atoms to delete in order to create defect surfaces on Perovskite (110) surfaces
</center>

Step 5: Repeat everything above for the (100) surface

On the (100) surface, you will run adsorptions and defects similar to the (110) done previously (including repeating the relaxed slab model to create a 2x2 surface!). See images below for examples of each setup. Create and relax trajectories matching those below in the directories corresponding from the organization tree.

<center><img src="../Images/rutile_100_oads.png" alt="window" style="width: 800px;"/><br>
Process for adsorbing O at 1ML coverage on the (111) surface
</center>

<center><img src="../Images/rutile_100_ohads.png" alt="window" style="width: 800px;"/><br>
Process for adsorbing OH at 1ML coverage on the (110) surface
</center>

<center><img src="../Images/rutile_100_vacs.png" alt="window" style="width: 800px;"/><br>
Schematic of Perovskite (110) surface with 0.25 ML O and OH Coverages
</center>

### Organization ###

Perovskites: Organization for the project is very important so that the data is accessible once the class is over. We will structure it to be like this for the Perovskite 110 surface. The same structure will apply for the (111) surface with 110 replaced by 111 instead in the paths below.

```bash
~/CBE544/FinalProject/perovskites/srago3/110/clean/
~/CBE544/FinalProject/perovskites/srago3/110/clean/vac/1ML/Sr
~/CBE544/FinalProject/perovskites/srago3/110/clean/vac/1ML/SrO
~/CBE544/FinalProject/perovskites/srago3/110/clean/vac/1ML/Ag
~/CBE544/FinalProject/perovskites/srago3/110/clean/vac/1ML/AgO2
~/CBE544/FinalProject/perovskites/srago3/110/clean/vac/0.25ML/Sr
~/CBE544/FinalProject/perovskites/srago3/110/clean/vac/0.25ML/SrO
~/CBE544/FinalProject/perovskites/srago3/110/clean/vac/0.25ML/Ag
~/CBE544/FinalProject/perovskites/srago3/110/clean/vac/0.25ML/AgO2
~/CBE544/FinalProject/perovskites/srago3/110/oads/1ML/Ag
~/CBE544/FinalProject/perovskites/srago3/110/oads/0.25ML/Ag
~/CBE544/FinalProject/perovskites/srago3/110/ohads/1ML/Ag
~/CBE544/FinalProject/perovskites/srago3/110/ohads/0.25ML/Ag
```
Rutile Oxides: Organization for the project is very important so that the data is accessible once the class is over. We will structure it to be like this for the Rutile (110) and (100) surfaces.

```bash
~/CBE544/FinalProject/rutiles/moo2/110/clean/
~/CBE544/FinalProject/rutiles/moo2/110/clean/vac/0.5ML/br-m
~/CBE544/FinalProject/rutiles/moo2/110/clean/vac/0.5ML/br-mo
~/CBE544/FinalProject/rutiles/moo2/110/clean/vac/0.5ML/br-mo2
~/CBE544/FinalProject/rutiles/moo2/110/clean/vac/0.5ML/br-o
~/CBE544/FinalProject/rutiles/moo2/110/clean/vac/0.5ML/cus-m
~/CBE544/FinalProject/rutiles/moo2/110/clean/vac/0.5ML/cus-mo
~/CBE544/FinalProject/rutiles/moo2/110/Oads-cus/0.5ML
~/CBE544/FinalProject/rutiles/moo2/110/OHads-cus
~/CBE544/FinalProject/rutiles/moo2/110/OHads-cus-br
~/CBE544/FinalProject/rutiles/moo2/110/OH-br
```

```bash
~/CBE544/FinalProject/rutiles/moo2/100/clean/
~/CBE544/FinalProject/rutiles/moo2/100/clean/vac/1ML/m
~/CBE544/FinalProject/rutiles/moo2/100/clean/vac/1ML/o
~/CBE544/FinalProject/rutiles/moo2/100/clean/vac/1ML/mo
~/CBE544/FinalProject/rutiles/moo2/100/clean/vac/1ML/mo2
~/CBE544/FinalProject/rutiles/moo2/100/oads/1ML/m
~/CBE544/FinalProject/rutiles/moo2/100/ohads/1ML/m
~/CBE544/FinalProject/rutiles/moo2/100/ohads/1ML/o
~/CBE544/FinalProject/rutiles/moo2/100/ohads/1ML/m_o
``` 

#### Jobs not reaching force convergence ####

Some jobs may run for the allotted time on Stampede2 and not reach the force convergence. To check whether your job has reached convergence you can go into the directory containing the output files and run the command ```$ tail opt.log```. This will display the last few lines of the log file. The last column displays the force calculated at that step. If the final line shows a force value less than 0.03, then it means the calculation has converged, cince we set that as the threshhold in relax.py. 

If the job has ended (due to failure, timeout, or some other reason) and this force value is not less than 0.03, it means the job needs to be rerun starting from the final structure. In order to do this, you will make a new directory (from inside the current directory) called 'extend'. 
```$ mkdir extend```
Then you can copy the files into this directory:
```$ cp opt.traj relax.py stampede.sub extend/```
Then move into the directory and rename opt.traj as init.traj:
```$ mv opt.traj init.traj```
Now you have everything that you need to start a new calculation from the output of the previous one, and you can submit with:
```$ sbatch stampede.sub```

###Jobs not being returned from $SCRATCH###

Abiding by the rules for running on Stampede2, we have to run the IO intensive jobs from the $SCRATCH partition. This is covered on this website where I introduced you to the submission script. Basically when the job begins, eveything that is needed to run it is copied to the $SCRATCH directory with the name of the jobid, then the job runs, then the files are moved back to the submission directory and the $SCRATCH directory is deleted. Sometimes, the job fails prematurely and the stampede.sub script is unable to complete resulting in files getting stranded in the $SCRATCH directory. If this happens, you will receive an email that the job failed, but when you go to the submission directory you will not see any output files (e.g. opt.traj, calcdir, etc...). The files are simply stranded on the $SCRATCH partition and you can manually move them back into the submission directory. Note the job ID from the file extension on out.### our err.###. You can move the files back to the submission directory with the command:
```$ mv $SCRATCH/######/* ./```
where you insert the jobid number in place of the '######' above. You can then proceed with extending the job as outlined above. The wildcard * in the above line of code instructs to move all files and directories in the directory $SCRATCH/######/ to the current directory that you are in.

<a name='deadlines'></a>

## Deadlines ##
1. Short update (few slides) on completed calculations: Tuesday 30 December during class (1 per group)
2. Final Presentation: Tuesday 7 December during class (1 per group)
3. Final Paper: Thursday 16 December by 5 PM (1 per group)

## Analysis ##

Do a detailed analysis of your system trying to identify trends in adsorption enery of OHx (x=0,1) species across similar surfaces of the materials in your group. You should at the least be able to generate a plot of average EO* vs average EOH* for the materials.

In addition to adsorption energy trends, you should be looking for trends with respect to changes in surface energy for adsorbing OHx species vs. for removing O and metal species from the surface in the form of defects. You will be tasked with plotting this data and then commenting on the observed behavior with respect to our hypothesis as stated at the top of this page.

### Requirements ###

At a minimum you should accomplish the following:

1. Complete the [HW5](../ASE/Getting_Started).
2. Complete all of the calculations outlined above for your material. This involves a set of adsorption and defect calculations on 2 distinct surface facets.
3. Analysis
4. Report (3~5 pages maximum)

### Presentation ##

The Final Presentation should include a summary of all of the calculations you were able to complete. Present the strucutres and adsorption energy trends. Also, take note of any major strucutral changes induced by the adsorbates.

Discuss any problems you had with running calculations. 

Each Group will also be required to ask questions and begin a disucssion about another groups work. Here are those assignments:

Ni group will ask Ti group questions.<br>
Ti group will ask Mg group questions.<br>
Mg group will ask Ni group questions.

<a name='report'></a>
### Final Report ###

The final report should be in the form of a 3-5 page long paper including figures and tables. Provide one report for each group. Please be succinct and organize it in the following way:

* Introduction (brief) - don’t write too much
* Calculation details
* Results and discussion including analysis.
* Conclusion (brief)

You are welcome to share data amongst your peers to discuss broader trends. 


