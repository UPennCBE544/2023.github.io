---
layout: page
mathjax: true
permalink: /Project/
---
## Course Project ##


1. [Introduction](#intro)
2. [Calculations](#calcs)
3. [Deadlines](#deadlines)
4. [Analysis](#analysis)
5. [Final Report](#report)


For the Final Project, you will be studying the nature of MXene edges. The students will work in 2 groups. Each group will probe a set of MXenes whose M belongs to a specific row in the periodic table. Each group will present their results in class that will be critiqued by the other groups. Finally, each group  will jointly write a final report on the combined data. T


Turn in your final report by emailing a PDF file to:

```
alevoj@seas.upenn.edu, yamilee@seas.upenn.edu
```
<a name='intro'></a>

## Introduction ##

Goal: Determine electronic properties of MXene edges.

Plan: Use DFT to calculate edges energies, oxygen adsorption energies on the edges(oxidation) and size effect of nanoparticles.

### Motivation ###

- In designing a new catalyst , you want to know as much as possible about a material
- MXenes make for good catalyst candidates because they are highly tunable
- Imagine we have a machine, a MXene modifier with 7 knobs. You can turn the knobs on any of these cranks and modify the chemistry
- A lot of work has been done on MXene basal planes
- However not much is know about the edges
- This may seem trivial since it is just one knob in our MXene diagram. However we know that in most material chemistry changes with sites. For example in MoS2 (2D material - catalyst for hydrodesulfurization) the active sites are on the edges not the basal planes. Meaning, the chemistry in this case varies greatly from basal plane to edge.
- Therefore we will be studying MXene edges to understand their chemistry and activity & down the line understand how that chemistry changes with respect to the basal plane

Edge Study Protocol - Picking which edges to study
    - Edges can be identified using their Miller Indexes (h, k , l)
    - One rule of thumb in modeling a system is to start with a simple model and add complexities later
    - Therefore we will be starting with the first sets of low Miller Indexes where          *-2 ≤ h, k, l ≤ 2*
        - To further simplify our cuts we will be setting l = 0. Because MXenes are 2D systems, to a first approximation, cuts in the z direction should not affect the edge cut. Which means 010 ≈ 011 ≈ 012 …
            
            **Note: MXenes are not completely 2D like graphene is so you may find small differences between the aforementioned edges. However since we will studying MXenes of at most 3 atoms, this is an appropriate approximation.**
            
        - This leaves us with possibility 25 edges where *-2 ≤ h, k ≤ 2 & l = 0*
        - The (0, 0, 0) edge does not exists - it’s just the bulk. We are now at 24 edges
        - After we account for cuts that are not unique such as (220), (330), (-1-10) that can all be reduced to (110).
        - We are left with (010), (110), (-210), (2-10), (1-20), (-120), (120), (210)
        - Because of the six fold symmetry of MXenes this list can be simplified to the unique MXenes (010), (110), (1-20), (120)
        - Every edge corresponds to an angle 0° *≤ ϴ ≤ 60° given by:*
        
                         *$ϴ =cos^{-1}√{ {(3(h+k)^2}/4((h+k)^2-hk))}$*
        
        - (010) is the Zigzag edge with an angle *ϴ = 30°*
        - (110) and (1-20) are the Armchair edges with angles *ϴ = 0° and 60°* respectively
        - (120) edge has an angle *ϴ = 10.89°*
        - The maximum angle size comes from the MXene’s 6-fold symmetry
        - In an effort to keep the model as simple as possible, we will only model the high symmetry edges (010), (110), (1-20)
<a name='calcs'></a>
## Calculations ##

To get started we need to get the files needed for the final project. Move into your CBE544 directory, download the project files, and open the tar file by running these commands:

```bash
cd
cd CBE544
wget https://upenncbe544.github.io/CBE544-2021.github.io/FinalProject.tar.gz
tar -zxvf FinalProject.tar.gz

```

In the FinalProject directory you should see a directory called scripts, which contains the scripts you will need for this project. You will also see two directories: rutiles and perovskites. Please check [Assignment](ProjectAssignments.md) in order to determine which material you will be working with. The bulks for all of these materials have already been optimized using the calculator settings that we determined from HW 5 for the perovskite SrTiO<sub>3</sub>. Move into the directory corresponding to your material.

### Things to keep in mind ###
For all of the calculations that you are running here on out, you will need to submit them using stampede.sub to the cluster. Do not run any of these jobs on the login node. 

Each calculation must be carried out in its own directory. All of the files necessary to carry out the calculation must also be in the directory at the time you submit the job or else it will not work. Usually this means the directory should have an ASE trajectory (init.traj), a relaxation script (relax.py), and a submit script (stampede.sub).

It is critical that you organize your directories consistently so that we can find the data later. See the [Organization](#organization) section above for additional guidance. 

### Task 1 ###

You will build surfaces from the optimized bulks for your assigned materials and run the relaxation. In order to build the surfaces, you will be using the surf_build.py script in the scripts directory. This script utilizes the build module in ASE to cut a surface from the optimized bulk trajectory. It reorients the surface so that the surface is perpendicular to the z-direction in the output file. Please skip ahead to the section heading relevant to your assigned material.

### Perovskites ###

You should already be in the directory corresponding to your material. In this directory make 2 new directories: 110 and 111. Go ahead and build the full directory tree for your material that matches the relevant one in the [Organization](#organization) section above. Please make it exactly like this tree. !!!Note that the example tree below is for SrAgO3 and you will need to rename things to match your assigned material depending on the transition metal in your system (throughout this tutorial write-up you will need to keep this in mind whenever I mention Ag!!! Move into the 110 directory. Move into the clean directory. Copy the relaxed bulk trajectory (opt.traj) from the bulk/relax directory corresponding to your material into your 110/clean directory. Copy the surf_build.py from the scripts directory into this directory as well. Now you will build a 4 layer (110) slab model using the surf_build.py script. Look at the surf_build.py script. It should look like the text below: 

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

<center><img src="../Images/perovskites_surfs.png" alt="window" style="width: 600px;"/><br>
Schematic of Perovskite Surfaces. Both surfaces should have a total of 20 atoms if you built them correctly (4 Sr atoms, 4 B site atoms, 12 oxygen atoms).
</center>


### Rutile Oxides ###

You should already be in the directory corresponding to your material. In this directory make 2 new directories: 110 and 100. Go ahead and build the full directory tree for your material that matches the relevant one in the [Organization](#organization) section above. Please make it exactly like this tree. !!!Note that the example tree below is for MoO2 and you will need to rename things to match your assigned material depending on the transition metal in your system (throughout this tutorial write-up you will need to keep this in mind whenever I mention Mo!!! Move into the newly created 110 directory. Make a directory called clean and move into it. Copy the relaxed bulk trajectory (opt.traj) from the bulk/relax directory corresponding to your material into your 110/clean directory. Copy the surf_build.py from the scripts directory into this directory as well. Now you will build a 4 layer (110) slab model using the surf_build.py script. Look at the surf_build.py script. It should look like the text below: 

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

<center><img src="../Images/rutile_surfs.png" alt="window" style="width: 600px;"/><br>
Schematic of Rutile Oxides Surfaces. The (110) surface should have a total of 24 atoms if you built it correctly (8 metal atoms and 16 oxygen atoms). The (100) surface should have a total of 15 atoms if you built it correctly (5 metal atoms and 10 oxygen atoms).
</center>


### Everyone ###
Once you have built both surface facets for you material, constrain the bottom two layers of the cell to the bulk lattice positions. To do this select the atoms that you want to constrain, click Tools -> Constraints -> Constrain Selected Atoms. The constrained atoms should now have dashed 'X's on them. These should match the images above. Make sure that you are constraining a stoichiometric number of atoms. For the Rutile Oxides this means you should be constraining an integer multiple of MO2 atoms. For the Perovskites you should be constraining an integer multiple of ABO<sub>3</sub> atoms.

Next, we need to make sure we have the appropriate vacuum set up between slabs. We will use 20 Angtroms of vacuum between slabs. To do this, copy the script titled vacuum.py into the directory with your trajectory and look at the script. It reads in a file called 'init.traj', centers it in the unit cell, and adds 10 Angstroms of vacuum on both sides of the slab (axis=2 refers to the z-axis) and then rewrites the file as 'init.traj'.

Finally, we can relax these surfaces to get the initial structures that we will use for all of the adsorption and defect calculations going forward. Copy the relax.py script into each surface folder. Copy the stampede.sub script into each surface folder. Then, from each directory make sure your stampede.sub script is copying and running the relevant files. Then run the relaxation using the command ```sbatch stampede.sub```. 

### Task 2 ###

IMPORTANT: You cannot begin Task 2 until the relaxations of the surfaces in Task 1 are complete. This is because you will be adsorbing species on the surfaces and comparing the change in energy as a result of the adsorbate-surface interaction. You must start from the relaxed clean surface in order to capture this consistently.

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

<center><img src="../Images/rutile_surf_descrpt_110.png" alt="window" style="width: 400px;"/><br>
Rutile oxide (110) surface with adsorption site labels.
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
Highlighted atoms to delete in order to create defect surfaces on Rutile (110) surfaces
</center>

Step 5: Repeat everything above for the (100) surface

On the (100) surface, you will run adsorptions and defects similar to the (110) done previously (including repeating the relaxed slab model to create a 2x2 surface!). See images below for examples of each setup. Create and relax trajectories matching those below in the directories corresponding from the organization tree.

<center><img src="../Images/rutile_100_oads.png" alt="window" style="width: 400px;"/><br>
Process for adsorbing O at 1ML coverage on the (100) surface
</center>

<center><img src="../Images/rutile_100_ohads.png" alt="window" style="width: 800px;"/><br>
Process for adsorbing OH at 1ML coverage on the (100) surface
</center>

<center><img src="../Images/rutile_100_vacs.png" alt="window" style="width: 800px;"/><br>
Schematic of Perovskite (100) surface with 1ML defect coverages
</center>

#### Jobs not reaching force convergence ####

Some jobs may run for the allotted time on Stampede2 and not reach the force convergence. To check whether your job has reached convergence you can go into the directory containing the output files and run the command ```$ tail opt.log```. This will display the last few lines of the log file. The last column displays the force calculated at that step. If the final line shows a force value less than 0.03, then it means the calculation has converged, since we set that as the threshhold in relax.py. 

If the job has ended (due to failure, timeout, or some other reason) and this force value is not less than 0.03, it means the job needs to be rerun starting from the final structure. In order to do this, you will make a new directory (from inside the current directory) called 'extend'. 
```$ mkdir extend```
Then you can copy the files into this directory:
```$ cp opt.traj relax.py stampede.sub extend/```
Then move into the directory and rename opt.traj as init.traj:
```$ mv opt.traj init.traj```
Now you have everything that you need to start a new calculation from the output of the previous one, and you can submit with:
```$ sbatch stampede.sub```

### Jobs not being returned from $SCRATCH ###

Abiding by the rules for running on Stampede2, we have to run the IO intensive jobs from the $SCRATCH partition. This is covered on this website where I introduced you to the submission script. Basically when the job begins, eveything that is needed to run it is copied to the $SCRATCH directory with the name of the jobid, then the job runs, then the files are moved back to the submission directory and the $SCRATCH directory is deleted. Sometimes, the job fails prematurely and the stampede.sub script is unable to complete resulting in files getting stranded in the $SCRATCH directory. If this happens, you will receive an email that the job failed, but when you go to the submission directory you will not see any output files (e.g. opt.traj, calcdir, etc...). The files are simply stranded on the $SCRATCH partition and you can manually move them back into the submission directory. Note the job ID from the file extension on out.### our err.###. You can move the files back to the submission directory with the command:
```$ mv $SCRATCH/######/* ./```
where you insert the jobid number in place of the '######' above. You can then proceed with extending the job as outlined above. The wildcard * in the above line of code instructs to move all files and directories in the directory $SCRATCH/######/ to the current directory that you are in.

<a name='deadlines'></a>

## Deadlines ##
1. Submit all calculations from Tasks 1 and 2 by Wednesday November 24
2. Short update (few slides) on completed calculations: Tuesday 30 November during class (1 per individual)
3. Final Presentation: Tuesday 7 December during class (1 per group)
4. Final Paper: Tuesday 14 December by 5 PM (1 per group)

<a name='analysis'></a>

## Analysis ##

Do a detailed analysis of your system trying to identify trends in adsorption energy of OHx (x=0,1) species across similar surfaces of the materials in your group. You should at the least be able to generate a plot of average EO* vs average EOH* for the materials.

In addition to adsorption energy trends, you should be looking for trends with respect to changes in surface energy for adsorbing OHx species vs. for removing O and metal species from the surface in the form of defects. You will be tasked with plotting this data and then commenting on the observed behavior with respect to our hypothesis as stated at the top of this page.

A set of electronic structure calculations will be required once the full set of relaxation calculations has been completed.

### Requirements ###

At a minimum you should accomplish the following:

1. Complete the [HW5](../ASE/Getting_Started).
2. Complete all of the calculations outlined above for your material. This involves a set of adsorption and defect calculations on 2 distinct surface facets.
3. Analysis
4. Report (3~5 pages maximum)

### Presentation ##

The Final Presentation should include a summary of all of the calculations you were able to complete. Present the strucutres and adsorption energy trends. Also, take note of any major structural changes induced by the adsorbates.

Discuss any problems you had with running calculations.

Each Group will also be required to ask questions and begin a discussion about another groups work. Here are those assignments:

Group 1 will ask Group 4 questions.<br>
Group 2 will ask Group 3 questions.<br>
Group 4 will ask Group 2 questions.<br>
Group 3 will ask Group 1 questions.

<a name='report'></a>
### Final Report ###

The final report should be in the form of a 3-5 page long paper including figures and tables. Provide one report for each group. Please be succinct and organize it in the following way:

* Introduction (brief) - don’t write too much
* Calculation details
* Results and discussion including analysis.
* Conclusion (brief)

You are encouraged to share data amongst your peers to discuss broader trends. 


