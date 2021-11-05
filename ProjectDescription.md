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

You should already be in the directory corresponding to your material. In this directory make 2 new directories: 110 and 111. Move into the newly created 110 directory. Make a directory called clean and move into it. Copy the relaxed bulk trajectory (opt.traj) from the bulk/relax directory corresponding to your material into your 110/clean directory. Copy the surf_build.py from the scripts directory into this directory as well. Now you will build a 4 layer (110) slab model using the surf_build.py script. Look at the surf_build.py script. It should look like the text below: 

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

You should already be in the directory corresponding to your material. In this directory make 2 new directories: 110 and 100. Move into the newly created 110 directory. Make a directory called clean and move into it. Copy the relaxed bulk trajectory (opt.traj) from the bulk/relax directory corresponding to your material into your 110/clean directory. Copy the surf_build.py from the scripts directory into this directory as well. Now you will build a 4 layer (110) slab model using the surf_build.py script. Look at the surf_build.py script. It should look like the text below: 

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

In the 100 directory repeat the same process but build the (100) surface and make it look exactly like the one below. You should now have an init.traj for both the (110) and the (100) surfaces in their own separate directories.

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

Now we want to set up and run relaxations for O and OH adsorptions at full coverage. We know from previous work that the favored adsorption site on these surfaces is above the B site cation. So let's start by setting up an adsorption calculation for O. In order to add an adsorbate, click on the atom that you want to add the adsorbate directly above. Press Ctrl+A or click Edit -> Add Atom and in the window that comes up type O in the top box. In the box below Position you can enter a number corresponding to the number of Angstroms above the highlighted Atom you want to place the adsorbate. For the first oxygen adsorbate try 2.3 A like in the image below. Look at the side view and make sure the adsorbed oxygen is in a similar position to the ones in my system. Repeat the same for all four adsorption sites until you have a system with four O adsorbates located above the B site atoms in the lattice.

<center><img src="../Images/perov_o_ads_full.png" alt="window" style="width: 800px;"/><br>
Process for adsorbing O at 1ML coverage on the (110) surface
</center>

Now in this directory copy the relax.py script and the stampede.sub script here from the scripts folder. Submit the job using sbatch stampede.sub. 

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

Step 4: Defects 1ML

Go to the directory clean/. Make a directory clean/vac/. Cd into vac/. Make two directories clean/vac/1ML/ and clean/vac/0.25ML/. Cd into 1ML. Make four directories: Sr, SrO, Ag, and AgO2. Copy the clean relaxed init.traj into each of these directories. For the Sr case, you will remove the 4 topmost Sr atoms, save the init.traj, and submit the job. For the SrO you will be removing the same 4 atoms plus the O atoms coordinated to them. Your final structures will resemble those shown below. Relax all of these.

<center><img src="../Images/rutile_surfs.png" alt="window" style="width: 800px;"/><br>
Schematic of Rutile Oxides Surfaces
</center>

Step 5: Defects 0.25 ML

Move into the clean/vac/0.25ML/ directory that you made in Step 4 above. Make four directories exactly the same as in Step 4: Sr, SrO, Ag, and AgO2. Repeat the same steps except instead of removing four atoms/units, you will remove only one. See images below for examples. Relax these.

<center><img src="../Images/rutile_surfs.png" alt="window" style="width: 800px;"/><br>
Schematic of Rutile Oxides Surfaces
</center>

On the (111) surface, you will repeat all of the same calculations as the (110) done previously (including repeating the slab model to create a 2x2 surface!). See images below for examples of each setup. 

### Rutile Oxides ###

On the (110) surface, adsorb H on the O and M sites at coverages of 1 adsorbate per surface and 1 adsorbate per site. Relax these systems. See the figures below for a guide to the adsorption sites.

On the (111) surface, adsorb H on the O and M sites at coverages of 1 adsorbate per surface and 1 adsorbate per site. Relax these systems. See the figures below for a guide to the adsorption sites.

### Task 3 ###

The final task will be creating surface defects and relaxing the systems. From these results, we will be able to calculate the defect formation energies. Please skip ahead to the section heading relevant to your assigned material.

Perovskites

On the (001)-AO terminated surface, create a surface layer A-site defect, O-site defect, and AO-defect. Also create a full surface O site defect (remove all surface layer O atoms), full surface A site defect. 

On the (001)-BO<sub>2</sub>-terminated surface, create a surface layer B-site defect, O-site defect, and BO<sub>2</sub>-defect. Also create a full surface O site defect (remove all surface layer O atoms) and full surface B site defect.

On the (111) surface, create a 

Rutile Oxides

On the (110) surface, create a metal defect, an O defect, and high coverage limit defects

### Organization ###

Organization for the project is very important so that the data is accesbile once the class is over. We will structure is to be something like this for the Perovskite 110 surface.

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


