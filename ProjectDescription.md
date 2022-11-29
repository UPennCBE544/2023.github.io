---
layout: page
mathjax: true
permalink: /Project/
---
## Course Project ##


1. [Introduction](#intro)
2. [Edge Study Protocol](#protocol)
3. [Calculations](#calcs)


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
- 
<a name='protocol'></a>

## Edge Study Protocol ##
   - Edges can be identified using their Miller Indexes (h, k , l)
- One rule of thumb in modeling a system is to start with a simple model and add complexities later
- Therefore we will be starting with the first sets of low Miller Indexes where            *-2 ≤ h, k, l ≤ 2*
    - To further simplify our cuts we will be setting l = 0. Because MXenes are 2D systems, to a first approximation, cuts in the z direction should not affect the edge cut. Which means 010 ≈ 011 ≈ 012 …
        
        **Note: MXenes are not completely 2D like graphene is so you may find small differences between the aforementioned edges. However since we will studying MXenes of at most 3 atoms, this is an appropriate approximation.**
        
    - This leaves us with possibility 25 edges where *-2 ≤ h, k ≤ 2 & l = 0*
    - The (0, 0, 0) edge does not exists - it’s just the bulk. We are now at 24 edges
    - After we account for cuts that are not unique such as (220), (330), (-1-10) that can all be reduced to (110).
    - We are left with (010), (110), (-210), (2-10), (1-20), (-120), (120), (210)
    - Because of the six fold symmetry of MXenes this list can be simplified to the unique MXenes (010), (110), (1-20), (120)
    - Every edge corresponds to an angle 0° *≤ ϴ ≤ 60°
    - (010) is the Zigzag edge with an angle *ϴ = 30°*
    - (110) and (1-20) are the Armchair edges with angles *ϴ = 0° and 60°* respectively
    - (120) edge has an angle *ϴ = 10.89°*
    - The maximum angle size comes from the MXene’s 6-fold symmetry
    - In an effort to keep the model as simple as possible, we will only model the high symmetry edges (010), (110), (1-20)
    
<a name='calcs'></a>
## Calculations ##


- Modeling
    - There are a number of terminations and termination combinations we can look at however we will only be considering 2 terminations: no termination & O termination. This choice is motivated through experiment
        - From Ti2C synthesis we know to a first approximation the basal planes can be terminated by O (Termination is synthesis dependent). Oxygen will therefore be one of our terminations
        - To get an idea of hour basal plane termination affects edge chemistry, it is import to look at how edges behave with respect to bare basal plane. Our second termination will be no termination
    - From there each student has 6 models
        - bare
            - (010)
            - (110)
            - (1-20)
        - O-term
            - (010)
            - (110)
            - (1-20)
            
- Running Calculations

   ```bash
   cd
   cd CBE544
   mkdir Final_Project
   cd Final_Project
   mkdir MXenes
   cd MXenes
   mkdir M2C
   cd M2C
   mkdir bare
   mkdir O-term
   cd bare

   ```
Please check [Assignment](ProjectAssignments.md) in order to determine which material you will be working with. 

- Edges are cuts made to the basal plane (relaxed). We will start by first relaxing the basal plane. Wherever you see ##M## replace it with your metal
    
    ```bash
    cp -R /home/x-yamilee/CBE544/Final_Project/Scripts ~/CBE544/Final_Project/
    cd /home/x-yourusername/CBE544/Final_Project/Scripts/
    vi qe.sub
    ```
- Change the email in *qe.sub*
    ```bash
    cd 
    cd CBE544/Final_Project/MXenes/M2C/bare
    cp /home/x-yamilee/CBE544/Final_Project/Final_Proj_Structures/M2C ./init.traj
    cp /home/x-yourusername/CBE544/Final_Project/Scripts/converging_scf.py ./
    cp /home/x-yourusername/CBE544/Final_Project/Scripts/qe.sub ./
    ```
    
-  submit the job

   ```bash
   sbatch qe.sub
   ```

- Repeat this step for *O-term*:

   ```bash
   cd ../O-term
   cp /home/x-yamilee/CBE544/Final_Project/Final_Proj_Structures/M2CO2 ./init.traj
   cp /home/x-yourusername/CBE544/Final_Project/Scripts/converging_scf.py ./
   cp /home/x-yourusername/CBE544/Final_Project/Scripts/qe.sub ./
   sbatch qe.sub
   
   ```

- Then following the exact steps you did for *bare*

### Directory Correction ###
   ```bash
   cd /home/x-yamilee/CBE544/Final_Project/MXenes/M2C
   mkdir Bulk
   mv bare Bulk
   mv O-term Bulk
   ```
- Edge Calculations
   ```bash
   cp /home/x-yamilee/CBE544/Final_Project/Scripts/edges_submission1.py /home/x-yourusername/CBE544/Final_Project/Scripts/edges_submission1.py
   vi edges_submission1.py
   ```
Inside edges_submission1.py, change the #M# & #username# variables to your assigned metal and username respectively. Then run edges_submission1.py.
   ```bash
   python edges_submission1.py
   ```
If the file runs with no error. Go back into edges_submission1.py and remove the '#' in front of line 40 and run edges_submission1.py again.
- Modeling Nanoparticles
Since the nanoparticles are large. We will be modeling static bare and O terminated nanoparticles. To do so, we first need to build the input structure for the nanoparticle.
```bash
   cd
   cd /home/x-yourusername/CBE544/Final_Project/MXenes/M2C/
   mkdir NP_bare_4
   mkdir NP_O-term_4
   cd /home/x-yourusername/CBE544/Final_Project/MXenes/M2C/Bulk/bare
   ase gui scf.out
   ```
 Follow the instructions in the powerpoint to cut the nanoparticle from the basal plane. Save the structure as init.traj in NP_bare_4/ 
 Copy the necessary files from the script and submit a static calculation for the nanoparticle. (Note: 4 at the end of NP_bare_4/ refers to the nanoparticle edge length)
 ```bash
   cd
   cp /home/x-yourusername/CBE544/Final_Project/Scripts/converging_scf.py /home/x-yourusername/CBE544/Final_Project/MXenes/M2C/NP_bare
   cp /home/x-yourusername/CBE544/Final_Project/Scripts/qe.sub /home/x-yourusername/CBE544/Final_Project/MXenes/M2C/NP_bare
   cd /home/x-yourusername/CBE544/Final_Project/MXenes/M2C/NP_bare
   ```
Before submitting the job, make sure mode in converging_scf.py is set to 'scf' and the k points are set to 1 1 1. Once you do, submit the job
   ```bash
   sbatch qe.sub
   ```
### Note ###
The instructions in the powerpoint are for a bare NP with 4 C long edges. To test for size effects, we will repeat all instructions for bare NP with 5 & 6 C long edges. To test for termination effects, we will repeat all the calculations usins the O terminated structures.

## Things to keep in mind ##
For all of the calculations that you are running here on out, you will need to submit them using qe.sub to the cluster. Do not run any of these jobs on the login node. 

Each calculation must be carried out in its own directory. All of the files necessary to carry out the calculation must also be in the directory at the time you submit the job or else it will not work. Usually this means the directory should have an ASE trajectory (init.traj), a relaxation script (converging_scf.py), and a submit script (qe.sub).

It is critical that you organize your directories consistently so that we can find the data later. 









