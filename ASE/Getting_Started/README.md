---
layout: page
mathjax: true
permalink: /ASE/Getting_Started/
---

# ASE Tutorials
1. [Introduction to ASE](../)
2. [Getting Started with DFT Calculations](../Getting_Started/)

____

## Getting Started with DFT Calculations ##

In the first exercise, we will be studying MXenes, how to determine their lattice constants, and O adsorption on the surface. For Homework 5, everyone will be studying the same system Ti<sub>2</sub>C. 

## Contents ##

1. [A Typical ASE Script](#a-typical-ase-script)
2. [Lattice Constant Determination](#lattice-constant-determination)
3. [Convergence with k-points](#convergence-with-k-points)
4. [Optimization](#optimization)


<a name='a-typical-ase-script'></a>

### A Typical ASE Script ###

ASE scripts can be run directly in the terminal (in the login node) or submitting to external nodes. Generally, you will be submitting jobs to external nodes and only small scripts will be run on the login node. By default, all output from any submitted script will be written *from the directory where the submission command was executed*, so make sure you are inside the calculation folder before running the submission command.

To start this tutorial and the exercises that follow, log on to Anvil and download the following:
```bash
wget https://upenncbe544.github.io/CBE544-2022.github.io/ASE/HW5_Compressed.tar.gz
tar -zxvf HW5_Compressed.tar.gz
cd HW5_Compressed
```

There are two files that are necessary to run jobs on the Anvil cluster. The first is `qe.sub`; this is the file that tells the scheduler how much time the job is allowed, how many processors it requires, and other pertinent information. First, notice the comments in the beginning. These include information such as how much time to allocate, the number of nodes required, what the names of the output and error files are, what the name of the job should be, and what your email is. 

```bash
                                                                                                                                                                     
#!/bin/bash
#SBATCH -J jobname #Job name
#SBATCH -N 1 #number of nodes used
#SBATCH --tasks-per-node=128 #number of tasks/node
#SBATCH -t 1:00:00 #Maximum job length
#SBATCH -o output #node Output file
#SBATCH -e error #node error file
#SBATCH --mail-user=your email #
#SBATCH --mail-type=ALL
#SBATCH -A EVE210010 #Allocation name, do not change

cd $SLURM_SUBMIT_DIR #Move to supply directory

module load openmpi #load openmpi, to prepare

python converging_scf.py

mpirun -np 120  /home/x-yamilee/q-e-qe-7.1/bin/pw.x -i scf.in > scf.out # 
```

The line ```python converging_scf.py``` picks the script you want to run. Therefore, you need to change the name of the file depending on which script you are running. We will be using this script later in this section for performing calculations to compute the lattice constant of bulk Ti<sub>2</sub>C.


Let's look at how a typical ASE script for geometry optimization is written. Open the `converging_scf.py` script, which will be used in a later section to perform create the input for a simple Ti<sub>2</sub>C relaxation. We import all the relevant ASE modules in for this calculation

```python
from ase import io
from espresso import espresso
from ase.optimize import BFGS

```

`from ase import io` imports the input/output commands for trajectory files, `from espresso import espresso` imports the Quantum Espresso calculator for the ASE interface.

An existing trajectory can be read in:

```python
mxene=io.read('init.traj') #read slab

```

Then, the Quantum ESPRESSO calculator is set up. All parameters related to the electronic structure calculation are included here. The following example shows typical parameters that we use in the group for calculations involving oxides.

```python
calc = espresso(pw=700,             #plane-wave cutoff
                dw=7000,                    #density cutoff
                xc='PBE',          #exchange-correlation functional
                kpts=(kx,ky,kz),       #k-point sampling;
                nbands=-20,             #20 extra bands besides the bands needed to hold valence electrons
                sigma=0.1,
                nosym=True,
                convergence= {'energy':1e-6,
                    'mixing':0.1,
                    'nmix':10,
                    'mix':4,
                    'maxsteps':500,
                    'diag':'david'
                    },  #convergence parameters for SCF
                 dipole={'status':True}, #dipole correction to account for periodicity in z
                 spinpol=False,
                 output = {'avoidio':False,
                    'removewf':True,
                    'wf_collect':False},
                 onlycreatepwinp = 'scf.in',
                 parflags='-npool 2 -nd 25',
                 outdir='calcdir')   #output directory for Quantum Espresso files
```




<a name='lattice-constant-determination'></a>

#### Lattice Constant Determination ####

Find the `lattice.py` script in the `lattice` folder. This script calculates the different energies of the system as a function of the lattice constant. Before you run this job, make sure you read the comments within to understand what it does.

Remember to add your email in the `qe.sub` file to receive notifications on the job! Submit the script by running:

```bash
sbatch qe.sub
``` 

The output trajectory `scf.out` contains information on the energy of the system with respect to the given lattice constant. To obtain the lattice constant that minimizes the energy, you will be writing a simple Python script to perform an Equation of State fit of the obtained energies as a function of the lattice constant. 

To proceed with writing this script, you will be modifying the example script provided here: [ASE-Equation of State](https://wiki.fysik.dtu.dk/ase/tutorials/eos/eos.html). Note that the sample script reads 5 configurations from the trajectory, but we have more configurations than that in our calculations. This script can be run on the login node directly. To execute the script you have written, use the command: 
```bash
python xyz.py
```
The output plot (`xyz.png`) should show the fitted energies as a function of the lattice volumes, with the volume corresponding to the minimum and the bulk modulus displayed on the top. Use this, and the fact that we have a cubic lattice to determine the DFT lattice constant.

**HW 5:** Show your Python script, Plot the Equation of State fit, and report the DFT lattice constant.

<a name='convergence-with-k-points'></a>

#### Convergence with k-Points ####
Next,you will be running the `kptconv.py` script in the `k-points` folder. Look through the script to understand what its doing. Run this script by submitting a job to an external node as discussed previously. Remember to change the name of the script to execute, in the `qe.sub` file. Upon completion, the script outputs a convergence plot and prints the total energies as a function of the k-points used in the calculation.

From the plot, and your understanding of concepts in DFT, suggest your pick for the k-points and the rationale behind your choice.

**HW 5:** Show the k-point convergence plot, your pick for the k-points, and your rationale.

<a name='optimization'></a>

#### Optimization ####
You will then be performing a geometry optimization on Ti<sub>2</sub>C. To proceed with this exercise, first take a look at the starting structure `init.traj` in the `relax` folder by using the GUI. You should see a 2x2x1 surface of Ti<sub>2</sub>C. You will be using this script for running the surface optimization calculations. Before submitting the job, please modify the following line (in addition to the script to run) in the `qe.sub` file:

```bash
#SBATCH --mail-user=miloue98@gmail.com
```
Take a look at the `scf.out` file after submitting the job. If you press `Esc`, the capital `G`, you should see:
```bash
=------------------------------------------------------------------------------=
   JOB DONE.
=------------------------------------------------------------------------------=
```
You can take a better look at the convergence criteria by doing `Esc` + `/converged`. You should see something that looks like this:
```bash
Total force =     0.000725     Total SCF correction =     0.000086
     SCF correction compared to forces is large: reduce conv_thr to get better values
     Energy error            =      8.8E-05 Ry
     Gradient error          =      5.3E-04 Ry/Bohr

     bfgs converged in  11 scf cycles and  10 bfgs steps
     (criteria: energy <  1.0E+00 Ry, force <  1.9E-03 Ry/Bohr)

     End of BFGS Geometry Optimization

     Final energy             =    -246.9201988489 Ry
```

This gives us the final energy in Rydbergs. 1 Ry = 13.605684 eV. If you want the energy in eV directly you can get it using ASE (python):
```python
from ase.io import read
final_traj = read('scf.out')
print(final_traj[-1].get_total_energy())
```
#### Adsorption ####
Finally, you will be calculating the adsorption energy of O on the Ti<sub>2</sub>C surface. Adsorption energy is given by:

$$
\Delta E_\mathrm{ads} = E_\mathrm{surface +ads}  - E_\mathrm{surface} - E_\mathrm{ads}
$$

To determine the most favorable site, you will adosrb O onto each of the 4 high symmetry sites and relax the structure. This will give the different total energies, that you can use to calculate adsorption energies. From there, you will be able to determine what is the most favorable site and whether or not Ti<sub>2</sub>C should be oxidized. 
Relax structures using instructions in the **Optimization** section.

**HW 5:** Report the converged energy of the optimized structure. 

**You must succesfully complete this task before proceeding to the Final Project**

