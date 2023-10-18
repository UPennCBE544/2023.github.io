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

In the first exercise, we will be studying alkaline oxides and how to determine their lattice constants, followed by surface relaxation of the (100) surface of the alkaline oxides. For Homework 5, everyone will be studying the same system MgO (100). 

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
wget https://upenncbe544.github.io/CBE544-2023.github.io/ASE/HW5_2023.tar.gz
tar -zxvf HW5_2023.tar.gz
cd HW5
```

There are two files that are necessary to run jobs on the Anvil cluster. The first is `anvil.sub`; this is the file that tells the scheduler how much time the job is allowed, how many processors it requires, and other pertinent information. First, notice the comments in the beginning. These include information such as how much time to allocate, the number of nodes required, what the names of the output and error files are, what the name of the job should be, and what your email is. 

```bash
#!/bin/bash

#SBATCH -J jobname   #Job name
#SBATCH -o out.%j #name of stdout output file
#SBATCH -e err.%j #name of stderr error file
#SBATCH -p wholenode  #queue type
#SBATCH -N 1 #no.of nodes
#SBATCH --ntasks-per-node 128 #no.of mpi tasks
#SBATCH -t 96:00:00 #maximum run time (hh:mm:ss)
#SBATCH --mail-user=abc@gmail.com #provide your email for notification
#SBATCH --mail-type=all #notify when job finishes
#SBATCH -A dmr180108  #Allocation (don't change this)


cd $SLURM_SUBMIT_DIR #Move to supply directory

mkdir -p $SCRATCH/$SLURM_JOBID
cp pw.in $SCRATCH/$SLURM_JOBID

cd $SCRATCH/$SLURM_JOBID

mpirun /home/x-syj1022/bin/qe-7.2/bin/pw.x -nd 4   <pw.in> pw.out

cp pw.out $SLURM_SUBMIT_DIR 
```

Let's look at how a typical ASE script for geometry optimization is written. The code below shows an example of a `relax.py` script, which will be used in a later section to perform create the input for a simple MgO relaxation. We import all the relevant ASE modules in for this calculation.

```python
from ase import io
from ase import Atoms
from espresso import espresso
from ase.optimize import BFGS
```

`from ase import io` imports the input/output commands for trajectory files, `from ase import Atoms` imports the Atoms object, useful editing and manipulating the system, `from espresso import espresso` imports the Quantum Espresso calculator for the ASE interface, and `from ase.optimize import BFGS` imports the BFGS algorithm to perform geometry optimization.

An existing trajectory can be read in:

```python
slab =  io.read('init.traj')   #read slab
slab.set_pbc([True,True,True])   #set periodic boundaries in all directions to true
```

Then, the Quantum ESPRESSO calculator is set up. All parameters related to the electronic structure calculation are included here. The following example shows typical parameters that we use in the group for calculations involving oxides.

```python
calc = espresso(pw=500,             #plane-wave cutoff
                dw=5000,                    #density cutoff
                xc='BEEF-vdw',          #exchange-correlation functional
                kpts=(4,4,1),	    #k-point sampling;
                nbands=-30,             #30 extra bands besides the bands needed to hold
                sigma=0.1,
                opt_algorithm = 'bfgs',
                fmax = 0.03,
                nstep=200,
                nosym=True,
                convergence= {'energy':1e-5,
                    'mixing':0.1,
                    'nmix':10,
                    'mix':4,
                    'maxsteps':2000,
                    'diag':'david'
                    },  #convergence parameters
                 dipole={'status':True}, #dipole correction to account for periodicity in z
                 output = {'avoidio':False,
                    'removewf':True,
                    'wf_collect':False},
                 spinpol=False,
                 parflags='-npool 1',
                 onlycreatepwinp='pw.in',
                 outdir='calcdir')   #output directory for Quantum Espresso files
```




<a name='lattice-constant-determination'></a>

#### Lattice Constant Determination ####

Find the `lattice.py` script in the `lattice` folder. This script calculates the different energies of the system as a function of the lattice constant. Before you run this job, make sure you read the comments within to understand what it does.

You should check the `traj` file (in this case it is `mgo-bulk.traj`) to make sure you have the correct starting structure that you work with. Also, make sure the `lattice.py` is reading the correct input file. Next, execute `lattice.py` by:

```bash
python lattice.py
```

This should generate a `pw.in` file and this is the input into the submission script `anvil.sub`. Remember to add your email in the `anvil.sub` file to receive notifications on the job! Submit the script by running:

```bash
sbatch anvil.sub
``` 

The output trajectory `pw.out` contains information on the energy of the system with respect to the given lattice constant. To obtain the lattice constant that minimizes the energy, you will be writing a simple Python script to perform an Equation of State fit of the obtained energies as a function of the lattice constant. 

To proceed with writing this script, you will be modifying the example script provided here: [ASE-Equation of State](https://wiki.fysik.dtu.dk/ase/tutorials/eos/eos.html). Note that the sample script reads 5 configurations from the trajectory, but we have more configurations than that in our calculations. This script can be run on the login node directly. To execute the script you have written, use the command: 
```bash
python xyz.py
```
The output plot (`xyz.png`) should show the fitted energies as a function of the lattice volumes, with the volume corresponding to the minimum and the bulk modulus displayed on the top. Use this, and the fact that we have a cubic lattice to determine the DFT lattice constant.

**HW 5:** Show your Python script, Plot the Equation of State fit, and report the DFT lattice constant.

<a name='convergence-with-k-points'></a>

#### Convergence with k-Points ####
Next,you will be running the `kptconv.py` script in the `k-points` folder. Look through the script to understand what its doing. Run this script by submitting a job to an external node as discussed previously. Remember now you run with `kptconv.py` not `lattice.py`. Upon completion, the script outputs a convergence plot and prints the total energies as a function of the k-points used in the calculation.

From the plot, and your understanding of concepts in DFT, suggest your pick for the k-points and the rationale behind your choice.

**HW 5:** Show the k-point convergence plot, your pick for the k-points, and your rationale.

<a name='optimization'></a>

#### Optimization ####
You will then be performing a geometry optimization on MgO. To proceed with this exercise, first take a look at the starting structure `mgo-100.traj` in the `relax` folder by using the GUI. You should see a 4x4x4 surface of MgO (100), with its bottom two layers fixed. You will be using this script for running the surface optimization calculations. Before submitting the job, please modify the following line (in addition to the script to run) in the `anvil.sub` file:

```bash
#SBATCH --mail-user=abc@gmail.com #provide your email for notification
```
Take a look at the `pw.out` file after submitting the job. If you press `Control` + `-`, followed by `Control` + `V`, you should see:
```bash
=------------------------------------------------------------------------------=
   JOB DONE.
=------------------------------------------------------------------------------=
```
You can take a better look at the convergence criteria by doing `Control` + `W` and search for `Final energy`. You should see something that looks like this:
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

This gives us the final energy in Rydbergs. 1 Ry = 13.605684 eV. If you want the energy in eV directly you can first run the command:

```python
python pwlog.py ./pw.out opt.traj
```
Next, in the python editor, enter the codes below:

```python
from ase.io import read
a = read('opt.traj')
print(a.get_total_energy())
```
#### Adsorption ####
In this part, you will be asked to plot the density of states (DOS) of the given structure. To do this, you simply just need to add a few lines in your `anvil.sub`:

```bash
/home/x-syj1022/apps/qe-7.2/bin/dos.x -in dos.in > dos.out  #DOS calculations
```
You also need to make sure that all the needed files are transferred to `$SCRATCH` and back once the job is done. Therefore, your `anvil.sub` now should look like this:

```bash
cd $SLURM_SUBMIT_DIR #Move to supply directory

mkdir -p $SCRATCH/$SLURM_JOBID
cp pw.in $SCRATCH/$SLURM_JOBID
cp dos.in $SCRATCH/$SLURM_JOBID

cd $SCRATCH/$SLURM_JOBID

mpirun /home/x-syj1022/apps/qe-7.2/bin/pw.x -nd 4   <pw.in> pw.out

/home/x-syj1022/apps/qe-7.2/bin/dos.x -in dos.in > dos.out  #DOS calculations

cp pw.out $SLURM_SUBMIT_DIR
cp dos.out $SLURM_SUBMIT_DIR
cp dos.dos $SLURM_SUBMIT_DIR
```

Note, in addition to `pw.in` as the input into `anvil.sub`, we now also need `dos.in`. This is already provided and should look like this:

```bash
&dos
   prefix='calc',
   outdir='.'
   Emin=-80.0
   Emax=20.0
   fildos='dos.dos'
/
```

Upon completion, the `dos.dos` file saves the data you need for the plot. You are then responsible for writing a script to generate the plot.

Finally, you will be calculating the adsorption energy of CO<sub>2</sub> on the MgO (100) surface. Adsorption energy calculation is given by:

$$
\Delta E_\mathrm{ads} = E_\mathrm{MgO+CO_{2}}  - E_\mathrm{MgO} - E_\mathrm{CO_{2}}
$$

You may use -1090.607 eV for E<sub>CO<sub>2</sub></sub>.

To adsorb an atom onto an oxygen, click on the oxygen you want to adsorb onto (for the example of MgO the surface is symmetric, therefore all the oxygens are equivalent). Then go to `Edit' -> 'Add atoms`. Alternatively, you can use `control`+`A`. Type in the symbol of element (e.g., C, O) and then select the relative coordinates. Finally, click on `Add` and the new atom should appear.

**HW 5:** Report the converged energy of the optimized structure. 

**You must succesfully complete this task before proceeding to the Final Project**

