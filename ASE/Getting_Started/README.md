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

To start this tutorial and the exercises that follow, log on to Stampede2 and download the following:
```bash
wget https://upenncbe544.github.io/CBE544-2022.github.io/ASE/HW5.tar.gz
tar -zxvf HW5.tar.gz
cd HW5
```

There are two files that are necessary to run jobs on the Anvil cluster. The first is `stampede.sub`; this is the file that tells the scheduler how much time the job is allowed, how many processors it requires, and other pertinent information. First, notice the comments in the beginning. These include information such as how much time to allocate, the number of nodes required, what the names of the output and error files are, what the name of the job should be, and what your email is. 

```bash
                                                                                                                                                                     
#!/bin/bash
#SBATCH -J vc-relax
#SBATCH -N 1
#SBATCH --tasks-per-node=128
#SBATCH -t 1:00:00
#SBATCH -o output
#SBATCH -e error
#SBATCH --mail-user=miloue98@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -A EVE210010

cd $SLURM_SUBMIT_DIR

module load openmpi

python converging_scf.py

mpirun -np 120  /home/x-yamilee/q-e-qe-7.1/bin/pw.x -i scf.in > scf.out
```

Stampede requires us to run jobs on the $SCRATCH partition to reduce heavy I/O on the $WORK partition. Therefore the second block of code sets up environment variables for the submission directory and the $SCRATCH directory.

Then, the lines ```cp ${submit_dir}/sto-bulk.traj ${SCRATCH_DIRECTORY}``` and ```cp ${submit_dir}/lattice.py ${SCRATCH_DIRECTORY}``` copy the relevant files to the $SCRATCH partition. You will need to edit these lines depending on the name of the python script and input files required.

The line ```python lattice.py``` picks the script you want to run. Therefore, you need to change the name of the file depending on which script you are running. We will be using this script later in this section for performing calculations to compute the lattice constant of bulk SrTiO<sub>3</sub> perovskite.


Let's look at how a typical ASE script for geometry optimization is written. Open the `relax.py` script, which will be used in a later section to perform a simple optimization on a (001) SrTiO<sub>3</sub> slab . We import all the relevant ASE modules in for this calculation

```python
from ase import io
from ase import Atoms
from espresso import espresso
from ase.optimize import QuasiNewton
```

`from ase import io` imports the input/output commands for trajectory files, `from ase import Atoms` imports the Atoms object, useful for editing and manipulating the system, `from espresso import espresso` imports the Quantum Espresso calculator for the ASE interface, and `from ase.optimize import QuasiNewton` imports the Quasi Newton algorithm to perform geometry optimization

An existing trajectory can be read in:

```python
slab =  io.read('001-bo2-sto.traj') #read slab
slab.set_pbc([True,True,True])     #set periodic boundaries in all directions to True
```

Then, the Quantum ESPRESSO calculator is set up. All parameters related to the electronic structure calculation are included here. The following example shows typical parameters that we use in the group for calculations involving oxides.

```python
calc = espresso(pw=500,             #plane-wave cutoff
                dw=5000,                    #density cutoff
                xc='PBE',          #exchange-correlation functional
                kpts=(5,5,1),       #k-point sampling;
                nbands=-20,             #20 extra bands besides the bands needed to hold valence electrons
                sigma=0.1,
                nosym=True,
                convergence= {'energy':1e-5,
                    'mixing':0.1,
                    'nmix':10,
                    'mix':4,
                    'maxsteps':500,
                    'diag':'david'
                    },  #convergence parameters for SCF
                 dipole={'status':True}, #dipole correction to account for periodicity in z
                 output = {'avoidio':False,
                    'removewf':True,
                    'wf_collect':False},
                 spinpol=False,
                 parflags='-npool 2',
                 outdir='calcdir')   #output directory for Quantum Espresso files
```

Finally, the Quantum ESPRESSO calculator is attached to the `slab` Atoms object, the relaxation calculation is run, and the total energy of the system is output in the log file. 

To submit any job on Stampede2, use:

```bash
sbatch stampede.sub
```

<a name='lattice-constant-determination'></a>

#### Lattice Constant Determination ####

Find the `lattice.py` script in the `lattice` folder. This script calculates the different energies of the system as a function of the lattice constant. Before you run this job, make sure you read the comments within to understand what it does.

Remember to add your email in the `stampede.sub` file to receive notifications on the job! Submit the script by running:

```bash
sbatch stampede.sub
``` 

The output trajectory `sto.traj` contains information on the energy of the system with respect to the given lattice constant. To obtain the lattice constant that minimizes the energy, you will be writing a simple Python script to perform an Equation of State fit of the obtained energies as a function of the lattice constant. 

To proceed with writing this script, you will be modifying the example script provided here: [ASE-Equation of State](https://wiki.fysik.dtu.dk/ase/tutorials/eos/eos.html). Note that the sample script reads 5 configurations from the trajectory, but we have more configurations than that in our calculations. This script can be run on the login node directly. To execute the script you have written, use the command: 
```bash
python xyz.py
```
The output plot (`xyz.png`) should show the fitted energies as a function of the lattice volumes, with the volume corresponding to the minimum and the bulk modulus displayed on the top. Use this, and the fact that we have a cubic lattice to determine the DFT lattice constant.

**HW 5:** Show your Python script, Plot the Equation of State fit, and report the DFT lattice constant.

<a name='convergence-with-k-points'></a>

#### Convergence with k-Points ####
Next, we will determine how well-converged the total energy is with respect to the number of k-points in each direction. First, take a look at the `resize.py` script which resizes the lattice to the DFT lattice constant computed in the previous exercise. Run this script directly in the login node to obtain the starting structure for the k-point calculations. Next,you will be running the `kptconv.py` script in the `k-points` folder. Look through the script to understand what its doing. Run this script by submitting a job to an external node as discussed previously. Remember to change the name of the script to execute, in the `stampede.sub` file. Upon completion, the script outputs a convergence plot and prints the total energies as a function of the k-points used in the calculation.

From the plot, and your understanding of concepts in DFT, suggest your pick for the k-points and the rationale behind your choice.

**HW 5:** Show the k-point convergence plot, your pick for the k-points, and your rationale.

<a name='optimization'></a>

#### Optimization ####
Finally, you will be performing a geometry optimization on the (001) BO2-terminated surface of SrTiO<sub>3</sub>. To proceed with this exercise, first take a look at the starting structure `001-bo2-sto.traj` in the `relax` folder by using the GUI. You should see a 2x2x4 surface of SrTiO<sub>3</sub>, with the bottom two layers fixed to the bulk positions. Next, take a look at the `relax.py` script discussed previously. You will be using this script for running the surface optimization calculations. Before submitting the job, please modify the following lines (in addition to the script to run) in the `stampede.sub` file:

```bash
#SBATCH -p normal #queue type
#SBATCH -N 1 #no.of nodes
#SBATCH -t 48:00:00 #run time (hh:mm:ss)
```
Take a look at the `opt.log` file after submitting the job, you should see something like this:
```bash
BFGSLineSearch:   0[  0]  12:27:27   -60932.887293       2.0460
BFGSLineSearch:   1[  2]  14:16:48   -60933.316386       0.8612
BFGSLineSearch:   2[  3]  14:59:04   -60933.508996       0.5930
BFGSLineSearch:   3[  4]  15:32:44   -60933.583239       0.5081
BFGSLineSearch:   4[  5]  16:07:11   -60933.612633       0.2711
BFGSLineSearch:   5[  7]  17:02:23   -60933.626425       0.1151
```
The optimization step is printed in the first column, the wall clock time in the second, total energy (in eV) in the third, and the forces (in eV/Å) in the fourth and final column. Report the converged energy once the job finishes. Note that the forces should have converged to a value less than the cut-off specified in the `relax.py` script. 

**HW 5:** Report the converged energy of the optimized structure. 

**You must succesfully complete this task before proceeding to the Final Project**

