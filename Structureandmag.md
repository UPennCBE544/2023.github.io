---
layout: page
mathjax: false 
permalink: /MagStruct/
---

# Magnetism and Structure 

## Magnetism 

## How to Edit Magnetism ##

To edit magnetism I would first make a new directory in the directly where you have done the original calculation (mkdir magmom). Copy into this directory the mag.traj, vasp-ase.sub, and relax.py or opt-ads.py script (whichever you used for the calculation). The main idea for this to change the magnetism to be in an absolute minimum value and from some prior calculation we believe that pattern to be Co having a spin of 2.8 on top and the bottom, and nearly 0 everywhere else. So what we will do is go into the directory for the new calculation and try to set the trajectory in this magnetism by opening the mag.traj editing the magnetism (select the atom, edit, modify and then put in the new magnetic moment and hit save) then saving the new trajectory as init.traj (or whatever you call you initial structure in your relaxation script). You can see this pattern on below with any doubts. 

### Clean LiCoO<sub>2</sub> Magnetism ###

<center><img src="../Images/LiCoO2mag1.png" alt="window" style="width: 400px;"/><br>
</center>
<center><img src="../Images/LiCoO2mag2.png" alt="window" style="width: 400px;"/><br>
104 Desired Magnetism
</center>

When you add an adsorbate there tends to be a large perturbation to the magnetism. This however does not mean that the absolute minimum chagned and you can see that our desired magnetism should still be about the same. 

<center><img src="../Images/LiCoO2adsmag2.png" alt="window" style="width: 400px;"/><br>
</center>
<center><img src="../Images/LiCoO2adsmag2.png" alt="window" style="width: 400px;"/><br>
104 Desired Magnetism with adsorbate
</center>

### Doped LiCoO<sub>2</sub> Magnetism ###

With an adsorbate we may see slight variations. Here is an example for Al doped LiCoO<sub>2</sub>. Here is an example of the magnetism with Al. You can see that there is a slight change when it is in row 1 and row2. These exact moments should be similar for Ti and Mg but Ni may have a different magnetic moment than 0. That being said these magnetic patterns are likely the minimums for the systems that you are running and should be what we aim for when rerunning the caclulations. 

<center><img src="../Images/AlLiCoO2mag1.png" alt="window" style="width: 400px;"/><br>
<img src="../Images/AlLiCoO2mag2.png" alt="window" style="width: 400px;"/><br>
<img src="../Images/AlLiCoO2mag3.png" alt="window" style="width: 400px;"/><br>
<img src="../Images/AlLiCoO2mag4.png" alt="window" style="width: 400px;"/><br>
104 Al Doped Desired Magnetism
</center>


## Strucutre

### Examples of Structural Changes to look for ##

#### LiCoO<sub>2</sub> without Al #### 
<center>
<img src="../Images/structure1.png" alt="window" style="width: 400px;"/><br>
<img src="../Images/structure2.png" alt="window" style="width: 400px;"/><br>
</center>
#### LiCoO<sub>2</sub> with Al #### 
<center>
<img src="../Images/structure3.png" alt="window" style="width: 400px;"/><br>
<img src="../Images/structure4.png" alt="window" style="width: 400px;"/><br>
<img src="../Images/structure5.png" alt="window" style="width: 400px;"/><br>
</center>
