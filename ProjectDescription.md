---
layout: page
mathjax: true
permalink: /Project/
---
## Course Project ##


1. [Introduction](#intro)
2. [Edge Study Protocol](#protocol)
3. [Calculations](#calcs)


For the Final Project, you will perform a comprehensive study on CO<sub>2</sub> adsorption on mineral surfaces and its transformation due to interaction with these surfaces. The students will be assigned to two groups. One group will study magnesium-based oxides and silicates, while the other group will study calcium-based oxides and silicates. Each group will present their results in class that will be critiqued by the other groups. Finally, each group will jointly write a final report on the combined data.


Turn in your final report by emailing a PDF file to:

```
alevoj@seas.upenn.edu, yingjies@seas.upenn.edu
```
<a name='intro'></a>

## Introduction ##

Goal: Determine the lowest CO<sub>2</sub> adsorption energies and the molecular configurations.

Plan: Use DFT to calculate CO<sub>2</sub> adsorption energies for each adsorption site and initial configuration (details are given below). Determine the molecular bond lengths and bond angle for each converged calculation. Plot the Density of States (DOS) of the configuration with the lowest adsorption energy. 

### Motivation ###

- In designing an efficient carbon mineralization material, knowing if it can uptake CO<sub>2</sub> is critial. In other word, we are urged to know if CO<sub>2</sub> can adsorb (adsorption thermodynamics) on it and react? 
- Reactive magnesium oxide cement (RMC) is widely studied for carbon capture processes. Hydrated MgO has stronger abilities to adsorb CO<sub>2</sub> and produce carbonated products. We will determine the structure of hydrated MgO through analyzing the stability phase diagram that relates Gibbs free energy of dissociative particle adsorption to water derivative chemical potential (details will be given below).
- There are many possible bulk cleavage planes for MgO and CaO. The most studied ones are (100), (110), and (111). The stability of these facets varies greatly, leading to significantly different adsorption phenomena of CO<sub>2</sub> on different facets.
- Mg and Ca containing (alumino)silicates naturally exist in mine waste tailings. The use of these tailings is considered an alternative avenue for mineralization because they are an inherently cost-effective feedstock that could generate revenue in addition to captuing CO<sub>2</sub>. The presence of Si adds complexity to the study due to the different Si coordination numbers.
- When CO<sub>2</sub> is hydrated, it can undergo protonation/deprotonation through interaction with water to form charged intermediates depending on the solution pH. These species typically are HCO<sub>3</sub><sup>-</sup> and CO<sub>3</sub><sup>2-</sup>. This motivates us to study how these species interaction with mineral surfaces is differed from molecular CO<sub>2</sub>.
- Therefore, we will be studying the chemistry of carbon mineralization and down the line understanding how that chemistry changes with respect to the many factors, such as cation, facet, presence of water, Si coordination, pH conditon, etc.
<a name='protocol'></a>

## Surface Facet (Cleavage Plane) ##

Both MgO and CaO are simple compounds with a cubic crystal. Cleavage is the tendency of a crystalline material to break along definite planes, creating smooth surfaces. In general, we compare the cleavage energy to determine whether a certain facet is stable. 

$$
\epsilon_{cleave}(i,j) = -\frac{1}{2A}(E_{slab}^{i}+E_{slab}^{j}-2E_{bulk})
$$
