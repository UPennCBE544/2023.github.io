---
layout: page
mathjax: false 
permalink: /Clusters/
---

# Getting Started
1. [Logging Into the Computing Clusters](../Clusters/)
2. [Basic UNIX](../UNIX/)
3. [Python](../Python/)

____

## Logging Into the Computing Clusters

Once you account on Chestnut has been activated. Follow the instructions and tests to make sure everything is set up properly and functional.

## Contents
1. [Installation](#installation)
2. [Logging On](#logging)
3. [First Time Logging On](#first-time)
4. [Making Sure Everything Works](#testing)

<a name='installation'></a>

## Installation

### Mac OSX
Download and install:

* [XQuartz](http://www.xquartz.org/)

To prevent X11 from timing out, open the terminal and type:

```bash
mkdir -p ~/.ssh
echo $'\nHost *\n ForwardX11Timeout 1000000\n' >>~/.ssh/config
```


### Windows

Download and install:

* [PuTTY](http://www.putty.org/)
* [Xming](http://sourceforge.net/projects/xming/) (Note: disable automatic installation of PuTTY with Xming. The above installer is a newer version)


### Linux (Debian-based, e.g. Ubuntu)
From the terminal
____

<a name='logging'></a>

## Logging onto the Clusters

Your chestnut username will be your username of logging into Penn systems.

Follow the instructions below for your system:

### Mac OSX

Open "Applications-> Utilities -> Terminal" or "Command+Space" to search Terminal using "spotlight search"
In a terminal:
```bash
ssh -X username@chestnut-login.seas.upenn.edu
```

### Windows 
Launch Xming. You will always need to have this open in order to forward graphical windows from the external clusters.

Start PuTTY, and:

* “Session” → “Host Name” `username@chestnut-login.seas.upenn.edu` for Chestnut
* “Connection” → “SSH” → “X11” check “Enable X11 forwarding”
* Back in “Session”, you can **save these settings for next time**.

You can start putty several times, if you need several terminal windows; only one instance of Xming needed.


### Linux ###

In a terminal:
```bash
ssh -X username@chestnut-login.seas.upenn.edu
```
____

<a name='first-time'></a>

### First Time Logging in ###

For the **first login** only, add the following line to your ~/.bashrc script command:

``` bash
source /scratch/alevoj1/Scripts/group_bashrc
````
In order to do this you must use a text editor in the terminal. VIM and Nano are two popular text editors. Introductions to both of these are on the fron page of the group website. So do either:

```bash
vi ~/.bashrc
```
OR
```bash
nano ~/.bashrc
```

and paste the source line there. Once that is done exit the text editor and run this command in the terminal.

```bash
source ~/.bashrc
```

**Change the permission of files:**

```bash
cd
cd work
mkdir CBE544
chmod g+r ./
```

Do all of your work in the directory you have called CBE544.

For example, Create a folder `hw5` under `CBE544` and perform all calculations of HW5 in `hw5`. If you have already started you jobs somewhere else, you can copy the entire folder to CBE544 once you have done all calculations (e.g. `cp -r folderpath ~/work/CBE544/hw5`).  However, you need to complete this step before sending us the report of HW5, and please include the path of `CBE544` folder in your report or email. You can obtain the path by 

```bash
cd
cd CBE544
pwd
```

<a name='first-time-cees'></a>
____

<a name='testing'></a>

## Making Sure Everything Works ##

Once you are logged into the terminal, run:

```bash
ase-gui
```

and make sure a graphical interface appears. Next, run Python in interactive mode by typing:

```bash
python
```

and make sure the following commands work:

```python
import ase
import numpy
```


