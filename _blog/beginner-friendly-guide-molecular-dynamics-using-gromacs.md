---
title: 'Beginner-friendly guide to molecular dynamics using Gromacs'
date: 2023-02-22
permalink: /blog/2023/02/beginner-friendly-guide-molecular-dynamics-using-gromacs
excerpt_separator: <!--more-->
toc: true
tags:
  - molecular dynamics
  - drug design
  - career
---

This tutorial will guide you step-by-step through the process of setting up and running a molecular dynamics simulation. It's assumed that you have a basic understanding of working with the command line, including basic commands and file navigation. Ideally, you should also have some experience using bash on a Linux system, since many of the instructions and commands in this tutorial are specific to this environment. Whether you're a beginner or an experienced user, this tutorial will provide you with valuable insights and practical knowledge to help you successfully conduct your molecular dynamics simulation.
<!--more-->

## 1.	Introduction
In this tutorial, we will be using a small protein called alpha-synuclein, which is involved in Parkinson's disease. This protein aggregates when in the presence of membranes and forms large clumps called Lewy bodies which are involved in neurodegeneration. While larger proteins would take too much time to simulate for a short tutorial session, we have chosen to use alpha-synuclein monomer (PDB ID: 1XQ8) as it is a small protein formed by 140 amino. For more information on alpha-synuclein, please refer to the [research paper](https://www.jbc.org/article/S0021-9258(19)62979-0/fulltext/).

| ![image1xa8](https://user-images.githubusercontent.com/7014404/220604350-7025c221-8108-4d86-8aaf-05ae21d9bdca.png) |
|:--:|
| <b>Figure 1: 3D structure of alpha-synuclein monomer formed by three regions: N-terminal (green), NAC region (yellow), and C-terminal (red).</b> |


## 2. Files 
All files you will need during this tutorial can be downloaded here.
The purpose of each file will be explained more closely throughout the tutorial.


| File  | Description  |
| ------------- | ------------- |
| [1xq8.pdb](https://github.com/yboulaamane/yboulaamane.github.io/blob/master/gromacs-files/1xq8.pdb)  | Alpha-synuclein structure file in PDB format  |
| [ions.mdp](https://github.com/yboulaamane/yboulaamane.github.io/blob/master/gromacs-files/ions.mdp) | Parameter file for ion addition  |
| [em.mdp](https://github.com/yboulaamane/yboulaamane.github.io/blob/master/gromacs-files/em.mdp) | Parameter file for energy minimization  |
| [eq.mdp](https://github.com/yboulaamane/yboulaamane.github.io/blob/master/gromacs-files/eq.mdp) | Parameter file for equilibration |
| [md.mdp](https://github.com/yboulaamane/yboulaamane.github.io/blob/master/gromacs-files/md.mdp) | Parameter file for production MD |
| [distancemap.sh](https://github.com/yboulaamane/yboulaamane.github.io/blob/master/gromacs-files/distancemap.sh) | Bash script for the calculation of distance maps |
| [generateFES.py](https://github.com/yboulaamane/yboulaamane.github.io/blob/master/gromacs-files/generateFES.py) | Python script for calculating free energies |
| [plotFES.gnu](https://github.com/yboulaamane/yboulaamane.github.io/blob/master/gromacs-files/plotFES.gnu) | Gnuplot file for plotting free energies  |



## 3.	Setup
This section will explain what is needed to set up an MD simulation with Gromacs, starting with protein coordinates in a PDB file.
The figures throughout this tutorial give an overview of the current position in the molecular dynamics workflow. They will expand by clicking on them.

### 3.1.	Introduction
It takes two ingredients for Gromacs to understand a molecular system: atom coordinates and topology information. The topology file (typically ending in *.top) contains information on how many molecules of which kind a system contains, of what kind their interactions are and the parameters that are used to compute the energy and forces. Gromacs can generate topologies for you. To define a topology, a force field needs to be specified. A force field is a collection of well balanced parameters for the energy and force computation.

| ![image](https://user-images.githubusercontent.com/7014404/220606936-54ec8647-707a-4d6d-bf82-0ed39f1937d5.png) |
|:--:|
| <b>Figure 2: Force field terms present in classical molecular mechanics force fields: harmonic bond length and angle bending potentials, periodic torsion angle potentials, Lennard-Jones potentials to describe van der Waals interactions and coulombic interactions. In addition, the solation effects can be included included implicitly through a continuum solvent model. (Figure taken from Wikipedia, Creative Commons Attribution-Share Alike 3.0 Unported license)</b> |

Several different families of force fields exist, differing in the way their parameters were derived. Within each family, often multiple versions exist, each of which includes corrections to the force field parameters, hopefully making the force field more accurate (within its domain of applicability). But newer versions do not necessarily perform better than their predecessors in all situations. Which force field is good for your purpose depends strongly on your system. For protein systems, the AMBER and CHARMM force fields are very often a good bet. If you are interested, you can have a look at this historical overview of force field evolution and a review of modern force fields. Protein force fields only contain parameters for standard amino acids, ions and sometimes for nucleic acids and some sugars. If your system contains ligands or unnatural amino acids you will need to find out how to parameterize molecules in accordance your specific choice of force field. In addition to the force field, a water model is also needed. Water is a tough molecule to simulate in classical MD simulations. This is due to the electronic structure of water that gives rise to complex interactions. Several different water models exist. A good overview of water models can be found on the web. It is usually a good idea to stick to a water model that is recommended for your choice of force field or a compatible one.

### 3.2.	PDB formatted input
Download the file 1xq8.pdb from the Files section of this tutorial. Open the file with VMD to have a look at the protein. Now generate a topology with the command below. Choose Amber99sb-ildn as a force field together with the TIP3P water model. While the TIP4P-EW water model performs better in combination with the Amber force fields, it consists of four instead of three particles. To limit the simulation size we will use the TIP3P water model here.

```
gmx pdb2gmx -f 1xq8.pdb -o 1xq8.gro -p 1xq8.top -i 1xq8_posre.itp -ignh
```
All Gromacs commands work similar. The main program is called gmx and starts every command. It will help you to remember that gmx help commands shows you a list of all available commands with a short description for each of them, while gmx help <command;> gives a more detailed description of the tool whose name you substituded for <command;>. (Gromacs versions before 5.0 consisted of many programs that have now been merged into gmx. For compatibility reasons it is therefore still possible to omit gmx and directly start with e.g., pdb2gmx, but this is not guaranteed to work in the future). The command generated three files. 1xq8.gro contains residue and atom names and the coordinates of all atoms. It could in addition contain velocities for the atoms, but at this point we haven’t generated any yet. The second file is the topology file (1xq8.top). Have a look at its contents using the less command. After some comments, this file first includes the parameters of the force field in use:
```
; Include forcefield parameters
#include "amber99sb-ildn.ff/forcefield.itp"
```

An include directive (#include) instructs Gromacs to substitute the contents of the specified file at that particular position in the topology file. The forcefield.itp file contains the force field parameters of the Amber99sb-ildn force field, which is contained in the Gromacs installation folder. The next section defines molecule types to make them known to Gromacs.
```
[ moleculetype ]
; Name nrexcl
Protein 3
```
Right now we only have a single protein in our system. If there were more, chain identifiers would be added to the molecule names. The atom section that follows lists all atoms of all residues with some properties like mass and charge. Then some sections follow that specify the bonds, angles, dihedrals and pair (also known as VdW) interactions.
The third file we have generated is 1xq8_posre.itp. It contains position restraints for all heavy atoms of the protein, i.e., harmonic potentials that couple the atoms to reference positions. This will be necessary during preparatory steps for the MD simulation to let the water equilibrate around the protein without the protein itself moving significantly, so it cannot get deformed by unequilibrated water molecules. After that, topologies for the water molecules and ions are included as well.
Finally the system is given a name and the number of each molecular species in the system is listed:
```
[ system ]
; Name
Protein in water
[ molecules ]
; Compound #mols
Protein 1
```
The latter table will be updated when we add water and ions. The rest of the topology file usually remains constant beyond this point. Consult the Gromacs manual for more information on the layout of topology files. Note that different pH conditions can only be considered in a static fashion, that is explicit protonation states need to be specified during this initial setup and remain constant during the simulation. The protonation states of most residues will be set corresponding to neutral pH, while the protonation states of histidines are chosen by pdb2gmx to optimize hydrogen bonding networks.
We have now generated a topological description of our protein.

### 3.3. Setting up the simulation box
The next step is to create a box around our protein that we can later fill with solvent. Gromacs can generate boxes of various shapes: triclinical cells, cubes, dodecahedra and octahedra. The triclinical cell has box vectors at right angles, but with different length, while the cubical box has in addition box vectors of equal length. The solute will be centerd in the simulation box with a minimum distance to any wall of 10 Ångstrom (1.0 nm). As the speed of the simulation decreases approximately linearly with the number of particles, we want the box to be as small as possible, but at the same time keep the solute away from its periodic images at all times. The triclinical cell does not work well, as its rectangular shape is problematic for non-globular proteins that might rotate, unless the center of mass rotation is taken care of. We will choose a dodecahedral cell, as it is close in shape to a sphere and minimizes the volume needed to keep the minimum solute to wall distance:
```
gmx editconf -f 1xq8.gro -o 1xq8_box.gro -bt dodecahedron -d 1.0
```
This updated the information about the box that contains our system in the final line of the *.gro file. (You can navigate to the end of a text file in the less program with the G key. g brings you back to the top) The numbers in this final line correspond to the box vector components (in nm) in the order: v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y). The latter 6 numbers may be omitted for rectangular boxes and Gromacs only supports boxes with v1(y)=v1(z)=v2(z)=0. Generate some of the other box shapes too and have a look at them with VMD by loading the corresponding *.gro file and typing pbc box in the VMD console. The dodecahedral and octahedral boxes appear not as you might have expected. The reason is that the boxes are always represented in a general, brick shaped unit cell representation for efficiency. Consult the section on periodic boundary conditions in the Gromacs manual to get an idea of how the different representations are transformed into each other. Learn more about the different box types (popup) offered by Gromacs.

### 3.4. Solvation
Now the box can be filled with solvent molecules. We have selected the TIP3P water model before, so we need to provide the coordinates of a three-particle water model. As there are no TIP3P coordinates in the Amber force field directories, we will use the coordinates from spc216.gro. This file will be found automatically by Gromacs, without specifying a detailed path. Note that we are still using the correct parameters of the TIP3P water model and only borrowing the coordinates from 216 water molecules, which were produced with the SPC water model.
```
gmx solvate -cp 1xq8_box.gro -cs spc216.gro -o 1xq8_solvated.gro -p 1xq8.top
```
Have a look again at the output file in VMD. The solvent box looks brick-shaped for the reason given in the section on box generation. We will later show how to transform the box into its proper octahedral/dodecahedral representation.

### 3.5. Adding ions
It is also possible to add some ions to mimic simulation conditions with specified salt concentrations. However, ions are only approximated as charged point particles and do not possess more detailed electronic structure, which might be relevant for ion binding sites and transition metals ligated to proteins. Though for NaCl the simple point charge description is usually sufficient. We add enough ions to neutralize the system and add an additional concentration of 100 mM NaCl.
For doing so, we need an *.mdp file that contains molecular dynamics parameters. In general, an *.mdp file contains information on the algorithms and options to use during a simulation with Gromacs. It is parsed by the Gromacs preprocessor (grompp) and combined with initial coordinates and velocities into a *.tpr file. *.tpr files contain the complete information that is needed for an energy minimization or any kind of MD simulation with gromacs and is the only required input file for Gromacs' computing engine mdrun.
So why do we need an *.mdp file at this step, as we only want to add ions to our system? Gromacs tries to insert the ions at optimal positions that minimize the energy of the system. Therefore, it needs to know which algorithms and parameters to use to compute the energy.
Have a look at the file ions.mdp (popup).Parse the contents of ions.mdp:
```
gmx grompp -f ions.mdp -po ions.out.mdp -c 1xq8_solvated.gro -p 1xq8.top -o ions.tpr -maxwarn 1
```

As a net charge of the system different from zero could lead to artefacts in combination with the Ewald method, a warning is issued if this should be the case and no files are produced by grompp. However, since the system will be neutralized at the next step, we can ignore this warning which is done with the flag -maxwarn 1. The resulting *.tpr file will be fed into a tool that acually generates the ions:
```
gmx genion -s ions.tpr -o 1xq8_ions.gro -p 1xq8.top -pname NA -nname CL -neutral -conc 0.1
```
When Gromacs asks which molecules to replace with ions, select the SOL group, which contains all atoms of the solvent (water) molecules. The ions will then be substituted for water molecules at positions that minimize the overall electrostatic energy.

## 4. Energy minimization
During the setup of our system, we might have introduced unnatural stress, for example by placing two atoms accitentally too close to each other. This might result in large forces acting on the particles that would blow up our system should we simply start a regular molecular dynamics simulation at this point. To cope with this problem, we will perform an energy minimization in which we relax the system to the closest local energy minimum. We will again first generate a run-input or .tpr file with grompp to combine the initial setup with the algorithms and parameters from an .mdp file:
```
gmx grompp -f em.mdp -po em.out.mdp -c 1xq8_ions.gro -p 1xq8.top -o em.tpr
```
The result is a .tpr file that is the only file needed to start our energy minimization. Some new parameters are introduced in the file em.mdp. Examine it and also have a look again at the mdp section of the Gromacs Online Reference to familiarize with all the parameters.
We now use a steepest descent algorithm to locate the closest local minimum and set some convergence criteria. In addition, we set the coordinate output frequency to every 10 steps with nstxout = 10. The energy minimization can be run with:
```
gmx mdrun -v -deffnm em >& em.out &
```
The -v switch will make mdrun create verbose output and report the computation progress. -deffnm defines a default filename prefix that is used for all in- and output files. The corresponding file endings will be appended for the different in- and output files, so they do not need to be specified individually. You will obtain several output files. The file em.log contains detailed information on several energy terms throughout the simulation and a performance analysis at the end. The file em.gro contains the coordinates of the final configuration of the system. The file em.trr contains coordinates (and optionally velocities in case of MD) and forces for each output frame. The em.edr file contains information on different energy terms and, if applicable for the simulation, information on the temperature, pressure, volume and some other quantities during the simulation. You can extract this information with:
```
gmx energy -f em.edr -o em_Epot.xvg
```
Follow the instructions to select the potential energy. You can plot the graph with the program xmgrace:
xmgrace em_Epot.xvg
You should see a monotonically decaying graph as the potential energy gets minimized. You can also follow the equilibration process with VMD. It is a good idea to ensure that all molecules are complete in the simulation box and do not cross the periodic boundaries to stick out on one end and enter again on the opposite. This would make a huge mess in VMD. You can transform the trajectory with:
```
gmx trjconv -f em.trr -s em.tpr -o em_centered.xtc -ur compact -pbc mol -center
gmx trjconv -f em.trr -s em.tpr -o em_centered.pdb -ur compact -pbc mol -center -dump 0
```
The -dump 0 flag creates a pdb file that contains the coordinates of the protein at time t=0. Select the protein group for centering and the whole system for output. The trajectory will be converted into another file format (*.xtc) that has less precision, but is much tighter compressed, resulting in smaller files. Now load the system into VMD with:
```
vmd em_centered.pdb em_centered.xtc
```
This will load the final configuration stored in em.gro. The information on bonds are automatically generated by VMD. You can optionally delete the coordinates by right clicking on the molecule in the main window and deleteing the single frame. This will make the molecule dissapear in the display window, but the information on atom types and bonds is still retained internally in VMD. This has the advantage that we will not see a jump from the final configuration (as loaded from em.gro) to the first frame of the trajectory. The trajectory can be loaded by selecting File -> New Molecule from the VMD main window. First select the already loaded molecule from the drop-down menu. Then browse to the trajectory file em_centered.xtc. The filetype will be detected automatically. Hit the load button.

## 5. Equilibration
Our energy minimized system now needs an equilibration MD simulation to further relax the system. The water molecules that we created as a crystal lattice can then rearrange around the solute. As we do not want the protein to move significantly during the equilibration phase, we will restrain it with harmonic forces to its initial position. This is done by defining the macro POSRES in the file eq.mdp and making use of 1xq8_posre.itp, which is achieved by the following statement in 1xq8.top:
```
; Include Position restraint file
#ifdef POSRES
#include "1xq8_posre.itp"
#endif
```
This causes grompp to insert the harmonic restraining forces into the *.tpr file. In the parameter file (eq.mdp), we choose the integrator md, which is acutally a leap-frog verlet integrator to compute the position updates at every timestep from the previous timestep. We set the timestep length to 2 femtoseconds. In addition to the parameters you already know from the energy minimization, we now use a thermostat and a barostat to make the temperature and pressure fluctuate around room temperature and atmospheric pressure. We instruct grompp to generate new velocities from a Boltzmann distribution at room temperature. Finally, we set constraints on all bonds. This will ensure that bonds will not deviate too much from their equilibrium positions, even if we increase the timestep beyond the point where it can resolve bond vibrations properly. As we do not expect bond vibrations to influence the overall dynamics of our system significantly, this speeds up the simulation quite a bit without sacrificying much accuracy. Consult the Gromacs manual for a detailed explanation. Generate the run input file:
```
gmx grompp -f eq.mdp -po eq.out.mdp -c em.gro -r em.gro -t em.trr -p 1xq8.top -o eq.tpr
```
Start the equilibration with:
```
gmx mdrun -v -deffnm eq >& eq.out &
```
The & at the end ensures that the command, i.e., the simulation will be run in the background. All output is written in the file eq.out. While the simulation is running you can monitor the progress with the command:
```
tail -f eq.out
```
This will display the end of the file and also print more data as the file grows. The same works as well for the log file (eq.log). Once the simulation is done, you should check that the results are not unphysical. Despite visual inspection you can again have a look at the energy output (eq.edr) with gmx energy to make sure that, for example, the temperature fluctuates around the desired value and the density is not decreasing, which would indicate that your system is blowing up. Once the simulation is done, convert it again to make all molecules whole and transformed to the central simulation box:
```
gmx trjconv -f eq.xtc -s eq.tpr -o eq_centered.xtc -ur compact -pbc mol -center
gmx trjconv -f eq.xtc -s eq.tpr -o eq_centered.pdb -ur compact -pbc mol -center -dump 0
```
The -dump 0 flag creates a pdb file that contains the coordinates of the protein at time t=0. Open the trajectory in VMD as you did after the energy minimization. You should see a protein that is hardly moving, at least not changing its overall shape, and freely moving water molecules and ions. To get a little more used to VMD try the following: Open the Representations menu that you find in the Graphics menu of the VMD Main Window. Change the Selection to protein and the Drawing Method to Licorice. Now switch to the Trajectory tab and set the Trajectory Smoothing Window Size to 5. In effect, each frame you see in the Graphics window is now the average over 5 frames and the protein movement will appear much smoother. Now create a new Representation. Enter the selection text same residue as (within 5 of protein). Set the Drawing Method to Lines and the Trajectory Smoothing Window Size again to 5 or higher. You can now observe the water molecules that are initially close to the protein (in a 5 Ångstrom shell around it) diffusing away. Alternatively you can set the option Update Selection Every Frame in the Representations window’s Trajectory tab. This will show you the water molecules that are close to the protein in each frame. If you are certain that your equilibration looks fine, move on to the production MD stage.

## 6. Production MD
Now we will set up and run the unrestrained production MD simulation, which we will analyse more thoroughly afterwards. You already know how to create the run input file, so this is slowly becoming a routine. The only things that we change are the removal of position restraints, usage of different thermo- and barostats, and increasing the simulation length to 1 ns. The reason we now use the Nose-Hoover thermostat and Parrinello-Rahman barostat is that they more accurately sample the isothermal-isobaric ensemble We did not use them for equilibration as they introduce large oscillations around their target values (T and p) if the system starts far away from equilibrium. Again, have a look at the mdp file (md.mdp) and familiarize yourself with the parameter settings. Run grompp:
```
gmx grompp -f md.mdp -po md.out.mdp -c eq.gro -t eq.trr -p 1xq8.top -o md.tpr
```
Followed by mdrun:
```
gmx mdrun -v -deffnm md >& md.out &
```
