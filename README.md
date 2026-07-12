# EL-Cage
EL-Cage takes trajectory data from molecular simulations and extracts cage-breaking (or hopping) events. Correlated motion can also be detected. The 'EL' refers to the Energy Landscape as the code is particularly effective when applied to locally minimised configurations. 'Cage' refers to cage-breaking, as before contributing to long-term diffusion, particles in solid or solid-like materials must break out of the nearest-neighbour cage formed by neighbouring particles.  

Originally developed for classical models of disordered materials, incuding the binary Lennard-Jones model glass-former [1], Lewis-Wahnström Ortho-Terphenyl (OTP) [2] and van Beest-Kramer-van Santen (BKS) silica [3]. 
FORTRAN and Python versions of this code (previously known as Cage_Break) were made available on the Cambridge Data Repository [4]: https://doi.org/10.17863/CAM.13646  

The FORTRAN code was written by Vanessa Ward (née de Souza) with minor changes and/or comments added by Sam Niblett and Myra Biedermann. The Python code was a re-write of the original FORTRAN code, created by Sam Niblett.

<img width="593" height="159" alt="image" src="https://github.com/user-attachments/assets/f22aa86f-eda5-45c8-9a8e-cad146470104" />  

The images above show a Cage-break or atomic hop in a supercooled binary Lennard-Jones model glassformer [1]. 

EL-Cage is an updated version of this code that allows for greater functionality and can be applied to analyse a wide range of different materials [5]. For simple systems, where there is only one type of ion hop, and self-correlations are only present in the form of back-and-forth correlations, the EL-Cage code provides a complete analysis of the simulation data. For more complex systems and/or for analysing distinct-ion correlations pertinent to conductivity, further post-processing of the extracted ion hopping data is required.

<img width="609" height="202" alt="image" src="https://github.com/user-attachments/assets/5574bc54-1720-4e78-ad34-4d5f3cc08ba6" />  

The images above show the multiple-step ion hops with back-and-forth motion for the Li3OCl anti-perovskite [5].

## Installation
EL-Cage can be installed from source:

```bash
git clone https://github.com/VKWard/EL-Cage.git
cd EL-Cage
make
```
Or the code can be copied, and compiled by simply running 'make' in the directory. By default, the Makefile is configured for the gfortran compiler (https://fortran-lang.org/learn/os_setup/install_gfortran/). 

Useful source files:

- *criteria.h* contains the various parameters used in the definition of cage breaks.<br/>
DISTAA, DISTAB, DISTBB are the distance criteris used to identify when two atoms are nearest neighbours (for A-A, A-B and B-B pairs of atoms, respectively).<br/>
DISTVAA, DISTVAB, DISTVBB are the 'movement cutoffs' described in [1-3].<br/>
REVCUT and SAMECUT are the squared distance cutoffs used to identify reversed cage-breaks (for direct reversals by cage-breaking rearrangements and indirect reversals by non-cage-breaks, respectively).<br/>
NCUTA and NCUTB are the numbers of neighbours which must change for an atom to count as breaking its cage (for A and B atoms, respectively).<br/>

- *io1.f* contains the unit numbers used for reading and writing files. This is useful for working out what the rest of the program is doing.<br/>

- *system.h* contains the various system parameters, including box lengths for the periodic system and the number of atoms of each type.<br/>

- *SETUP* contains information about the number of trajectory configurations. This is used to allocate array dimensions.<br/>

## Execution
A trajectory data file _pathcoords.1_ containing all the (minimised) coordinates of the system is required. Once unzipped, this file is data in the form x1 y1 z1 with no headers. The program uses 'unwrapped' coordinates*.<br/>
System details are added in the header file, system.h, this file should be edited before compilation. EL_Cage currently only uses two types of atom, A and B (although multiple different atoms can be grouped together as A and/or B).<br/>
It is the _criteria.h_ file that determines how ion hops are identified. Ion hops for both A and B species can be identified, using A, B or A and B neighbours as a reference. Once this file has been updated, the code can be compiled and run as follows:
```bash
make
./CageBreak >outfile 
```

*_The additional removebox.f program can be used to convert data into this format, from LAMMPS output - one file containing the coordinates of each configuration, with 9 header lines including timestep, number of atoms and boxlengths. Coordinates are read from oldfile.xyz, and written to newfile.xyz. Using the example data this code can be compiled and run as follows:_
```bash
gfortran -o removebox.o removebox.f
./removebox.o 
```
_newfile.xyz can be copied to the EL-Cage input file, pathcoords.1 <br/>
removebox.f also provides three files _boxlengths_, _nrecs_ and _ionnumbers_, with information about periodic box lengths, the number of records in the trajectory and the numbers of different ions, in the format required for the main EL-Cage code. These lines can be copied into system.h line 1, SETUP line 2 and system.h line 3 respectively.<br/> 
The version of removebox.f here only extracts the first two types of ions/atoms, e.g. For Li3OCl, the Li and O ions remain but the Cl ions are removed._ 

## Output Files

The main information about the ion hops is contained within 3 files:<br/>

- *nnchange.1* contains information about the ion hops identified through nearest neighbour changes. The data are: step number, atom number, reversal count (ony relevant for simple systems), number of neighbour changes, neighbours that changed (a minus sign indicates a lost neighbour, a positive value is a neighbour gained).
- *cbposition.1* contains two lines for each ion hop event, the starting (line 1) and end (line 2) coordinates for each hop. The data are: step number, atom number, x, y, z.
- *rjump.1* contains the hopping distances for each hop. The data are: hop distance, ion number.

The _outfile_ contains a summary of the analysis including the total square diplacement calculated by summing hopping events (Total r2: all steps) and an estimate of the total square displacement when back-and-forth events are excluded (Total r2: irreversible steps).

Example output is shown in the 'output' subdirectory.

## Output without a neighbour analysis

The EL_Cage code can also be used to identify ion hops by the distance moved.<br /> 
The criterion used is DISTVAA in the criteria.h file.<br/> 
This method is used if ABSDIST is set to true in switch.h.<br/>
See the required switch.h file and output for this method in the subdirectory 'output.distonly'.<br/>
The simple distance criterion can work well in ordered systems.


## References
[1] V. K. de Souza and D. J. Wales, J. Chem. Phys., 129, 164507 (2008). _Energy Landscapes for Diffusion: Analysis of Cage-Breaking Processes._  
https://doi.org/10.1063/1.2992128  
[2] S. P. Niblett, V. K. de Souza, J. D. Stevenson and D. J. Wales, J. Chem. Phys. 145, 024505 (2016). _Dynamics of a Molecular Glass Former: Energy Landscapes for Diffusion in Ortho-Terphenyl._  
https://doi.org/10.1063/1.4954324  
[3] S. P. Niblett, M. Biedermann, D. J. Wales and V. K. de Souza, J. Chem. Phys., 147, 152726 (2017). _Pathways for diffusion in the potential energy landscape of the network glass former SiO2._  
https://doi.org/10.1063/1.5005924  
[4] Cambridge Data Repository, Research Data supporting [3]. https://doi.org/10.17863/CAM.13646  
[5]  V. K. Ward and K.E. Johnston, submitted for publication 2026. _A rule that’s made to be broken? Reframing the Arrhenius law for ion transport in solid electrolytes._



