# ohess
Implementation of OHESS for calculating elastic properties in VASP

OHESS: Optimized High-Efficiency Strain matrix Sets 
Formalism (and definitions for many variable names) given in 
'High-efficiency calculation of elastic constants enhanced by the optimized strain - matrix sets' - Zhong-Li Liu, 2020

Creating the strain sets is done in the c++ code, the python code is a set of helper functions for post-processing (not necessary to apply OHESS).  

Instructions:
(1) Build the cpp file (on linux, you might use 'g++ -o makeStrainSets -std=c++11 MakeStrainSets.cpp')
(2) Create directories ("./ makeStrainSets -i CONTCAR -g G -d D -n N")
  CONTCAR = the (VASP-5+) POSCAR-style file to apply OHESS to
  G = the crystal group index (0 <= G < 9)
  D = delta max (0.0 < D <= about 0.01 to reasonably guarentee that you stay within the elastic regime)
  N = the number of deltas to aply for each deformation matrix.  (2 <= N <= 5 for reasonable compromise between accuracy and efficiency)
(3) Run VASP with appropriate INCAR tags, making sure that ISIF = 2 (fixed cell, write stress components)

* Mapping of crystal group names to crystal group indices:
* 0 = Cubic
* 1 = Hexagonal
* 2 = Rhombohedral / Trigonal I
* 3 = Rhombohedral / Trigonal II
* 4 = Tetragonal I
* 5 = Tetragragonal II
* 6 = Orthorhombic
* 7 = Monoclinic
* 8 = Triclinic

Comments:
-To get a consistent set of k-points, you should break symmetry before applying OHESS strains to a material.
-Relaxing the initial structure to a strict(er) force tolerance allows one to choose delta max to be smaller (and thus can save computation time)
 I've found that getting the intial forces to 1 meV/angstrom allows delta max = 0.002 (with N=5).  Especially consider this if you plan on doing phonon
 calculations anyway
-Restarting all strained jobs from the initial wavefunctions and charge densities is safe and efficient

To use postprocessing:
---shell script---
(1) Don't mess with the names of the created directories
(2) Don't delete POSCAR, CONTCAR, or OUTCAR from anywhere
(3) Use writeOutput.sh "./writeOutput.sh".  From the current working directory down, the only directories that exist should be ones that were made with the ohess code
    (unless they don't contain POSCAR/CONTCAR and OUTCAR, in which case it isnt a problem)
(4) 'output.csv' will be created, which is the file that the python code is designed to read
---python---
(5) Make a .py file and include the header file as an import.
(6) The main function is used like so: sets = headerOhessElas.GetOhessDict(infileLoc="output.csv", indexToMakeKey=3)
    in which sets is a dict of ohessSets whose keys are the indexToMakeKey'th (3rd in this case) portion of the absolute path to a given deformation set
    
      ex: you've done OHESS on TiN and Ti2N.  Your directory paths are /home/calcs/ohess/TiN/... and /home/calcs/ohess/Ti2N/...
          Setting indexToMakeKey=4 would give you a dict of 2 ohess sets with keys=[TiN, Ni2N].
          Setting indexToMakeKey=3 would give you a dict of 1 ohess set with key=TiN (and a bunch of garbage results + possibly a code crash).
        
   If you are only looking at one material, indexToMakeKey isn't that important.  Your ohess data will be the only value in the dict.  
(7) Each ohess instance has two important properties: The Cij matrix (ohess.cijMatrix) and the maxErr (ohess.maxErr).  Both have GPa as units
    -The Cij matrix is filled in thorugh symmetry considerations
    -The max error is an estimate based off of how far the initial (theoretically zero-strain) model's stress is from zero.  Your calculations are no more accurate
     than this value.  

Comments:
-The Sij matrix is not calculated here, but it can be found by inverting the Cij matrix.  
-The Cij matrix is zero-indexed, so to get C23 you'd do something like ohess.cijMatrix[1][2].  



