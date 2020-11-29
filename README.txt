SpikeMutator: COVID-19 Spike Mutator & Stability Energy Plotter
(c) 2020 Mohamed R. Smaoui, Hamdi Yahyaoui
Dept of Computer Science, College of Sciences, Kuwait University, Kuwait
mohamed.smaoui@ku.edu.kw
http://spikemutator.com

SpikeMutator is a program for generating PDB structures of desired mutations in the Spike protein of the SARS-CoV-2 virus.
The mutations are applied on the 6VXX PDB template and undergo Molecular Dynamics Energy Minimization to relax any steric clashing.
The program utilizes a dipolar-solvent water model to compute the solvation energy of the mutated structure and computes the Lennard-Jones 
and Coulomb potentials to determine the stability of the generated structures.

SpikeMutator runs as a command-line application on linux machines.


INSTALLATION:

	SpikeMutator utilizes the output of several programs to compute stability. In particular, it requires
	the prior installation of SCWRL4, AQUASOL, PDB2PQR, GROMACS, gfortran, and Python2.7.

	Step 1: 
		The first 3 software (SCWRL4, AQUASOL, and PDB2PQR) can be found in the "extra_tools" directory. To use SCWRL4, you must register 
		for a free academic license at:
		http://dunbrack.fccc.edu/scwrl4/license/index.html

		unzip all 3 zip filest
		$ cd extra_tools
		$ tar -xvf AquaSol_Complexes.tgz
		$ unzip pdb2pqr-1.7.zip
		$ unzip SCWRL4.zip
	

	Step 2:
	
		Install the AquaSol program that calculates Solvation energy execute the following 2 commands
		
		$ cd AquaSol_Complexes
		$ cd Aqua; sudo make clean
		$ cd ../src; sudo make clean
		$ cd ../Aqua; sudo make
		$ cd ../src; sudo make
	

	Step 3:
		
		Install SCWRL4.0 (after obtaining the license from http://dunbrack.fccc.edu/scwrl4/license/index.html)

		$ cd SCWRL4.0
		$ ./install_Scwrl4_Linux

		Make sure you manually enter the full path to the SCWRL4.0 directory and update your license information

	
	Step 4:
		
		Open the SpikeMutator.py file and edit lines 4-8 to update the paths of SCWRL4, AQUASOL, PDB2PQR, and GROMACS.

		If you don't have GROMACS installed, you can download the latest version of GROMACS from: 
		http://www.gromacs.org/Downloads
		


USAGE:

	USAGE 1: Mutate a single residue

		$ python SpikeMutator.py -p position -m AAmutation
			-p position: the residue (amino acid) position you want to mutate
			-m AAmutation: the amino acid mutation you want to replace position p with.
			   This is a letter. ex. Alanine is A

		This commands performs a single point mutation on the spike protein trimer structure, generates a pdb file of the mutant trimer, and returns
		the Coulomb, LJ, and Solvation energies of the mutant spike. 


	USAGE 2: Mutate multiple residues simultaneously

		$ python SpikeMutator.py -p p1,p2,p3,...,pn -m m1,m2,m3,...,mn 
			-p p1,p2,p3,...,pn: A set of integer positions representing the amino acid 
			   positions you wish to mutate
			-m m1,m2,m3,...,mn: A set of amino acid point mutations you wish to replace
			   p1,p2,p3,...,pn with, respectively. Each mutation should be a letter. ex Alanine is A

		This commands performs simultaneous multiple single point mutations on the spike protein structure, and
		returns the Coulomb, LJ, and Solvation energies of the mutant spike.  

