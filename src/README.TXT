##############__HW5 In Course: Distributed and parallel computing__##############

Made by: Tamir Reznik 
	I.D: 308346170

Brief:
------
		In this project I implemented shear sort:
		The input is list of Cuboid's located in text file "cuboid.dat", every line in file represent Cuboid detailes -> id, length, width, height 
		The sort criteria is the volume of each cuboid and incase the volume of 2 cuboid's are the same then it compare the heights
		In this parallel implementation i used mpi library to parallel, reorder procesors in cartesian way and communication between processes
whats in project: 

3 files: 

	- HW5.c				source file to compile and execute (execute as parallel program ->  compile via mpicc -lm and run via mpiexec!)
	
	- HW5.h				header file used in HW5.c (necessary for execute HW5.c)
	
	- cuboids.dat		Example input file with list of 16 cuboid's
	
	- README.TXT 			<---#--->
	
How to run:
	
	two ways to run the file (ubuntu enviorment):
	
	1) use the terminal:
			- mpicc [C_FILE_NAME] -lm -o [EXECUTABLE_NAME] 		|command to compile c files that are not sequential!
			
			- mpiexe -np [#NUMBER_OF_PROCESSES] [EXECUTABLE_FILE_NAME] to run the parallel program file
			
	2) use the eclipse
	
		if not sequential file:
			- Build project 
			- configure execution settings (number of processes, name of executable,arguments, etc.)
			- run
			###NOTE!! - number of cuboid lines in file must be equal to number of processes!!###
			