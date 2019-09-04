#!/bin/sh

#module load openmpi-x86_64
module load openmpi-1.10-x86_64

# Run MPI program through Ethernet eth0

mpirun -np 28 /localdisk/software/hyphy/hyphy-2.2.4/HYPHYMPI USEPATH=/dev/null BASEPATH=/localdisk/software/hyphy/lib/hyphy /localdisk/home/PATH_TO_FOLDER_/SCUEALfullG.bf

