#!/bin/bash
#! How many whole nodes should be allocated?
#SBATCH --nodes=10
#! How many (MPI) tasks will there be in total? (<= nodes*56)
#SBATCH --ntasks=20
#! How much wallclock time will be required?
#SBATCH --time=1000:00:00
#SBATCH --partition=compute
FILE1=nvt.gro
FILE2=nvt.cpt
FILE3=topol.top 
FILE4=index.ndx
if test -f "$FILE1"; then
       if test -f "$FILE2"; then
 		if test -f "$FILE3"; then
                       if test -f "$FILE4"; then
		                        gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr	
		       else 
   	 		   echo "$FILE4 does not exist."
               		   exit
		       fi
	        else 
                    echo "$FILE3 does not exist."
               	    exit
		fi
       else 
           echo "$FILE2 does not exist."
           exit
       fi
else
    echo "$FILE1 does not exist."
    exit
fi

FILE=npt.tpr
if test -f "$FILE"; then
	gmx mdrun -deffnm npt -v
else
    echo "$FILE does not exist."
    exit
fi
