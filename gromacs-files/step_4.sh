#!/bin/bash
FILE=nvt.tpr
if test -f "$FILE"; then
	gmx mdrun -deffnm nvt -v
else 
    echo "$FILE does not exist."
                exit
fi

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

FILE1=npt.gro
FILE2=npt.cpt
FILE3=topol.top 
FILE4=index.ndx
if test -f "$FILE1"; then
       if test -f "$FILE2"; then
 		if test -f "$FILE3"; then
                       if test -f "$FILE4"; then
		                       gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr	
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

FILE=md.tpr
if test -f "$FILE"; then
	gmx mdrun -deffnm md -v
else
    echo "$FILE does not exist."
    exit
fi 
   	 	echo "$FILE does not exist."
                exit

