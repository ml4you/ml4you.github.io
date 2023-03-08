#!/bin/bash
FILE=em.gro
if test -f "$FILE"; then
	gmx make_ndx -f em.gro -o index.ndx < Choice_4.txt
else 
   	 	echo "$FILE does not exist."
                exit
fi


FILE1=em.gro
FILE2=topol.top
FILE3=index.ndx
if test -f "$FILE1"; then
       if test -f "$FILE2"; then
 		if test -f "$FILE3"; then
		gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr	
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

FILE=nvt.tpr
if test -f "$FILE"; then
	gmx mdrun -deffnm nvt -v
else 
    echo "$FILE does not exist."
                exit
fi
