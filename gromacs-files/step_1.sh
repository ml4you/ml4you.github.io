#!/bin/bash
for i in 1 2 
do
cd $i
FILE=complex.gro
if test -f "$FILE"; then
	gmx editconf -f complex.gro -o newbox.gro -bt dodecahedron -d 1.0
else 
   	 	echo "$FILE does not exist."
                exit
fi

FILE=newbox.gro
if test -f "$FILE"; then
	gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
else 
   	 	echo "$FILE does not exist."
                exit
fi
FILE1=ions.mdp
FILE2=solv.gro
FILE3=topol.top
if test -f "$FILE1"; then
       if test -f "$FILE2"; then
 		if test -f "$FILE3"; then
			gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
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

FILE1=ions.tpr
FILE3=solv_ions.gro
FILE3=topol.top
if test -f "$FILE1"; then
       if test -f "$FILE2"; then
 		if test -f "$FILE3"; then
			gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral < Choice.txt
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
FILE1=em.mdp 
FILE3=solv_ions.gro
FILE3=topol.top
if test -f "$FILE1"; then
       if test -f "$FILE2"; then
 		if test -f "$FILE3"; then
			gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
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
FILE=em.tpr
if test -f "$FILE"; then
	gmx mdrun -v -deffnm em
else 
   	 	echo "$FILE does not exist."
                exit
fi	
FILE=nad.gro
if test -f "$FILE"; then
	gmx make_ndx -f nad.gro -o index_nad.ndx < Choice_2.txt
else 
   	 	echo "$FILE does not exist."
                exit
fi

FILE1=nad.gro
FILE2=index_nad.ndx
if test -f "$FILE1"; then
 		if test -f "$FILE2"; then
			gmx genrestr -f nad.gro -n index_nad.ndx -o posre_nad.itp -fc 1000 1000 1000 < Choice_3.txt
		else 
   	 		echo "$FILE2 does not exist."
               		 exit
		fi
else
echo "$FILE1 does not exist."
               		 exit
fi
