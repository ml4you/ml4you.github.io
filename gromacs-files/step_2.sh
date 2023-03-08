#!/bin/bash
FILE=topol.top
if test -f "$FILE"; then
	cat topol.top > topol_backup.top
else 
   	 	echo "$FILE does not exist."
                exit
fi


FILE=topol.top
if test -f "$FILE"; then
	python3 ADD.py topol.top > topol_2.top
fi


FILE=topol_2.top
if test -f "$FILE"; then
	mv topol_2.top  topol.top 
else 
   	 	echo "$FILE does not exist."
                exit
FILE=em.gro
if test -f "$FILE"; then
	gmx make_ndx -f em.gro -o index.ndx < Choice_4.txt
else 
   	 	echo "$FILE does not exist."
                exit
fi
fi
