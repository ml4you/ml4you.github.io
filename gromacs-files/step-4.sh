#!/bin/bash
FILE=nvt.tpr
if test -f "$FILE"; then
	gmx mdrun -deffnm nvt -v
else 
    echo "$FILE does not exist."
                exit
fi
