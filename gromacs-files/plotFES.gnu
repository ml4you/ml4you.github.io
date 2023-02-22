set terminal png font "arial,20" size 1200,1200
set output 'FES_rgyr_rmsd.png'
set title "Free Energy Surface"
unset key
set pm3d map
set vi map
set xrange [0.7:1.5]
set yrange [0:0.8]

set xlabel "Rgyr [nm]"
set ylabel "RMSD [nm]"
set cblabel "{/Symbol D}G [kcal/mol]"

splot 'FES_rgyr_rmsd.dat' # Inputdatei