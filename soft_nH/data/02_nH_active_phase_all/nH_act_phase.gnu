#
set terminal\
    postscript\
    landscape\
    enhanced\
    color\
    blacktext\
    solid
# enter the name of the model file without ".txt" extension
modelname='nH_act_phase_i90_wd1p5e22'
# enter value of the corresponding nHwd value
nHwd=1.5e+22        
#
set out modelname.'.eps'
#
set size 0.71,0.607
#
set multiplot
#
set lmargin 7
set rmargin 0
set bmargin 0
set tmargin 0
#
set key at 0.5,25.5 center samplen 1 width -1.8 height 0.5 spacing 1.25\
           box lc rgb "black" lw 1
#
set size 0.7,0.5
set origin 0.007,0.095
#
set xrange [-0.05:2.05]
set yrange [21.5:26.3]
#
set xlabel "Orbital phase" font "Times-Roman,17"
set ylabel "Log(n_H [cm^{-2}])" font "Times-Roman,17"
#
set xtics -1,0.2
set ytics 0,1
#
set mxtics 2
set mytics 2
#
#
set arrow 7 from 1,21.5 to 1,26.3 nohead lt 0 lw 1
#
#----------------------------
#
plot \
'nH_act_phase_measured.txt' u 1:(log10($2)) :(log10($2+$3)):(log10($2-$4)) \
          t "n_H measured" w errorbars pt 7 ps 0.9 lc rgb "black" lw 1, \
'' u ($1+1):(log10($2)):(log10($2+$3)):(log10($2-$4)) \
          t"" w errorbars pt 7 ps 0.9 lc rgb "black" lw 1, \
modelname.'.txt' u 1:(log10($2))\
               t"RG wind" w l lt 0 lc rgb "#ff0066" lw 4, \
modelname.'.txt' u ($1+1):(log10($2))\
               t"" w l lt 0 lc rgb "#ff0066" lw 4, \
modelname.'.txt' u 1:(log10($2+nHwd))\
               t"WD wind + RG wind" w l lc rgb "#cc99ff" lw 4, \
modelname.'.txt' u ($1+1):(log10($2+nHwd))\
               t"" w l lc rgb "#cc99ff" lw 4
#
quit
#
