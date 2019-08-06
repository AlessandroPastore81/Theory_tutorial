#Print the first n eigenstates for proton and neutron
#Number of states ns
ns = 5
set terminal pdf       
set output 'wfn.pdf'
set xlabel "r[fm]"
do for [j=0:ns]{
s=j+1
set multiplot
set ylabel "Neutron Wave-function n=".s.''
plot 'neutron_wfs.dat' i j pt 2 lt rgb "red" title 'Neutron Wave function n='.s.'','wfs_test.dat' i j pt 1 lt rgb "purple" title 'Test function n='.s.'', '' u 1:(0) w l
unset multiplot}
reset

set terminal pdf       
set output 'wfp.pdf'
set xlabel "r[fm]"
do for [j=0:ns]{
s=j+1
set multiplot
set ylabel "Proton Wave-function n=".s.''
plot "protons_wfs.dat" i j pt 2 lt rgb "red" title 'Proton Wave function n='.s.''
unset multiplot}
reset

set terminal pdf       
set output 'DensityN.pdf'
set xlabel "r[fm]"
set multiplot
set ylabel "Neutron Particle density [fm^3]"
plot "DensityN.dat" u 1:2 pt 2 lt rgb "red" title 'Neutron Particle Density'
unset multiplot
set multiplot
set ylabel "Neutron Kinetic density [fm^-5]"
plot "DensityN.dat" u 1:3 pt 2 lt rgb "red" title 'Neutron Kinetic Density'
unset multiplot
set multiplot
set ylabel "Neutron Spin-orbit density [fm^-4]"
plot "DensityN.dat" u 1:4 pt 2 lt rgb "red" title 'Neutron Spin-orbit Density'
unset multiplot
reset


set terminal pdf       
set output 'DensityP.pdf'
set xlabel "r[fm]"
set multiplot
set ylabel "Proton Particle density [fm^-3]"
plot "DensityP.dat" u 1:2 pt 2 lt rgb "red" title 'Proton Particle Density'
unset multiplot
set multiplot
set ylabel "Proton Kinetic density [fm^-5]"
plot "DensityP.dat" u 1:3 pt 2 lt rgb "red" title 'Proton Kinetic Density'
unset multiplot
set multiplot
set ylabel "Proton Spin-orbit density [fm^-4]"
plot "DensityP.dat" u 1:4 pt 2 lt rgb "red" title 'Proton Spin-orbit Density'
unset multiplot
reset


set terminal pdf       
set output 'PotentialN.pdf'
set xlabel "r[fm]"
set multiplot
set ylabel "Neutron V(r) Potential"
plot "neutron_potential.dat" u 1:2 pt 2 lt rgb "red" title 'V(r)'
unset multiplot
set multiplot
set ylabel "Neutron Spin-orbit Potential"
plot "neutron_potential.dat" u 1:3 pt 2 lt rgb "red" title 'Vso(r)'
unset multiplot
set multiplot
set ylabel "Neutron Total Potential"
plot "neutron_potential.dat" u 1:($2+$3) pt 2 lt rgb "red" title 'Vtot(r)'
unset multiplot
reset

set terminal pdf       
set output 'PotentialP.pdf'
set xlabel "r[fm]"
set multiplot
set ylabel "Proton V(r) Potential"
plot "proton_potential.dat" u 1:2 pt 2 lt rgb "red" title 'V(r)'
unset multiplot
set multiplot
set ylabel "Proton Spin-orbit Potential"
plot "proton_potential.dat" u 1:3 pt 2 lt rgb "red" title 'Vso(r)'
unset multiplot
set multiplot
set ylabel "Proton Coulomb Potential"
plot "proton_potential.dat" u 1:4 pt 2 lt rgb "red" title 'Vc(r)'
unset multiplot
set multiplot
set ylabel "Proton Total Potential"
plot "proton_potential.dat" u 1:($2+$3+$4) pt 2 lt rgb "red" title 'Vtot(r)'
unset multiplot
reset

set terminal pdf       
set output 'spe.pdf'
set multiplot
set key off
unset xtics
set ylabel "Energies [MeV]"
plot 'neutron_singleparticles.dat' using (-1):($1):(2):(0) with vectors nohead lw 4 lt rgb "red"
unset multiplot
reset

