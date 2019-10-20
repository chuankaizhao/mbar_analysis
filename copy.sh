#!/bin/bash

addnum=3
maxnum=40
prefix=PYL2-HAB1-APO-REUS-SEP

for j in `seq 0 ${addnum}`
do
  mkdir ${j}
  cd ${j}
  cp ../../../../../more_sampling/03-reus-equilibrium/02-equilibration/window_${j}/*.colvars.traj ${prefix}.job0.${j}.colvars.traj
  cp ../../../../../more_sampling/03-reus-equilibrium/02-equilibration/window_${j}/*.log          ${prefix}.job0.${j}.history
  linebreak=`awk "/  500000  / {print NR}" ${prefix}.job0.${j}.colvars.traj`
  total=`wc -l < ${prefix}.job0.${j}.colvars.traj`
  tail -$(( ${total} - ${linebreak} )) ${prefix}.job0.${j}.colvars.traj > ${prefix}.job0.${j}.sort.colvars.traj
  grep "ENERGY:" ${prefix}.job0.${j}.history | tail -12000 > temp.history
  awk '{print $2, $9, $13, $14}' temp.history > tempNew.history
  awk 'NR % 2 == 0' tempNew.history > ${prefix}.job0.${j}.sort.history
  rm temp.history tempNew.history
  cd ..
done

for j in `seq 0 ${maxnum}`
do
id=$(( $j + $addnum + 1 ))
mkdir ${id}
cp ../${j}/${prefix}.job0.${j}.sort.colvars.traj ${id}/${prefix}.job0.${id}.sort.colvars.traj
cp ../${j}/${prefix}.job0.${j}.sort.history ${id}/${prefix}.job0.${id}.sort.history
done
