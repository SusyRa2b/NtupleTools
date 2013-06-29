#!/bin/sh

#first compile the macro
root -b -l -q signalEff2012_writetxt_compile.C

files=`ls eventcounts.*.T1ttcc.root`

#loop over files
for j in $files;
do 

  filename=$(basename "$j")
  filename="${filename%.*}"
  root -b -l -q "signalEff2012_combinebins.C(\"$j\")" >& eventcounts2x2.$filename.log &

done
