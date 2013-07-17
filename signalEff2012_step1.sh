#!/bin/sh

#first compile the macro to count events
root -b -l -q signalEff2012_compile.C

#nominal T1bbbb version
#files=`ls /cu4/ra2b/reducedTrees/v68_1/ORIGINALS/reducedTree.*.*T1bbbb*.root`

#files=`ls /cu4/ra2b/reducedTrees/v68_3/ORIGINALS/reducedTree.*.*T5tttt*.root`
#files=`ls /cu6/joshmt/reducedTrees/v68_3_pawandeep/reducedTree.*.*T1t1t*.root`
#files=`ls /cu4/ra2b/reducedTrees/v68_3/ORIGINALS/reducedTree.*.*T7btw*.root`
#files=`ls /cu5/joshmt/reducedTrees/v68_4_t14t/reducedTree.*.*T1tttt*.root`
#files=`ls /cu5/joshmt/reducedTrees/v68_6/reducedTree.*.*T1ttcc*.root`
files=`ls /cu6/joshmt/reducedTrees/v69_2/reducedTree.*.*T1bbbb*.root`

#must not be in afs
#####these were for cases where there are many input files
#outdirlog=/home/joshmt/logfiles_6June2013
#outdirroot=/home/joshmt/rootfiles_6June2013

#mkdir -p $outdirlog
#mkdir -p $outdirroot
####end of many input file stuff

#count events in each file
for j in $files;
do 

  filename=$(basename "$j")
  filename="${filename%.*}"
#event counts in all bins
  root -b -l -q "signalEff2012_T1bbbb_step1.C(\"$j\")" >& eventcounts.${filename}_a.log &
#joining b-tag bins
  root -b -l -q "signalEff2012_T1bbbb_step1b.C(\"$j\")" >& eventcounts.${filename}_b.log &
#with PU systematic shift (nominal sample only)
  root -b -l -q "signalEff2012_T1bbbb_step1c.C(\"$j\")" >& eventcounts.${filename}_c.log &
#with pdf variations (nominal only)
  root -b -l -q "signalEff2012_T1bbbb_step1d.C(\"$j\")" >& eventcounts.${filename}_d.log &
#isr +1 (nominal sample only)
  root -b -l -q "signalEff2012_T1bbbb_step1e.C(\"$j\")" >& eventcounts.${filename}_e.log &
#isr -1 (nominal sample only)
  root -b -l -q "signalEff2012_T1bbbb_step1f.C(\"$j\")" >& eventcounts.${filename}_f.log &

  wait
#wait for all of the above to finish before moving on; fancier behavior would be nice but would take new tools (gnu parallel?)

#special hack for the case where a made a gazillion output files
#  tar cvfz $outdirlog/${filename}_logs.tgz ev*log
#  rm ev*log
#  tar cvfz $outdirroot/${filename}_root.tgz ev*root
#  rm ev*root
####end of special hack

done

#from cmd prompt:
#[dellcmscornell] ~/work/private/cfaAnalysis/CMSSW_5_2_5/src/NtupleTools/BasicLoopCU $ root -b -l -q 'signalEff2012_T1bbbb_step1.C("/cu4/ra2b/reducedTrees/v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.TTbarJets.root")' > & eventcounts.reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.TTbarJets.root.log &
