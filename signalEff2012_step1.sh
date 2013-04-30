#!/bin/sh

#first compile the macro to count events
root -b -l -q signalEff2012_compile.C

#nominal T1bbbb version
#files=`ls /cu4/ra2b/reducedTrees/v68_1/ORIGINALS/reducedTree.*.*T1bbbb*.root`
#pMSSM1 version
files=`ls /cu4/ra2b/reducedTrees/v68_3/ORIGINALS/reducedTree.*.*T5tttt*.root`

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

done

#from cmd prompt:
#[dellcmscornell] ~/work/private/cfaAnalysis/CMSSW_5_2_5/src/NtupleTools/BasicLoopCU $ root -b -l -q 'signalEff2012_T1bbbb_step1.C("/cu4/ra2b/reducedTrees/v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.TTbarJets.root")' > & eventcounts.reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.TTbarJets.root.log &
