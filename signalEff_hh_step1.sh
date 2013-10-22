#!/bin/sh

#first compile the macro to count events
root -b -l -q signalEff_hh_compile.C

files=`ls /cu6/joshmt/reducedTrees/v71_5b/reducedTree.*.SMS-TChiHH*.root`

#JES0_JER0_PFMETTypeI_METunc0_PUunc0_hpt20

for j in $files;
do 

  filename=$(basename "$j")
  filename="${filename%.*}"

#for all samples
#event counts in all bins
  root -b -l -q "signalEff_hh_step1.C(\"$j\",\"nominal\")" >& eventcounts.${filename}_a.log &
#joining b-tag bins
  root -b -l -q "signalEff_hh_step1.C(\"$j\",\"joinbtag\")" >& eventcounts.${filename}_b.log &

#find the "nominal" sample and do extra event counting there
  if [[ $j = *JES0_JER0_PFMETTypeI_METunc0_PUunc0_hpt20* ]]
    then
#with PU systematic shift (nominal sample only)
      root -b -l -q "signalEff_hh_step1.C(\"$j\",\"pushift\")" >& eventcounts.${filename}_c.log &
#with pdf variations (nominal only)
      root -b -l -q "signalEff_hh_step1.C(\"$j\",\"pdfs joinbtag\")" >& eventcounts.${filename}_d.log &
#isr +1 (nominal sample only)
      root -b -l -q "signalEff_hh_step1.C(\"$j\",\"isrup joinbtag\")" >& eventcounts.${filename}_e.log &
#isr -1 (nominal sample only)
      root -b -l -q "signalEff_hh_step1.C(\"$j\",\"isrdown joinbtag\")" >& eventcounts.${filename}_f.log &

#btag +1 (nominal sample only)
#      root -b -l -q "signalEff_hh_step1.C(\"$j\",\"btagup\")" >& eventcounts.${filename}_e.log &
#btag -1 (nominal sample only)
#      root -b -l -q "signalEff_hh_step1.C(\"$j\",\"btagdown\")" >& eventcounts.${filename}_f.log &

  fi

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
#[dellcmscornell] ~/work/private/cfaAnalysis/CMSSW_5_2_5/src/NtupleTools/BasicLoopCU $ root -b -l -q 'signalEff_hh_step1.C("/cu4/ra2b/reducedTrees/v66_10/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.TTbarJets.root")' > & eventcounts.reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.TTbarJets.root.log &
