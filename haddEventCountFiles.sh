#!/bin/sh

#goal:
#take a list of eventcounts files and hadd similar ones

#example:

#want to hadd together:
#eventcounts.Isr0.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-MadGraph_T1bbbb*.root
#into
#eventcounts.Isr0.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.T1bbbbMG.root

sampleid='SMS-MadGraph_T5tttt'
sampleidout='T5tttt'
#sampleid='pMSSM12_MCMC1_120mh130_batch4_'
#sampleidout='pMSSM4'
finaldir='eventCounts-v68_3/ORIGINALS'

mkdir -p $finaldir

#this is for the samples with ISR (madgraph)
stub[0]='eventcounts.Isr0.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0'
stub[1]='eventcounts.mergebbins.Isr0.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0'
stub[2]='eventcounts.mergebbins.IsrDown.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0'
stub[3]='eventcounts.mergebbins.IsrUp.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0'
stub[4]='eventcounts.mergebbins.withpdfs.Isr0.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0'
stub[5]='eventcounts.pusyst.Isr0.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0'
stub[6]='eventcounts.mergebbins.Isr0.CSVM_PF2PATjets_JESup_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0'
stub[7]='eventcounts.mergebbins.Isr0.CSVM_PF2PATjets_JESdown_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0'

#this is for samples without ISR (pythia)
#stub[0]='eventcounts.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0'
#stub[1]='eventcounts.mergebbins.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0'
#stub[2]='eventcounts.mergebbins.withpdfs.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0'
#stub[3]='eventcounts.pusyst.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0'
#stub[4]='eventcounts.mergebbins.CSVM_PF2PATjets_JESup_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0'
#stub[5]='eventcounts.mergebbins.CSVM_PF2PATjets_JESdown_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0'

for t in "${stub[@]}"
do
echo $t
hadd ${t}.${sampleidout}.root ${t}.${sampleid}*.root
mv ${t}.${sampleid}*.root $finaldir
done

