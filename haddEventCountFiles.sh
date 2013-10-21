#!/bin/sh

#goal:
#take a list of eventcounts files and hadd similar ones

#example:

#want to hadd together:
#eventcounts.Isr0.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.SMS-MadGraph_T1bbbb*.root
#into
#eventcounts.Isr0.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.T1bbbbMG.root

sampleid='SMS-TChiHH'
sampleidout='TChiHH'
#sampleid='pMSSM12_MCMC1_120mh130_batch4_'
#sampleidout='pMSSM4'
finaldir='TChiHH'
inputdir='.'

mkdir -p $finaldir

#this is for the samples with ISR (madgraph)
stub[0]='eventcounts.Isr0.JES0_JER0_PFMETTypeI_METunc0_PUunc0_hpt20'
stub[1]='eventcounts.mergebbins.Isr0.JES0_JER0_PFMETTypeI_METunc0_PUunc0_hpt20'
stub[2]='eventcounts.mergebbins.IsrDown.JES0_JER0_PFMETTypeI_METunc0_PUunc0_hpt20'
stub[3]='eventcounts.mergebbins.IsrUp.JES0_JER0_PFMETTypeI_METunc0_PUunc0_hpt20'
stub[4]='eventcounts.mergebbins.withpdfs.Isr0.JES0_JER0_PFMETTypeI_METunc0_PUunc0_hpt20'
stub[5]='eventcounts.pusyst.Isr0.JES0_JER0_PFMETTypeI_METunc0_PUunc0_hpt20'
stub[6]='eventcounts.mergebbins.Isr0.JESup_JER0_PFMETTypeI_METunc0_PUunc0_hpt20'
stub[7]='eventcounts.mergebbins.Isr0.JESdown_JER0_PFMETTypeI_METunc0_PUunc0_hpt20'

for t in "${stub[@]}"
do
echo $t
hadd ${inputdir}/${t}.${sampleidout}.root ${inputdir}/${t}.${sampleid}*.root
mv ${inputdir}/${t}.${sampleid}*.root $finaldir
done

