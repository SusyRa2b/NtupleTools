#!/bin/sh

#final argument should be true for ISR samples and false for others
root -b -l -q "signalEff2012_writetxt_run.C(\"counts\",\"T5tttt\",3,true)"
root -b -l -q "signalEff2012_writetxt_run.C(\"JES\",\"T5tttt\",3,true)"
#needed for samples with ISR, not for others
root -b -l -q "signalEff2012_writetxt_run.C(\"ISR\",\"T5tttt\",3,true)"
root -b -l -q "signalEff2012_PDF_details.C++(\"CTEQMSTW\",\"T5tttt\",3,true)" 
root -b -l -q "signalEff2012_PDF_details.C++(\"CTEQNNPDF\",\"T5tttt\",3,true)" 
root -b -l -q "signalEff2012_PDF_details_writetxt.C++(\"T5tttt\")"

#this line is usually redundant, but not if we skip step 2
#root -b -l -q signalEff2012_writetxt_compile.C
#true/false is an ISR switch
#root -b -l -q "signalEff2012_writetxt_run.C(\"counts\",\"pMSSM4\",3,false)"
#root -b -l -q "signalEff2012_writetxt_run.C(\"JES\",\"pMSSM4\",3,false)"
#root -b -l -q "signalEff2012_PDF_details.C++(\"CTEQMSTW\",\"pMSSM4\",3,false)" 
#root -b -l -q "signalEff2012_PDF_details.C++(\"CTEQNNPDF\",\"pMSSM4\",3,false)" 
#root -b -l -q "signalEff2012_PDF_details_writetxt.C++(\"pMSSM4\")"

#not doing this again (needed to do it once to evaluate MET and JER systematics)
#root -b -l -q "signalEff2012_writetxt_run.C(\"MET\",\"T1bbbb\")"
#root -b -l -q "signalEff2012_writetxt_run.C(\"JER\",\"T1bbbb\")"




