#!/bin/sh

root -b -l -q signalEff_writetxt_compile.C

#bool args are useISR,rawCounts
#for hh we want true,false
root -b -l -q "signalEff_writetxt_run.C(\"counts\",\"TChiHH\",true,false)"
root -b -l -q "signalEff_writetxt_run.C(\"JES\",\"TChiHH\",true,false)"
root -b -l -q "signalEff_writetxt_run.C(\"JER\",\"TChiHH\",true,false)"
root -b -l -q "signalEff_writetxt_run.C(\"ISR\",\"TChiHH\",true,false)"
#pdf uncertainties
###these need some cleaning up
#root -b -l -q "signalEff2012_PDF_details.C++(\"CTEQMSTW\",\"TChiHH\")" 
#root -b -l -q "signalEff2012_PDF_details.C++(\"CTEQNNPDF\",\"TChiHH\")" 
#root -b -l -q "signalEff2012_PDF_details_writetxt.C++(\"TChiHH\")"





