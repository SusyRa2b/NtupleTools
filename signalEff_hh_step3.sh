#!/bin/sh

root -b -l -q signalEff_writetxt_compile.C

root -b -l -q "signalEff_writetxt_run.C(\"counts\",\"TChiHH\")"
root -b -l -q "signalEff_writetxt_run.C(\"JES\",\"TChiHH\")"
root -b -l -q "signalEff_writetxt_run.C(\"JER\",\"TChiHH\")"
root -b -l -q "signalEff_writetxt_run.C(\"ISR\",\"TChiHH\")"
#pdf uncertainties
###these need some cleaning up
#root -b -l -q "signalEff2012_PDF_details.C++(\"CTEQMSTW\",\"TChiHH\")" 
#root -b -l -q "signalEff2012_PDF_details.C++(\"CTEQNNPDF\",\"TChiHH\")" 
#root -b -l -q "signalEff2012_PDF_details_writetxt.C++(\"TChiHH\")"





