#!/bin/sh

root -b -l -q signalEff_writetxt_compile.C

root -b -l -q "signalEff_writetxt_run.C(\"counts\",\"TChiHH\")"
root -b -l -q "signalEff_writetxt_run.C(\"JES\",\"TChiHH\")"
root -b -l -q "signalEff_writetxt_run.C(\"JER\",\"TChiHH\")"
root -b -l -q "signalEff_writetxt_run.C(\"ISR\",\"TChiHH\")"
root -b -l -q "signalEff_writetxt_run.C(\"PU\",\"TChiHH\")"

#pdf uncertainties
##for hh search, the 3 is just a dummy argument. doesn't mean anything but needs to be there.
root -b -l -q "signalEff2012_PDF_details.C++(\"CTEQMSTW\",\"TChiHH\",3,true)" 
root -b -l -q "signalEff2012_PDF_details.C++(\"CTEQNNPDF\",\"TChiHH\",3,true)" 
root -b -l -q "signalEff2012_PDF_details_writetxt.C++(\"TChiHH\")"
