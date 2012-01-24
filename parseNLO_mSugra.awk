##goal: turn sanjay's files from here:
#http://uaf-2.t2.ucsd.edu/~spadhi/slha/xsection/2012/m0_m12_40_m500_1/
#into the C++ friendly format
#Harold had a tool parseNLO.py but it seems to barf on the latest output from Sanjay for some reason
#(which is weird, because the format is almost unchanged)

#usage:
#awk -F'|' -f parseNLO_mSugra.awk msugra_m0_m12_40_m500_1_NLO_1.0.txt > NLOxsec_tanb40_10.txt

#just to throw out the first line
/scale/ {
  m0str="";
  m12str="";
#first look for m0=1234,
  match_result =  match($2,"m0=[0-9][0-9][0-9]?[0-9]?,");
  if (match_result == 0) {
    print "Problem! No match found!";
    abort();
  }
  else {
    m0str=substr($2, match_result+3,RLENGTH-4);
  }

#first look for m1/2=1234,
  match_result =  match($2,"m1/2=[0-9][0-9][0-9]?[0-9]?,");
  if (match_result == 0) {
    print "Problem! No match found!";
    abort();
  }
  else {
    m12str=substr($2, match_result+5,RLENGTH-6);
  }

  print m0str,m12str,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12;

}
