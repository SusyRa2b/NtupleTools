#! /bin/env python
# usage: python parseNLO.py
# This will generate a new file called NLOxsec_parsed.txt in a format
# easy to read for mini-files code.
#
# Generated file contains columns of
# m0 m1/2 ng ns nn ll sb ss tb bb gg sg

###This is a script from Harold
###I used in the Summer11 but it seems to barf on the new files for 2012.
###I can't figure out why, but I wrote a new awk script to replace it: parseNLO_mSugra.awk

f = file("NLOxsec.txt","r")
#f = file("xsec.txt","r")
new = open("NLOxsec_parsed.txt", "w")


for line in f.readlines():
  fields = line.split('|')
  point_fields = fields[1].split(",")
 
 
  if "(" in point_fields[0]:
     pfield = point_fields[0]
     bracket = pfield[pfield.find("(")+1:pfield.find(")")]
     point_fields[0] = pfield[pfield.find(")")+1:]
     point_fields.append(bracket)
 
  i = 0 
  for pfield in point_fields:
   if i < 2: 
    keyval = pfield.split("=")
    print keyval[1].strip().rstrip(),
    new.write(keyval[1].strip().rstrip() + " "),
    i = i + 1
  for nlo in fields[2:-1]:
    print nlo.strip().rstrip(),
    new.write(nlo.strip().rstrip() + " "),
  print "\n",
  new.write("\n"),
   # new.write(keyval[0])

f.close()
new.close()
