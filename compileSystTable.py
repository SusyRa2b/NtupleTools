#need to load the 5 search region files
#from each file pull the second column of each line

name1BL = 'signalSystDetails.LM9.ge1bLoose.dat'
f1BL = open(name1BL,'r')

name2BL = 'signalSystDetails.LM9.ge2bLoose.dat'
f2BL = open(name2BL,'r')

name1BT = 'signalSystDetails.LM9.ge1bTight.dat'
f1BT = open(name1BT,'r')

name2BT = 'signalSystDetails.LM9.ge2bTight.dat'
f2BT = open(name2BT,'r')

name3B = 'signalSystDetails.LM9.ge3bLoose.dat'
f3B = open(name3B,'r')

rows = []
values1BL = []
values2BL = []
values1BT = []
values2BT = []
values3B = []

#on the first try this worked fine except the 'number' field
#included a \n character. The rstrip seems to do the trick to get rid of it

for line in f1BL:
    mypair = line.split("&")
    label = mypair[0]
    number = mypair[1]
    rows.append(label)
    values1BL.append(number.rstrip())

for line in f2BL:
    mypair = line.split("&")
    number = mypair[1]
    values2BL.append(number.rstrip())

for line in f1BT:
    mypair = line.split("&")
    number = mypair[1]
    values1BT.append(number.rstrip())

for line in f2BT:
    mypair = line.split("&")
    number = mypair[1]
    values2BT.append(number.rstrip())

for line in f3B:
    mypair = line.split("&")
    number = mypair[1]
    values3B.append(number.rstrip())

print " & 1BL & 2BL & 1BT & 2BT & 3B \\\\"

while len(rows)!=0:
    print rows.pop(0),"&",values1BL.pop(0),"&",values2BL.pop(0),"&",values1BT.pop(0),"&",values2BT.pop(0),"&",values3B.pop(0),"\\\\"
