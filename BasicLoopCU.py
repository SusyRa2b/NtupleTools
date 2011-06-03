#!/usr/bin/env python
# -----------------------------------------------------------------------------
#  File:        BasicLoopCU.py
#  Description: Analyzer for ntuples created by TheNtupleMaker
#  Created:     Fri Jun  3 16:38:31 2011 by mkntanalyzer.py
#  Author:      Ben Kreis
# -----------------------------------------------------------------------------
from ROOT import *
from string import *
from BasicLoopCUlib import *
import os, sys, re
# -----------------------------------------------------------------------------
# -- Procedures and functions
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
def main():

	# Get number of events
	nevents = stream.size()
	print "Number of events:", nevents


	# Notes:
	#
	# 1. Use
	#   ofile = outputFile(cmdline.outputfile, stream)
	#
	# to skim events to output file in addition to writing out histograms.
	#
	# 2. Use
	#   ofile.addEvent(event-weight)
	#
	# to specify that the current event is to be added to the output file. If
	# omitted, the event-weight is taken to be 1.
	#
	# 3. Use
	#    ofile.count(cut-name, event-weight)
	#
	# to keep track, in the count histogram, of the number of events passing
	# a given cut. If omitted, the event-weight is taken to be 1. If you want
	# the counts in the count histogram to appear in a given order, specify
	# the order, before entering the event loop, as in the example below
	# 
	#   ofile.count("NoCuts", 0)
	#	ofile.count("GoodEvent", 0)
	#	ofile.count("Vertex", 0)
	#	ofile.count("MET", 0)
	
	ofile = outputFile(cmdline.outputfilename)

	# -------------------------------------------------------------------------
	# Define histograms
	# -------------------------------------------------------------------------
	setStyle()



	# -------------------------------------------------------------------------
	# Loop over events
	# -------------------------------------------------------------------------
	for entry in xrange(nevents):
		stream.read(entry)

		# -- Event selection

	stream.close()
	ofile.close()
# -----------------------------------------------------------------------------
main()
