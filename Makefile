#------------------------------------------------------------------------------
# Description: Makefile to build executable BasicLoopCU
# Created:     Fri Jun  3 16:38:31 2011 by mkanalyzer.py
#
#               available switches:
#
#                 debugflag  (e.g., debugflag=-ggdb [default])
#                 cppflags
#                 cxxflags
#                 optflag
#                 verbose    (e.g., verbose=1)
#                 withcern   (e.g., withcern=1  expects to find CERN_LIB)
# Author:      Ben Kreis
#------------------------------------------------------------------------------
ifndef ROOTSYS
$(error *** Please set up Root)
endif
withroot:=1
#------------------------------------------------------------------------------
ifndef program
# default program name
program := BasicLoopCU
endif

cppsrcs	:= $(wildcard *.cpp)

cxxsrcs	:= $(wildcard *.cxx)

cxxhead	:= $(patsubst %.cxx,%.h,$(cxxsrcs))

ccsrcs  := $(wildcard *.cc)

# filter out all main programs except the one to be built
# 1. search files for main(...) and write list of files to .main
$(shell grep -H "main[(].*[)]" $(cxxsrcs) $(cppsrcs) $(ccsrcs)|cut -f1 -d: > .main)
# 2. send list back to Makefile
main	:= $(shell cat .main)
# 3. remove subset of files (including the file to be built) from main
main	:= $(filter-out $(program).cc $(program).cpp $(program).cxx treestream.cc,$(main))
# 4. remove the set main from the set of all files in the directory
cxxsrcs	:= $(filter-out $(main),$(cxxsrcs))
cppsrcs	:= $(filter-out $(main),$(cppsrcs))
ccsrcs	:= $(filter-out $(main),$(ccsrcs))

cxxobjs	:= $(patsubst %.cxx,tmp/%.o,$(cxxsrcs))
cppobjs	:= $(patsubst %.cpp,tmp/%.o,$(cppsrcs))
ccobjs	:= $(patsubst %.cc,tmp/%.o,$(ccsrcs))

objects	:= $(ccobjs) $(cppobjs) $(cxxobjs) tmp/cxxDict.o
say     := $(shell echo "Program: $(program)" >& 2)
#------------------------------------------------------------------------------
ifdef GCC_DIR
GCC_BIN_PREFIX	:= $(GCC_DIR)/bin/
else
GCC_BIN_PREFIX	:=
endif
C++	    := $(GCC_BIN_PREFIX)g++
LDSHARED:= $(GCC_BIN_PREFIX)g++
C++VER	:= $(shell $(C++) --version)
COMP	:= $(word 1, $(C++VER))
CTYPE	:= $(word 2, $(C++VER))
CVER	:= $(word 3, $(C++VER))
say 	:= $(shell echo "$(COMP) $(CTYPE) $(CVER)" >& 2)
#------------------------------------------------------------------------------
ifdef verbose
	AT =
else
	AT = @
endif
#------------------------------------------------------------------------------
# Products to compile/link against
#------------------------------------------------------------------------------
ifdef withcern
	ifndef CERN_LIB
		ifdef CERN_DIR
			CERN_LIB:= $(CERN_DIR)/lib
		else
			say:=$(error CERN_LIB must point to CERN lib directory)
		endif
	endif
	cernlib	:= -L$(CERN_LIB) -lpacklib -lmathlib -lkernlib
endif

roofitsys = $(shell echo $(ROOFITSYS))
rootsys = $(shell echo $(ROOTSYS))

cmssw_base = $(CMSSW_BASE)

ifdef withroot
	rootcpp	:= $(shell root-config --cflags)
	rootlib	:= $(shell root-config --glibs) -L $(rootsys)/lib -lTreePlayer -lMinuit -L $(roofitsys)/lib -lRooFitCore 
endif
#------------------------------------------------------------------------------
# Switches/includes
# debug flag is on by default
#------------------------------------------------------------------------------
debugflag:=-ggdb

ifndef optflag
	optflag:=-O2
endif

CPPFLAGS:= -I. $(rootcpp) $(cppflags) -I$(roofitsys)include 
CXXFLAGS:= -c -pipe $(optflag) -fPIC -Wall $(cxxflags) $(debugflag)
LDFLAGS	:= $(ldflags) $(debugflag)
LIBS	:= $(libs) $(rootlib) $(cernlib) 
#------------------------------------------------------------------------------
# Rules
#------------------------------------------------------------------------------
bin:	$(program)

$(program)	: $(objects)
	@echo "---> Linking $@"
	$(AT)$(LDSHARED) $(LDFLAGS) $(objects) $(LIBS) -o $@
	@echo ""

$(cxxobjs)	: tmp/%.o : %.cxx
	@echo "---> Compiling $<" 
	$(AT)$(CXX)	$(CXXFLAGS) $(CPPFLAGS) $< -o $@

$(cppobjs)	: tmp/%.o : %.cpp
	@echo "---> Compiling $<" 
	$(AT)$(CXX)	$(CXXFLAGS) $(CPPFLAGS) $< -o $@

$(ccobjs)	: tmp/%.o : %.cc
	@echo "---> Compiling $<" 
	$(AT)$(CXX)	$(CXXFLAGS) $(CPPFLAGS) $< -o $@

tmp/cxxDict.o   : tmp/cxxDict.cc
	@echo "---> Compiling $<" 
	$(AT)$(CXX)	$(CXXFLAGS) $(CPPFLAGS) $< -o $@

tmp/cxxDict.cc  : $(cxxhead)
	@echo "---> Creating $@" 
	$(AT)rootcint -f $@ -c $(CXXFLAGS) $(CPPFLAGS) -p $^

# Define clean up rule
clean   	:
	rm -rf tmp/*.o tmp/cxxDict.* $(program)
