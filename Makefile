RELVERSION  = $(shell cat .release)

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

MakefileFullPath = $(abspath $(lastword $(MAKEFILE_LIST)))
MakefileDirFullPath = $(shell dirname $(MakefileFullPath))
INSTALLDIR = $(MakefileDirFullPath)/install.$(RELVERSION)/

CXX  = g++
CXX += -I./
CXX += -I./obj/
CXX += -I./src/

CXXFLAGS  = -g -Wall -fPIC -Wno-deprecated
CXXFLAGS += $(ROOTCFLAGS)
CXXFLAGS += $(ROOTLIBS)
CXXFLAGS += $(ROOTGLIBS)
CXXFLAGS += -std=c++14
CXXFLAGS += -fconcepts

OUTLIB = ./obj/

#----------------------------------------------------#

all: makedir principalfaces

makedir:
	mkdir -p $(OUTLIB);

.PHONY: printmakehelp_and_reminder
printmakehelp_and_reminder: principalfaces.cpp Makefile
	$(info  /******************************************************************************/)
	$(info  * task --> printmakehelp_and_reminder: principalfaces.cpp Makefile          *)
	$(info  * $$@ ----> $@                                         *)
	$(info  * $$< --------------------------------> $<                   *)
	$(info  * $$^ --------------------------------> $^          *)
	$(info  /******************************************************************************/)

principalfaces: principalfaces.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS)

clean:
	rm -f principalfaces
	rm -f *~
	rm -f $(OUTLIB)*.o
