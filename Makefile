# $Id: Makefile,v 1.3 2010/12/08 11:49:30 rahatlou Exp $
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

SRCDIR=src/
BINFILES =  $(SRCDIR)analysis.cc $(SRCDIR)redntpApp.cc

PROGRAMS = $(patsubst $(SRCDIR)%.cc,%,$(BINFILES))


# --- External configuration ---------------------------------
CC         = g++ -Wno-deprecated
CCFLAGS    =  -g 
MFLAGS     = -MM
INCLUDES   =
WORKDIR    = tmp/
LIBDIR     = $(WORKDIR)
OBJDIR=$(WORKDIR)/objects/
# -------------------------------------------------------------

ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')
ROOFIT_LIBS := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIB= | sed 's|LIB=||')
ROOFIT_LIBS += $(shell cd $(CMSSW_BASE); scram tool info roofit | grep LIB= | sed 's|LIB=||') 


INCLUDES += -I.  -I.. -I$(ROOTSYS)/include  -I$(ROOFIT_INCLUDE)/ -I$(SRCDIR) 
ROOTSYS  ?= ERROR_RootSysIsNotDefined

EXTRALIBS  :=  -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit -lRIO -lTMVA

# CC files excluding the binaries
CCFILES=$(filter-out $(BINFILES),$(wildcard $(SRCDIR)*.cc))

# List of all object files to build
OOLIST=$(patsubst %.cc,%.o,$(CCFILES))
OLIST=$(patsubst $(SRCDIR)%.o,$(OBJDIR)/%.o,$(OOLIST))

# Implicit rule to compile all classes
$(OBJDIR)/%.o : $(SRCDIR)%.cc
	@echo "Compiling $<"
	@mkdir -p $(OBJDIR)
	@$(CC) $(CCFLAGS) -c $< -o $@ $(INCLUDES)


$(PROGRAMS) : $(OLIST)
	@echo "Linking $@"
	@$(CC) $(CCFLAGS)  $(INCLUDES) $(OLIST) \
	$(ROOTLIBS) $(EXTRALIBS) -o $(WORKDIR)/$@   $(patsubst %,$(SRCDIR)%.cc,$@)

default : analysis 

all : ${PROGRAMS}

clean:
	rm -Rf $(WORKDIR)/*
	@#rm -f $(LIBFILE)
	@rm -Rf *.o

veryclean:
	rm -Rf $(WORKDIR)

#.PHONY : shlib lib default clean
