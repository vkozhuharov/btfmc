# $Id: GNUmakefile,v 1.1 2014/01/22 15:35:03 veni Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := exampleN02
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

ROOTCONFIG   := root-config
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)

include $(G4INSTALL)/config/architecture.gmk
CPPFLAGS  += $(ROOTCFLAGS) -g -fPIC
EXTRALIBS := $(ROOTLIBS)

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
