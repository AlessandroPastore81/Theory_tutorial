#SHELL = tcsh
# add .f and .o to the list of acceptable suffixes
.SUFFIXES: .f90 .f .o .tar .gz $(SUFFIXES)

VERSION = 1.00
#
# V0.01: Beginning of the project...
#

#variables can be defined for use with later rules
#  (eg turn on debugging)
TAR = tar cvf
GZIP = gzip --best
CAT = cat
#
BUILDER = gfortran
FFLAGS = -Wall -O3
OPT = 
LFLAGS = 

#
F90 = $(BUILDER)
LD  = $(BUILDER)
AR = ar
CP = cp
# RANLIB is useless on some systems (HPUX), replace it by ls or file.
RANLIB = ranlib

# Sources and utilities
SRC =   ws.f90 ws_param.f90 ws_mod.f90 ws_density.f90 ws_fields.f90 \
        boundary.f90

FILES = $(SRC) Makefile forces.param

# define your objects
OBJS = ws_param.o ws_mod.o ws_density.o ws_fields.o boundary.o

OBJMAIN = ws.o
MODS = ws_param.mod ws_auxiliary.mod densita.mod fields.mod rspace.mod

MAIN = ws
ARC = tot$(VERSION).tar
DIRS = out
TOT = ws_tot.f90
LIB = liblr.a

all: $(MAIN) $(DIRS)
lib: $(LIB)
tot: $(TOT)

# this rule defines how your executable is to be made
# and is the default target
$(MAIN): $(OBJS) $(OBJMAIN)
	$(LD) -o $(MAIN) $(OBJS) $(OBJMAIN) $(LFLAGS)

%.o: %.f90
	$(F90) $(FFLAGS) -o $@ -c $< $(OPT)

$(TOT): $(SRC) Makefile
	rm -f $(TOT)
	$(CAT) $(SRC) > $(TOT)

$(LIB): $(OBJS)
	$(AR) vr $(LIB) $(OBJS)
	$(RANLIB) $(LIB)

out:
	

# Dependences
 
 
 ws_param.o: ws_param.f90
 boundary.o: boundary.f90  ws_param.o
 ws_mod.o: ws_param.o ws_mod.f90
 ws_density.o: ws_density.f90 ws_param.o
 ws_fields.o: ws_fields.f90 ws_density.o  ws_param.o
 ws.o: ws.f90  ws_param.o ws_mod.o ws_density.o ws_fields.o \
      boundary.o

# house keeping

clean:
	rm -f *.a *.o *.mod *.d $(MAIN) $(MAIN2)

clean_all:
	rm -f *.a *.o *.mod $(MAIN) $(MAIN2) *~

