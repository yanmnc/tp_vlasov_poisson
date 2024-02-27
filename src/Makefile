#**************************************************************
#  Copyright Euratom-CEA
#  Authors : 
#     Virginie Grandgirard (virginie.grandgirard@cea.fr)
#
#  Code plasma beam : solving of the 2D Vlasov-Poisson
#    system for plasma beam study
#  
#  This software is governed by the CeCILL-B license 
#  under French law and abiding by the rules of distribution 
#  of free software.  You can  use, modify and redistribute 
#  the software under the terms of the CeCILL-B license as 
#  circulated by CEA, CNRS and INRIA at the following URL
#  "http://www.cecill.info". 
#**************************************************************

F90 = gfortran
BINDIR = ../build
EXEC = $(BINDIR)/VlasovPoisson.exe

SRCDIR = .
OBJDIR = ../build

OBJS = $(addprefix $(OBJDIR)/, prec_const.o globals.o read_data.o \
        utils.o geometry.o fft_NRF90.o LUsolving_tridiag.o \
        spline1d.o spline1d_natural_types.o spline1d_natural.o \
        spline1d_periodic_types.o spline1d_periodic.o interpolation.o \
        fdistribu2d.o efield.o advec1D_SL.o energy.o ascii_saving.o physics.o \
	RHS_operators.o resol_steps.o VlasovPoisson.o)

SOURCES = $(OBJS:.o=.f90) $(OBJS:.o=.F90)

.SUFFIXES: .f90 .F90 .o

debug: $(EXEC)

$(EXEC): $(OBJS) | $(BINDIR)
	$(F90) $(OBJS) $(LIBS) -o $@

$(BINDIR):
	mkdir -p $(BINDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.f90 | $(BINDIR)
	$(F90) -J$(OBJDIR) -c $(F90FLAGS) $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.F90 | $(BINDIR)
	$(F90) -J$(OBJDIR) -c $(F90FLAGS) $< -o $@

clean:
	@rm -f $(OBJDIR)/*
	@rmdir $(OBJDIR)
