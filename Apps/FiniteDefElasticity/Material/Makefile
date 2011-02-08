# $Id: Makefile,v 1.1 2010-12-18 20:03:54 kmo Exp $

LIB	= ../libMaterial.a

SRC	= amat3d.f cons2d.f cons3d.f eigs3d.f ieps3d.f liel3d.f mdma3d.f \
	  monh3d.f pcst3d.f push3d.f stnh3d.f vfnh3d.f wlog3d.f rprin0.f
OBJS    = $(SRC:.f=.o)

FC	= gfortran
FFLAGS	= -O3
#FFLAGS	= -g

all: feninc $(LIB)

clean:
	-rm $(OBJS)

$(LIB): $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)

feninc:
	ln -s . include
	ln -s . feninc
