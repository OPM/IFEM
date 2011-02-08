# $Id: Makefile,v 1.13 2011-02-04 17:31:40 kmo Exp $
#===============================================================================
#
# File: Makefile
#
# Created: Jan 12 2010
#
# Author: Knut Morten Okstad / SINTEF
#
# Description:
# Makefile for spline-based isogeometric finite element simulation applications.
#
#===============================================================================

# Comment in this line if PETSc and SLEPc is used
#include ${SLEPC_DIR}/conf/slepc_common_variables

SRCS0	= SIMLinEl2D.C SIMLinEl3D.C \
	  FiniteDefElasticity/NonlinearElasticity.C \
	  FiniteDefElasticity/NonlinearElasticityTL.C
SRCS1	= main_LinEl3D.C $(SRCS0)
SRCS2	= main_LinElStep.C $(SRCS0)
OBJS1	= $(SRCS1:.C=.o) FiniteDefElasticity/stiff_TL.o
OBJS2	= $(SRCS2:.C=.o) FiniteDefElasticity/stiff_TL.o
APP1	= linel3D
APP2	= linelStep

MYLIB	= ../libSplineFEM.a
GOTLIB	= -L../GoTools/lib -lGoTrivariate -lGoToolsCore
BLASLIB	= -larpack -llapack -lblas
SLULIB	= -lsuperlu
#SLULIB	= -L../SuperLU_MT/lib -lsuperlu -lpthread
#SAMLIB	= -L../SAM -lSPR -lASM -lEIG -lMAT -lMDR -lUTL
#SAMGLIB	= -L../SAMG/samg -lamg_coo
# Comment in this line if PETSc and SLEPc is used
#PETSCOPT = -DHAS_PETSC $(CC_INCLUDES)
# Comment in this line if parallel PETSc and SLEPc are used
#PETSCOPT = -DHAS_PETSC -DPARALLEL_PETSC $(CC_INCLUDES)
VTFLIB	= -L../GLviewExpressWriter/lib -lVTFExpressAPI
LIBS	= $(MYLIB) $(GOTLIB) $(VTFLIB) \
	  $(SAMLIB) $(SAMGLIB) $(SLULIB) $(SLEPC_LIB) $(BLASLIB)

FC      = gfortran
CXX	= g++

#CDEBUG	= -g -DINDEX_CHECK
CDEBUG	= -O2

FFLAGS  = -O3 -funroll-loops
CXXFLAGS = $(CDEBUG) -Dreal=double -DUSE_CBLAS \
	-I. -I../Integrands -I../SIM -I../ASM -I../LinAlg -I../Utility \
	-IFiniteDefElasticity $(PETSCOPT)

all: $(APP1) $(APP2)

$(APP1): $(OBJS1) $(MYLIB) Makefile
	$(CXX) -o $@ $(OBJS1) $(LIBS) -lgfortran

$(APP2): $(OBJS2) $(MYLIB) Makefile
	$(CXX) -o $@ $(OBJS2) $(LIBS) -lgfortran

../libSplineFEM.a:
	cd ..; make

clean:
	rm -f $(OBJS1) $(OBJS2) $(OBJS3)
	sed '1,/^# DO NOT DELETE/w Makefile.tmp' < Makefile > /dev/null
	mv Makefile.tmp Makefile

depend:
	sed '1,/^# DO NOT DELETE/w Makefile.tmp' < Makefile > /dev/null
	$(CXX) -M $(CXXFLAGS) *.C >> Makefile.tmp
	mv Makefile Makefile.bak
	mv Makefile.tmp Makefile

# DO NOT DELETE
