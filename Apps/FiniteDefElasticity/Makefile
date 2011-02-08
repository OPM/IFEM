# $Id: Makefile,v 1.7 2011-02-08 09:32:18 kmo Exp $
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

CSRC	= main_NonLinEl.C NonlinearElasticityTL.C NonlinearElasticityUL.C \
	  NonlinearElasticityULMX.C NonlinearElasticityULMixed.C \
	  NonlinearElasticity.C NeoHookeElasticity.C \
	  SIMFiniteDefEl.C ../SIMLinEl2D.C ../SIMLinEl3D.C
FSRC	= accKM.f accKMx.f stiff_TL.f
OBJS	= $(CSRC:.C=.o) $(FSRC:.f=.o)

APP	= nonlinel

CXX	= g++
FC	= gfortran

MYLIB	= ../../libSplineFEM.a libMaterial.a
GOTLIB	= -L../../GoTools/lib -lGoTrivariate -lGoToolsCore
VTFLIB	= -L../../GLviewExpressWriter/lib -lVTFExpressAPI
BLASLIB	= -larpack -llapack -lblas
SLULIB	= -lsuperlu
#SLULIB	= -L../../SuperLU_MT/lib -lsuperlu_mt_PTHREAD -lpthread
#SAMGLIB	= -L../../SAMG/samg -lamg_coo
LIBS	= $(MYLIB) $(GOTLIB) $(VTFLIB) $(SAMGLIB) $(SLULIB) $(BLASLIB)

#CDEBUG	= -g -DINDEX_CHECK=2 -DINT_DEBUG=2
CDEBUG	= -O2

FFLAGS  = -O3 -funroll-loops
CXXFLAGS = $(CDEBUG) -Dreal=double -DUSE_CBLAS -DUSE_FTNMAT \
	 -I. -I.. -I../../Integrands -I../../SIM -I../../ASM \
	 -I../../LinAlg -I../../Utility

all: $(APP)

$(APP): $(OBJS) $(MYLIB) Makefile
	$(CXX) -o $@ $(OBJS) $(LIBS) -lgfortran

../../libSplineFEM.a:
	cd ../..; make

libMaterial.a:
	cd Material; make

clean:
	rm -f $(OBJS)
	sed '1,/^# DO NOT DELETE/w Makefile.tmp' < Makefile > /dev/null
	mv Makefile.tmp Makefile

depend:
	sed '1,/^# DO NOT DELETE/w Makefile.tmp' < Makefile > /dev/null
	$(CXX) -M $(CXXFLAGS) *.C >> Makefile.tmp
	mv Makefile Makefile.bak
	mv Makefile.tmp Makefile

tar: clean
	tar cfv FDSolve.tar Makefile $(FSRC) $(CSRC) *.h \
	Material/Makefile Material/*.[fh]
	gzip FDSolve.tar

# DO NOT DELETE
