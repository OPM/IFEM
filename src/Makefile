# $Id: Makefile,v 1.49 2011-02-08 13:10:06 rho Exp $
#===============================================================================
#
# File: Makefile
#
# Created: Jan 12 2010
#
# Author: Knut Morten Okstad / SINTEF
#
# Description:
# Makefile for the spline-based isogeometric finite element toolbox, IFEM.
#
#===============================================================================

# Comment in this line if PETSc and SLEPc are used
#include ${SLEPC_DIR}/conf/slepc_common_variables

FSRC	= Eig/eig_drv1.f Eig/eig_drv2.f Eig/eig_drv3.f \
	  Eig/eig_drv4.f Eig/eig_drv5.f Eig/eig_drv6.f
CSRC	= SIM/SIMinput.C SIM/SIMbase.C SIM/SIM3D.C SIM/SIM2D.C SIM/SIM1D.C \
	  SIM/NonLinSIM.C ASM/Lagrange.C \
	  ASM/AlgEqSystem.C ASM/ASMbase.C ASM/ASMmxBase.C ASM/ASMstruct.C \
	  ASM/ASMs3D.C ASM/ASMs3DLag.C ASM/ASMs3DSpec.C \
	  ASM/ASMs3Dmx.C ASM/ASMs3DmxLag.C \
	  ASM/ASMs2D.C ASM/ASMs2DLag.C ASM/ASMs2DSpec.C \
	  ASM/ASMs2Dmx.C ASM/ASMs2DmxLag.C \
	  ASM/ASMs1D.C ASM/ASMs1DLag.C ASM/ASMs1DSpec.C \
	  Integrands/AnalyticSolutions.C Integrands/AnalyticSolutionsStokes.C \
	  Integrands/Poisson.C \
	  Integrands/LinearElasticity.C Integrands/Elasticity.C Integrands/Stokes.C \
	  Integrands/StabilizedStokes.C Integrands/StabilizedNavierStokes.C \
	  Integrands/NavierStokesG2.C Integrands/NavierStokesG2CN.C \
	  Integrands/NavierStokesG2MP.C Integrands/NavierStokesG2GenTheta.C \
	  Integrands/ChorinVelPred.C Integrands/ChorinVelPredBE.C \
	  Integrands/ChorinVelPredBDF2.C \
	  Integrands/ChorinPressCorr.C Integrands/ChorinVelCorr.C \
	  LinAlg/SystemMatrix.C LinAlg/DenseMatrix.C LinAlg/SPRMatrix.C \
	  LinAlg/SparseMatrix.C LinAlg/PETScMatrix.C LinAlg/LinSolParams.C \
	  LinAlg/LinAlgInit.C LinAlg/MatVec.C LinAlg/SAM.C \
	  Eig/EigSolver.C ASM/SAMpatch.C ASM/SAMpatchPara.C ASM/GlbNorm.C \
	  Utility/CoordinateMapping.C Utility/MPC.C \
	  Utility/GaussQuadrature.C Utility/Legendre.C \
	  Utility/ElementBlock.C Utility/VTF.C \
	  Utility/Vec3Oper.C Utility/Tensor.C Utility/Function.C \
	  Utility/Utilities.C Utility/Profiler.C Utility/Functions.C \
	  Utility/AnaSol.C

OBJS	= $(CSRC:.C=.o) $(FSRC:.f=.o)
HDRS	= $(CSRC:.C=.h) SIM/SIMenums.h SIM/Property.h \
	  ASM/TimeDomain.h ASM/ElmMats.h ASM/ElmNorm.h \
	  ASM/GlobalIntegral.h ASM/LocalIntegral.h ASM/IntegrandBase.h \
	  Utility/MPCLess.h Utility/Vec3.h LinAlg/matrix.h

LIB	= libSplineFEM.a

CXX	= g++
FC	= gfortran

GOTOPT	= -IGoTools/include
VTFOPT	= -DHAS_VTFAPI -IGLviewExpressWriter/include
BLASOPT	= -DUSE_CBLAS
SLUOPT	= -DHAS_SUPERLU -I/usr/include/superlu
#SLUOPT	= -DHAS_SUPERLU_MT -ISuperLU_MT/SRC
#SAMOPT	= -DHAS_SPR -DUSE_F77SAM
#SAMGOPT	= -DHAS_SAMG -DSAMG_UNIX_LINUX -DSAMG_LCASE_USCORE -ISAMG
# Comment in this line if PETSc and SLEPc are used
#PETSCOPT = -DHAS_PETSC $(CC_INCLUDES)
# Comment in this line if parallel PETSc and SLEPc are used
#PETSCOPT = -DHAS_PETSC -DPARALLEL_PETSC $(CC_INCLUDES)

#CDEBUG	= -g -DINDEX_CHECK -DSP_DEBUG=3
CDEBUG	= -O2

#FFLAGS	= -g
FFLAGS	= -O3 -funroll-loops
CXXFLAGS = $(CDEBUG) -Dreal=double -DepsZ=1.0e-12 -DPROFILE_LEVEL=3 \
	   -ISIM -IASM -IEig -IIntegrands -ILinAlg -IUtility \
	   $(GOTOPT) $(SAMOPT) $(SAMGOPT) $(SLUOPT) $(PETSCOPT) $(BLASOPT) \
	   $(VTFOPT)

all: $(LIB)

libs:
	if [ -d ARPACK ]; cd ARPACK; make lib DIRS='$(PWD)/ARPACK/SRC $(PWD)/ARPACK/UTIL'
	if [ -d SAM ]; then cd SAM/ASM; make; cd ../SPR; make; fi
	if [ -d SuperLU ]; then cd SuperLU/SRC; make double; fi
	if [ -d SuperLU_MT ]; then cd SuperLU_MT/SRC; make double; fi
	if [ -d SAMG ]; then cd SAMG/samg; make; fi
	if [ -d GLviewExpressWriter ]; then \
	cd GLviewExpressWriter/source; make install; fi

$(LIB): $(OBJS) Makefile
	$(AR) $(ARFLAGS) $@ $(OBJS)

clean:
	rm -f $(OBJS)
	sed '1,/^# DO NOT DELETE/w Makefile.tmp' < Makefile > /dev/null
	mv Makefile.tmp Makefile

depend:
	sed '1,/^# DO NOT DELETE/w Makefile.tmp' < Makefile > /dev/null
	$(CXX) -M $(CXXFLAGS) $(CSRC) >> Makefile.tmp
	mv Makefile Makefile.bak
	mv Makefile.tmp Makefile
	for dir in ASM Eig Integrands LinAlg SIM Utility; do \
	sed -i "/^[A-Za-z].*: $$dir/s/.*:/$$dir\/&/" Makefile; done

tar: clean
	cd Apps; make clean
	tar cfv SplineFEM.tar README Makefile $(FSRC) $(CSRC) $(HDRS) \
	Apps/Makefile Apps/*.[Ch]
	gzip SplineFEM.tar

# DO NOT DELETE
