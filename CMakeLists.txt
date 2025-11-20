cmake_minimum_required(VERSION 3.22)

set(IFEM_VERSION_MAJOR 0)
set(IFEM_VERSION_MINOR 9)
set(IFEM_VERSION_PATCH 0)
set(IFEM_VERSION "${IFEM_VERSION_MAJOR}.${IFEM_VERSION_MINOR}.${IFEM_VERSION_PATCH}")
set(IFEM_ABI_VERSION ${IFEM_VERSION_MAJOR}.${IFEM_VERSION_MINOR})

project(IFEM
        VERSION ${IFEM_VERSION}
        LANGUAGES C CXX Fortran
)

set(IFEM_PATH ${PROJECT_SOURCE_DIR})

# Generate header with version info
configure_file(src/IFEM.h.in IFEM.h @ONLY
)

# Add local modules
list(APPEND CMAKE_MODULE_PATH
  ${PROJECT_SOURCE_DIR}/cmake
  ${PROJECT_SOURCE_DIR}/cmake/Modules
  ${PROJECT_SOURCE_DIR}/cmake/Scripts
)

include(IFEMOptions)
include(IFEMBuildTypes)
include(IFEMLTO)

add_library(IFEM)

target_compile_features(IFEM PUBLIC cxx_std_17)
set_target_properties(
  IFEM
  PROPERTIES
    VERSION
      ${IFEM_VERSION}
    SOVERSION
      ${IFEM_ABI_VERSION}
)
ifem_lto(TARGET IFEM)

target_compile_definitions(IFEM PUBLIC $<$<CONFIG:Debug>:INDEX_CHECK=2>)
if(VERBOSE_DEBUG GREATER 0)
  target_compile_definitions(IFEM PRIVATE $<$<CONFIG:Debug>:SP_DEBUG=${VERBOSE_DEBUG}>)
  target_compile_definitions(IFEM INTERFACE $<$<CONFIG:Debug>:INT_DEBUG=${VERBOSE_DEBUG}>)
endif()

include(IFEMAddTargets)
include(IFEMCompilerFlags)
include(IFEMDoxy)
include(IFEMPortability)
include(IFEMSetupDependencies)
include(IFEMSetupVariables)

target_sources(
  IFEM
  PRIVATE
    3rdparty/expreval/except.cpp
    3rdparty/expreval/expr.cpp
    3rdparty/expreval/func.cpp
    3rdparty/expreval/funclist.cpp
    3rdparty/expreval/node.cpp
    3rdparty/expreval/parser.cpp
    3rdparty/expreval/vallist.cpp
    src/IFEM.C
    src/ASM/AlgEqSystem.C
    src/ASM/ASM1D.C
    src/ASM/ASM2D.C
    src/ASM/ASM3D.C
    src/ASM/ASMbase.C
    src/ASM/ASMLagBase.C
    src/ASM/ASMmxBase.C
    src/ASM/ASMs1D.C
    src/ASM/ASMs1DC1.C
    src/ASM/ASMs1DLag.C
    src/ASM/ASMs1DSpec.C
    src/ASM/ASMs2D.C
    src/ASM/ASMs2DC1.C
    src/ASM/ASMs2DIB.C
    src/ASM/ASMs2DLag.C
    src/ASM/ASMs2Dmx.C
    src/ASM/ASMs2DmxLag.C
    src/ASM/ASMs2Drecovery.C
    src/ASM/ASMs2DSpec.C
    src/ASM/ASMs2DTri.C
    src/ASM/ASMs3D.C
    src/ASM/ASMs3DLag.C
    src/ASM/ASMs3Dmx.C
    src/ASM/ASMs3DmxLag.C
    src/ASM/ASMs3Drecovery.C
    src/ASM/ASMs3DSpec.C
    src/ASM/ASMstruct.C
    src/ASM/ASMsupel.C
    src/ASM/ASMu1DLag.C
    src/ASM/ASMu2DLag.C
    src/ASM/ASMutils.C
    src/ASM/BasisFunctionCache.C
    src/ASM/BDFMats.C
    src/ASM/BlockElmMats.C
    src/ASM/DomainDecomposition.C
    src/ASM/dslbln.f
    src/ASM/dsolv3.f
    src/ASM/DualField.C
    src/ASM/ElmMats.C
    src/ASM/Field.C
    src/ASM/Fields.C
    src/ASM/FiniteElement.C
    src/ASM/GlbForceVec.C
    src/ASM/GlbL2projector.C
    src/ASM/GlbNorm.C
    src/ASM/HasGravityBase.C
    src/ASM/HHTMats.C
    src/ASM/IBGeometries.C
    src/ASM/ImmersedBoundaries.C
    src/ASM/IntegrandBase.C
    src/ASM/Lagrange.C
    src/ASM/LagrangeField2D.C
    src/ASM/LagrangeField3D.C
    src/ASM/LagrangeFields2D.C
    src/ASM/LagrangeFields3D.C
    src/ASM/NewmarkMats.C
    src/ASM/ReactionsOnly.C
    src/ASM/SAMpatch.C
    src/ASM/SplineField2D.C
    src/ASM/SplineField3D.C
    src/ASM/SplineField.C
    src/ASM/SplineFields1D.C
    src/ASM/SplineFields2D.C
    src/ASM/SplineFields2Dmx.C
    src/ASM/SplineFields3D.C
    src/ASM/SplineFields3Dmx.C
    src/ASM/SplineInterpolator.C
    src/Eig/eig_drv1.f
    src/Eig/eig_drv2.f
    src/Eig/eig_drv3.f
    src/Eig/eig_drv4.f
    src/Eig/eig_drv5.f
    src/Eig/eig_drv6.f
    src/Eig/EigSolver.C
    src/LinAlg/DenseMatrix.C
    src/LinAlg/DiagMatrix.C
    src/LinAlg/LinAlgInit.C
    src/LinAlg/LinSolParams.C
    src/LinAlg/MatVec.C
    src/LinAlg/ProcessAdm.C
    src/LinAlg/SAM.C
    src/LinAlg/SparseMatrix.C
    src/LinAlg/SPRMatrix.C
    src/LinAlg/SystemMatrix.C
    src/SIM/AdaptiveSetup.C
    src/SIM/AdaptiveSIM.C
    src/SIM/EigenModeSIM.C
    src/SIM/ForceIntegrator.C
    src/SIM/FunctionSum.C
    src/SIM/GenAlphaSIM.C
    src/SIM/HHTSIM.C
    src/SIM/ModelGenerator.C
    src/SIM/MultiStepSIM.C
    src/SIM/NewmarkNLSIM.C
    src/SIM/NewmarkSIM.C
    src/SIM/NodeVecFunc.C
    src/SIM/NonLinSIM.C
    src/SIM/SIM1D.C
    src/SIM/SIM2D.C
    src/SIM/SIM3D.C
    src/SIM/SIMadmin.C
    src/SIM/SIMargsBase.C
    src/SIM/SIMbase.C
    src/SIM/SIMdependency.C
    src/SIM/SIMgeneric.C
    src/SIM/SIMinput.C
    src/SIM/SIMmodal.C
    src/SIM/SIMmultiCpl.C
    src/SIM/SIMoptions.C
    src/SIM/SIMoutput.C
    src/SIM/SIMsolution.C
    src/SIM/SIMsupel.C
    src/SIM/TimeStep.C
    src/SIM/XMLInputBase.C
    src/Utility/AnaSol.C
    src/Utility/BDF.C
    src/Utility/Chebyshev.C
    src/Utility/ControlFIFO.C
    src/Utility/CoordinateMapping.C
    src/Utility/DataExporter.C
    src/Utility/ElementBlock.C
    src/Utility/ExprFunctions.C
    src/Utility/FieldFunctions.C
    src/Utility/ftnfiles.f
    src/Utility/Functions.C
    src/Utility/GaussQuadrature.C
    src/Utility/HDF5Base.C
    src/Utility/HDF5Reader.C
    src/Utility/HDF5Restart.C
    src/Utility/HDF5Writer.C
    src/Utility/IFEM_math.C
    src/Utility/LagrangeInterpolator.C
    src/Utility/Legendre.C
    src/Utility/LogStream.C
    src/Utility/MPC.C
    src/Utility/Profiler.C
    src/Utility/ScopedLogger.C
    src/Utility/SplineUtils.C
    src/Utility/StbImage.C
    src/Utility/StringUtils.C
    src/Utility/Tensor4.C
    src/Utility/Tensor.C
    src/Utility/TensorFunction.C
    src/Utility/TextureProperties.C
    src/Utility/ThreadGroups.C
    src/Utility/TractionField.C
    src/Utility/TriangleQuadrature.C
    src/Utility/Utilities.C
    src/Utility/Vec3Oper.C
    src/Utility/VTF.C
)

if(ISTL_FOUND)
  target_sources(
    IFEM
    PRIVATE
      src/LinAlg/ISTLMatrix.C
      src/LinAlg/ISTLSolParams.C
  )
endif()

if(PETSc_FOUND)
  target_sources(
    IFEM
    PRIVATE
      src/ASM/SAMpatchPETSc.C
      src/LinAlg/PETScMatrix.C
      src/LinAlg/PETScPCPerm.C
      src/LinAlg/PETScSchurPC.C
      src/LinAlg/PETScSolParams.C
  )
endif()

# TODO: Remove list when we can drop support for cmake 3.22
set(IFEM_PUBLIC_HEADERS
  ${PROJECT_BINARY_DIR}/IFEM.h
  3rdparty/stb_image.h
  src/ASM/AlgEqSystem.h
  src/ASM/ASMbase.h
  src/ASM/ASMenums.h
  src/ASM/ASMmxBase.h
  src/ASM/ASMs2Dmx.h
  src/ASM/ASMs3Dmx.h
  src/ASM/BasisFunctionCache.h
  src/ASM/BasisFunctionVals.h
  src/ASM/ASMstruct.h
  src/ASM/ASMunstruct.h
  src/ASM/BDFMats.h
  src/ASM/BlockElmMats.h
  src/ASM/ElmMats.h
  src/ASM/HHTMats.h
  src/ASM/NewmarkMats.h
  src/ASM/ElmNorm.h
  src/ASM/Field.h
  src/ASM/Fields.h
  src/ASM/GlbForceVec.h
  src/ASM/FiniteElement.h
  src/ASM/GlbNorm.h
  src/ASM/GlbL2projector.h
  src/ASM/GlobalIntegral.h
  src/ASM/IntegrandBase.h
  src/ASM/HasGravityBase.h
  src/ASM/ImmersedBoundaries.h
  src/ASM/Interface.h
  src/ASM/Integrand.h
  src/ASM/Lagrange.h
  src/ASM/LocalIntegral.h
  src/ASM/SAMpatch.h
  src/ASM/SplineField2D.h
  src/ASM/SplineField3D.h
  src/ASM/SplineFields1D.h
  src/ASM/SplineFields2D.h
  src/ASM/SplineFields3D.h
  src/ASM/TimeDomain.h
  src/ASM/ASMs1D.h
  src/ASM/ASMs2D.h
  src/ASM/ASMs3D.h
  src/ASM/ASM1D.h
  src/ASM/ASM2D.h
  src/ASM/ASM3D.h
  src/ASM/ASMs1DLag.h
  src/ASM/ASMs2DLag.h
  src/ASM/ASMs3DLag.h
  src/ASM/ASMu1DLag.h
  src/ASM/ASMu2DLag.h
  src/ASM/ASMLagBase.h
  src/ASM/ASMutils.h
  src/ASM/DomainDecomposition.h
  src/ASM/ItgPoint.h
  src/ASM/ReactionsOnly.h
  src/LinAlg/BLAS.h
  src/LinAlg/DenseMatrix.h
  src/LinAlg/DiagMatrix.h
  src/LinAlg/ISTLMatrix.h
  src/LinAlg/ISTLSolParams.h
  src/LinAlg/ISTLSupport.h
  src/LinAlg/LAPack.h
  src/LinAlg/LinAlgenums.h
  src/LinAlg/LinAlgInit.h
  src/LinAlg/LinSolParams.h
  src/LinAlg/matrix.h
  src/LinAlg/matrixnd.h
  src/LinAlg/MatVec.h
  src/LinAlg/PETScMatrix.h
  src/LinAlg/PETScPCPerm.h
  src/LinAlg/PETScSchurPC.h
  src/LinAlg/PETScSolParams.h
  src/LinAlg/PETScSupport.h
  src/LinAlg/print_tol.h
  src/LinAlg/ProcessAdm.h
  src/LinAlg/readIntVec.h
  src/LinAlg/SAM.h
  src/LinAlg/SparseMatrix.h
  src/LinAlg/SPRMatrix.h
  src/LinAlg/SystemMatrix.h
  src/SIM/AdaptiveSetup.h
  src/SIM/AdaptiveSIM.h
  src/SIM/EigenModeSIM.h
  src/SIM/ForceIntegrator.h
  src/SIM/FunctionSum.h
  src/SIM/GenAlphaSIM.h
  src/SIM/HHTSIM.h
  src/SIM/ModelGenerator.h
  src/SIM/MultiStepSIM.h
  src/SIM/NewmarkNLSIM.h
  src/SIM/NewmarkSIM.h
  src/SIM/NodeVecFunc.h
  src/SIM/NonLinSIM.h
  src/SIM/Property.h
  src/SIM/SIM1D.h
  src/SIM/SIM2D.h
  src/SIM/SIM3D.h
  src/SIM/SIMadmin.h
  src/SIM/SIMargsBase.h
  src/SIM/SIMbase.h
  src/SIM/SIMdependency.h
  src/SIM/SIMdummy.h
  src/SIM/SIMenums.h
  src/SIM/SIMgeneric.h
  src/SIM/SIMinput.h
  src/SIM/SIMmodal.h
  src/SIM/SIMmultiCpl.h
  src/SIM/SIMoptions.h
  src/SIM/SIMoutput.h
  src/SIM/SIMsolution.h
  src/SIM/SIMsupel.h
  src/SIM/TimeStep.h
  src/SIM/TopologySet.h
  src/SIM/XMLInputBase.h
  src/Utility/AnaSol.h
  src/Utility/BDF.h
  src/Utility/Catch2Support.h
  src/Utility/Chebyshev.h
  src/Utility/ControlFIFO.h
  src/Utility/CoordinateMapping.h
  src/Utility/DataExporter.h
  src/Utility/ElementBlock.h
  src/Utility/ExprFunctions.h
  src/Utility/FieldFunctions.h
  src/Utility/Function.h
  src/Utility/Functions.h
  src/Utility/GaussQuadrature.h
  src/Utility/HDF5Base.h
  src/Utility/HDF5Reader.h
  src/Utility/HDF5Restart.h
  src/Utility/HDF5Writer.h
  src/Utility/IFEM_math.h
  src/Utility/LagrangeInterpolator.h
  src/Utility/Legendre.h
  src/Utility/LogStream.h
  src/Utility/MPC.h
  src/Utility/MPCLess.h
  src/Utility/Point.h
  src/Utility/Profiler.h
  src/Utility/ScopedLogger.h
  src/Utility/SplineUtils.h
  src/Utility/StbImage.h
  src/Utility/StringUtils.h
  src/Utility/Tensor4.h
  src/Utility/TensorFunction.h
  src/Utility/Tensor.h
  src/Utility/TextureProperties.h
  src/Utility/ThreadGroups.h
  src/Utility/TractionField.h
  src/Utility/TriangleQuadrature.h
  src/Utility/Utilities.h
  src/Utility/Vec3.h
  src/Utility/Vec3Oper.h
  src/Utility/VTF.h
)

if(CMAKE_VERSION GREATER_EQUAL 3.23)
  target_sources(
    IFEM
    PUBLIC
    FILE_SET
      HEADERS
    BASE_DIRS
      ${PROJECT_BINARY_DIR}
      3rdparty
      src/ASM
      src/LinAlg
      src/SIM
      src/Utility
    FILES
      ${IFEM_PUBLIC_HEADERS}
  )
endif()

target_include_directories(IFEM PUBLIC
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/ASM>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/Eig>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/LinAlg>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/SIM>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/Utility>
)

target_include_directories(IFEM PRIVATE
  3rdparty
  3rdparty/autodiff
  3rdparty/expreval
)

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.23)
  target_sources(
    IFEM
    PRIVATE
    FILE_SET
      private_headers
    TYPE
      HEADERS
    BASE_DIRS
      src/ASM
    FILES
      src/ASM/ASMs2DC1.h
      src/ASM/ASMs1DSpec.h
      src/ASM/ASMs2DC1.h
      src/ASM/ASMs2DIB.h
      src/ASM/ASMs2DmxLag.h
      src/ASM/ASMs2DSpec.h
      src/ASM/ASMs2DTri.h
      src/ASM/ASMs3DmxLag.h
      src/ASM/ASMs3DSpec.h
      src/ASM/ASMsupel.h
      src/ASM/DualField.h
      src/ASM/FieldBase.h
      src/ASM/IBGeometries.h
      src/ASM/LagrangeField2D.h
      src/ASM/LagrangeField3D.h
      src/ASM/LagrangeFields2D.h
      src/ASM/LagrangeFields3D.h
      src/ASM/SAMpatchPETSc.h
      src/ASM/SplineField.h
      src/ASM/SplineFields2Dmx.h
      src/ASM/SplineFields3Dmx.h
      src/ASM/SplineInterpolator.h
  )
endif()

if(LRSpline_FOUND)
  target_sources(
    IFEM
    PRIVATE
      src/ASM/LR/ASMLRSpline.C
      src/ASM/LR/ASMu2D.C
      src/ASM/LR/ASMu2DC1.C
      src/ASM/LR/ASMu2DIB.C
      src/ASM/LR/ASMu2Dmx.C
      src/ASM/LR/ASMu2Dnurbs.C
      src/ASM/LR/ASMu2Drecovery.C
      src/ASM/LR/ASMu3D.C
      src/ASM/LR/ASMu3Dmx.C
      src/ASM/LR/ASMu3Drecovery.C
      src/ASM/LR/GlobalNodes.C
      src/ASM/LR/LRSplineField2D.C
      src/ASM/LR/LRSplineField3D.C
      src/ASM/LR/LRSplineField.C
      src/ASM/LR/LRSplineFields2D.C
      src/ASM/LR/LRSplineFields2Dmx.C
      src/ASM/LR/LRSplineFields3D.C
      src/ASM/LR/LRSplineFields3Dmx.C
  )


  #TODO: Remove list when we can drop support for cmake 3.22
  set(IFEM_LR_HEADERS
    src/ASM/LR/ASMLRSpline.h
    src/ASM/LR/ASMu2DC1.h
    src/ASM/LR/ASMu2D.h
    src/ASM/LR/ASMu2DIB.h
    src/ASM/LR/ASMu2Dmx.h
    src/ASM/LR/ASMu3D.h
    src/ASM/LR/ASMu3Dmx.h
    src/ASM/LR/GlobalNodes.h
    src/ASM/LR/LRSplineField2D.h
    src/ASM/LR/LRSplineField3D.h
    src/ASM/LR/LRSplineField.h
    src/ASM/LR/LRSplineFields2D.h
    src/ASM/LR/LRSplineFields2Dmx.h
    src/ASM/LR/LRSplineFields3D.h
    src/ASM/LR/LRSplineFields3Dmx.h
  )

  if(CMAKE_VERSION GREATER_EQUAL 3.23)
    target_sources(
      IFEM
      PUBLIC
      FILE_SET
        HEADERS
      FILES
        ${IFEM_LR_HEADERS}
    )
  else()
    list(APPEND IFEM_PUBLIC_HEADERS ${IFEM_LR_HEADERS})
  endif()
endif()

target_compile_definitions(IFEM PUBLIC
  Real=double
  epsZ=1.0e-12
  PROFILE_LEVEL=3
)

include(IFEMDoxy)
ifem_add_doc_target(APP ifem DOX IFEM)

enable_testing()
include(IFEMTesting)

IFEM_add_test_app(
    NAME
      IFEM
    SOURCES
      src/ASM/Test/TestASMs1D.C
      src/ASM/Test/TestASMs2D.C
      src/ASM/Test/TestASMs2Dmx.C
      src/ASM/Test/TestASMs3D.C
      src/ASM/Test/TestASMs3Dmx.C
      src/ASM/Test/TestASMsupel.C
      src/ASM/Test/TestASMu1D.C
      src/ASM/Test/TestASMu2DLag.C
      src/ASM/Test/TestBlockElmMats.C
      src/ASM/Test/TestDomainDecomposition.C
      src/ASM/Test/TestFiniteElement.C
      src/ASM/Test/TestLagrangeField.C
      src/ASM/Test/TestLagrangeFields.C
      src/ASM/Test/TestMatlabPatch.C
      src/ASM/Test/TestPiolaMapping.C
      src/ASM/Test/TestSplineField.C
      src/ASM/Test/TestSplineFields.C
      src/Eig/Test/TestEigSolver.C
      src/LinAlg/Test/TestDiagMatrix.C
      src/LinAlg/Test/TestLinSolParams.C
      src/LinAlg/Test/TestMatrix.C
      src/LinAlg/Test/TestMatVec.C
      src/LinAlg/Test/TestSAM.C
      src/LinAlg/Test/TestSparseMatrix.C
      src/SIM/Test/TestInitialConditions.C
      src/SIM/Test/TestModelGenerator.C
      src/SIM/Test/TestNewmark.C
      src/SIM/Test/TestNonLinSIM.C
      src/SIM/Test/TestSIM.C
      src/SIM/Test/TestSolution.C
      src/SIM/Test/TestXMLInputBase.C
      src/Utility/Test/TestAnaSol.C
      src/Utility/Test/TestBDF.C
      src/Utility/Test/TestChebyshev.C
      src/Utility/Test/TestControlFIFO.C
      src/Utility/Test/TestCoordinateMapping.C
      src/Utility/Test/TestElementBlock.C
      src/Utility/Test/TestExprAutoDiff.C
      src/Utility/Test/TestFunctions.C
      src/Utility/Test/TestGaussQuadrature.C
      src/Utility/Test/TestImageReader.C
      src/Utility/Test/TestLegendre.C
      src/Utility/Test/TestLinearFunc.C
      src/Utility/Test/TestScopedLogger.C
      src/Utility/Test/TestSplineUtils.C
      src/Utility/Test/TestStringUtils.C
      src/Utility/Test/TestTensor4.C
      src/Utility/Test/TestTensor.C
      src/Utility/Test/TestThreadGroups.C
      src/Utility/Test/TestTractionField.C
      src/Utility/Test/TestTriangleQuadrature.C
      src/Utility/Test/TestUtilities.C
      src/Utility/Test/TestVec3Oper.C
    WORKDIR
      ${IFEM_PATH}
    LIBRARIES
      IFEM
)

if(TARGET IFEM-test)
  target_include_directories(
    IFEM-test
    PRIVATE
      ${PROJECT_BINARY_DIR}
      3rdparty
      3rdparty/expreval
      3rdparty/autodiff
      src/ASM
      src/Eig
      src/LinAlg
      src/SIM
      src/Utility
  )
  if(CMAKE_VERSION GREATER_EQUAL 3.23)
    target_sources(
      IFEM-test
      PRIVATE
      FILE_SET
        HEADERS
      FILES
        src/ASM/Test/ASM2DTests.h
        src/ASM/Test/ASMCube.h
        src/ASM/Test/ASMSquare.h
        src/ASM/Test/ASMTrap.h
        src/ASM/Test/ASMTrap3.h
        src/ASM/Test/FieldTests.h
        src/LinAlg/Test/MatrixTests.h
    )
  endif()

  if(PETSc_FOUND)
    target_sources(IFEM-test PRIVATE src/LinAlg/Test/TestPETScMatrix.C)
  endif()

  if(ISTL_FOUND)
    target_sources(IFEM-test PRIVATE src/LinAlg/Test/TestISTLMatrix.C)
  endif()

  if(SPR_FOUND)
    target_sources(IFEM-test PRIVATE src/LinAlg/Test/TestSPRMatrix.C)
  endif()

  if(LRSpline_FOUND)
    target_sources(
      IFEM-test
      PRIVATE
        src/ASM/LR/Test/TestASMu2D.C
        src/ASM/LR/Test/TestASMu2Dmx.C
        src/ASM/LR/Test/TestASMu3D.C
        src/ASM/LR/Test/TestASMu3Dmx.C
        src/ASM/LR/Test/TestGlobalNodes.C
        src/ASM/LR/Test/TestLRSplineField.C
        src/ASM/LR/Test/TestLRSplineFields.C
    )
    if(CMAKE_VERSION GREATER_EQUAL 3.23)
      target_sources(
        IFEM-test
        PRIVATE
        FILE_SET
          HEADERS
        FILES
          src/ASM/LR/Test/ASMuCube.h
          src/ASM/LR/Test/ASMuSquare.h
          src/ASM/LR/Test/ASMuTrap.h
          src/ASM/LR/Test/ASMuTrap3.h
      )
    endif()
  endif()

  if(HDF5_FOUND)
    target_sources(IFEM-test PRIVATE src/Utility/Test/TestFieldFunctions.C)
    if(LRSpline_FOUND)
      target_sources(IFEM-test PRIVATE src/Utility/Test/TestFieldFunctionsLR.C)
    endif()
  endif()
endif()


IFEM_add_test_app(
  NAME
    IFEM_noblas
  SOURCES
    src/LinAlg/Test/NoBlas/TestMatrixFallback.C
  WORKDIR
    ${IFEM_PATH}
  LIBRARIES
    IFEM
)

IFEM_add_test_app(
    NAME
      IFEM-MPI
    SOURCES
      src/ASM/Test/MPI/TestDomainDecomposition.C
    WORKDIR
      ${IFEM_PATH}
    PARALLEL
      4
    LIBRARIES
      IFEM
)

if(TARGET IFEM-MPI-test)
  if(PETSC_FOUND)
    target_sources(IFEM-MPI-test PRIVATE src/LinAlg/Test/MPI/TestPETScMatrix.C)
    if(ARPACK_FOUND AND ZOLTAN_FOUND)
      target_soruces(IFEM-MPI-test PRIVATE src/Eig/Test/MPI/TestEigSolver.C)
    endif()
  endif()
  if(ISTL_FOUND)
    target_sources(IFEM-MPI-test PRIVATE src/LinAlg/Test/MPI/TestISTLMatrix.C)
  endif()
endif()

add_subdirectory(Apps/Common)

if(IFEM_AS_SUBMODULE OR IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
else()
  ifem_add_check_target()
endif()

export(TARGETS IFEM IFEMAppCommon FILE IFEMTargets.cmake)

include(GNUInstallDirs)

target_include_directories(IFEM INTERFACE
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/IFEM>
)

include(CMakePackageConfigHelpers)
configure_package_config_file(cmake/Templates/IFEMConfig-install.cmake.in
                              ${PROJECT_BINARY_DIR}/IFEMConfig-install.cmake
                              INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/IFEM
)
configure_package_config_file(cmake/Templates/IFEMConfig.cmake.in
                              ${PROJECT_BINARY_DIR}/IFEMConfig.cmake
                              INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/IFEM
)

install(
  FILES
    ${PROJECT_BINARY_DIR}/IFEMConfig-install.cmake
  RENAME
    IFEMConfig.cmake
  DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake/IFEM
)

install(
  FILES
    cmake/IFEMAddTargets.cmake
    cmake/IFEMDoxy.cmake
    cmake/IFEMFindDependencies.cmake
    cmake/IFEMLTO.cmake
    cmake/IFEMOptionsDownstream.cmake
    cmake/IFEMTesting.cmake
    cmake/Modules/FindARPACK.cmake
    cmake/Modules/FindCBLAS.cmake
    cmake/Modules/FindISTL.cmake
    cmake/Modules/FindPETSc.cmake
    cmake/Modules/FindSAMG.cmake
    cmake/Modules/FindSLEPc.cmake
    cmake/Modules/FindSPR.cmake
    cmake/Modules/FindSuiteSparse.cmake
    cmake/Modules/FindSuperLU.cmake
    cmake/Modules/FindTinyXML2.cmake
    cmake/Modules/FindVTFWriter.cmake
  DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake/IFEM
)

if(CMAKE_VERSION GREATER_EQUAL 3.23)
  install(
      TARGETS
        IFEM
      EXPORT
        IFEM
      DESTINATION
        ${CMAKE_INSTALL_LIBDIR}
      FILE_SET
        HEADERS
      DESTINATION
        ${CMAKE_INSTALL_INCLUDEDIR}/IFEM
  )
else()
  install(
      TARGETS
        IFEM
      EXPORT
        IFEM
      DESTINATION
        ${CMAKE_INSTALL_LIBDIR}
  )
  install(
    FILES
      ${IFEM_PUBLIC_HEADERS}
    DESTINATION
      ${CMAKE_INSTALL_INCLUDEDIR}/IFEM
  )
endif()

install(
  EXPORT
    IFEM
  DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake/IFEM
  FILE
    IFEMTargets.cmake
)

install(FILES src/IFEM-test.C DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/IFEM/src)

# documentation and license
install(
  FILES
    doc/petsc.pdf
    doc/sim-input.pdf
    doc/Tutorials/GettingStarted.pdf
  DESTINATION
    ${CMAKE_INSTALL_DOCDIR}
)
