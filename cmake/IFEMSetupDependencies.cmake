# Setup dependencies and link these to the IFEM library

include(IFEMFindDependencies)

target_link_libraries(IFEM PUBLIC
  CBLAS::CBLAS
  LAPACK::LAPACK
  GoTools::GoToolsCore
  GoTools::GoTrivariate
  ARPACK::ARPACK
  tinyxml2::tinyxml2
  ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}
)

target_link_directories(IFEM PUBLIC ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})

if(IFEM_USE_SUPERLU OR IFEM_USE_SUPERLU_MT)
  if(IFEM_USE_SUPERLU_MT AND TARGET SuperLU::SuperLUMT)
    target_link_libraries(IFEM PUBLIC SuperLU::SuperLUMT)
    message(STATUS "Using SuperLU-MT")
  elseif(TARGET SuperLU::SuperLU)
    target_link_libraries(IFEM PUBLIC SuperLU::SuperLU)
    message(STATUS "Using SuperLU")
  endif()
endif()

if(IFEM_USE_LRSPLINES)
  if(TARGET LRSpline::LRSpline)
    target_link_libraries(IFEM PUBLIC LRSpline::LRSpline)
    target_compile_definitions(IFEM PUBLIC HAS_LRSPLINE=1)
  endif()
endif()

if(IFEM_USE_MPI)
  target_link_libraries(IFEM PUBLIC MPI::MPI_C MPI::MPI_CXX)
  target_compile_definitions(IFEM PUBLIC HAVE_MPI=1)
endif()

if(IFEM_USE_PETSC)
  if(TARGET PETSc::PETSc)
    target_link_libraries(IFEM PUBLIC PETSc::PETSc)
  endif()
endif()

if(IFEM_USE_SLEPC AND TARGET PETSc::PETSc)
  if(TARGET SLEPc::SLEPc)
    target_link_libraries(IFEM PUBLIC SLEPc::SLEPc)
  endif()
endif()

if(IFEM_USE_HDF5)
  if(TARGET hdf5::hdf5)
    target_link_libraries(IFEM PUBLIC hdf5::hdf5)
    target_compile_definitions(IFEM PUBLIC HAS_HDF5=1)
  endif()
endif()

if(NOT IFEM_USE_SPR MATCHES "OFF")
  if(TARGET SPR::SPR)
    target_link_libraries(IFEM PUBLIC SPR::SPR)
  endif()
endif()

if(IFEM_USE_SAMG)
  if(TARGET SAMG::SAMG)
    target_link_libraries(IFEM PUBLIC SAMG::SAMG)
  endif()
endif()

if(IFEM_USE_VTFWRITER)
  if(TARGET VTFWriter::VTFx)
    target_link_libraries(IFEM PUBLIC VTFWriter::VTFx)
  elseif(TARGET VTFWriter::VTF)
    target_link_libraries(IFEM PUBLIC VTFWriter::VTF)
  endif()
endif()

if(IFEM_USE_OPENMP)
  if(TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(IFEM PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(IFEM PUBLIC USE_OPENMP=1)
  endif()
endif()

if(IFEM_USE_UMFPACK)
  if(TARGET SuiteSparse::umfpack)
    target_link_libraries(IFEM PUBLIC SuiteSparse::umfpack)
    target_compile_definitions(IFEM PUBLIC HAS_UMFPACK=1)
  endif()
endif()

if(IFEM_USE_ISTL)
  if(TARGET Dune::ISTL)
    target_link_libraries(IFEM PUBLIC Dune::ISTL)
  endif()
endif()

if(IFEM_USE_CEREAL)
  if(TARGET cereal::cereal)
    target_link_libraries(IFEM PUBLIC cereal::cereal)
    target_compile_definitions(IFEM PUBLIC HAS_CEREAL=1)
    message(STATUS "Cereal serialization support enabled")
  endif()
endif()

if(IFEM_USE_ZOLTAN)
  if(TARGET trilinos_zoltan)
    target_include_directories(IFEM PUBLIC ${Zoltan_INCLUDE_DIRS})
    target_link_libraries(IFEM PUBLIC trilinos_zoltan)
    target_compile_definitions(IFEM PUBLIC HAS_ZOLTAN=1)
    message(STATUS "Zoltan support enabled")
  endif()
endif()

if(IFEM_USE_TRACY)
  if(TARGET Tracy::TracyClient)
    target_link_libraries(IFEM PUBLIC Tracy::TracyClient)
    target_compile_definitions(IFEM PUBLIC HAS_TRACY=1)
    message(STATUS "Tracy support enabled")
  endif()
endif()
