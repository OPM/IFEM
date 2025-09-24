# Options for the IFEM library

option(IFEM_USE_OPENMP         "Compile with OpenMP support?"         ON)
option(IFEM_USE_LRSPLINES      "Compile with LR-splines support?"     ON)
option(IFEM_USE_SUPERLU        "Compile with SuperLU support?"        ON)
option(IFEM_USE_SUPERLU_MT     "Compile with SuperLU_MT support?"     OFF)
option(IFEM_USE_PETSC          "Compile with PETSc support?"          ON)
option(IFEM_USE_MPI            "Compile with MPI support?"            OFF)
option(IFEM_USE_ISTL           "Compile with dune-istl support?"      ON)
option(IFEM_USE_SLEPC          "Compile with SLEPc support?"          OFF)
option(IFEM_USE_SAMG           "Compile with SAMG support?"           OFF)
option(IFEM_USE_HDF5           "Compile with HDF5 support?"           ON)
option(IFEM_USE_VTFWRITER      "Compile with VTFWriter support?"      ON)
option(IFEM_USE_UMFPACK        "Compile with UMFPACK support?"        ON)
option(IFEM_USE_CEREAL         "Compile with cereal support?"         ON)
option(IFEM_USE_ZOLTAN         "Compile with zoltan support?"         OFF)
option(IFEM_USE_TRACY          "Enable tracer profiler?"              OFF)

set(IFEM_USE_SPR "OFF" CACHE STRING "Compile with SPR support?")
set_property(CACHE IFEM_USE_SPR PROPERTY STRINGS ON I32 I64 OFF)

include(IFEMOptionsDownstream)
