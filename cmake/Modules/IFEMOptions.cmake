OPTION(IFEM_USE_OPENMP         "Compile with OpenMP support?"         ON)
OPTION(IFEM_USE_LRSPLINES      "Compile with LR-splines support?"     ON)
OPTION(IFEM_USE_SUPERLU        "Compile with SuperLU support?"        ON)
OPTION(IFEM_USE_SUPERLU_MT     "Compile with SuperLU_MT support?"     OFF)
OPTION(IFEM_USE_PETSC          "Compile with PETSc support?"          ON)
OPTION(IFEM_USE_MPI            "Compile with MPI support?"            OFF)
OPTION(IFEM_USE_ISTL           "Compile with dune-istl support?"      ON)
OPTION(IFEM_USE_SLEPC          "Compile with SLEPc support?"          OFF)
OPTION(IFEM_USE_SPR            "Compile with SPR support?"            OFF)
OPTION(IFEM_USE_SAMG           "Compile with SAMG support?"           OFF)
OPTION(IFEM_USE_HDF5           "Compile with HDF5 support?"           ON)
OPTION(IFEM_USE_VTFWRITER      "Compile with VTFWriter support?"      ON)
OPTION(IFEM_USE_UMFPACK        "Compile with UMFPACK support?"        ON)
OPTION(IFEM_USE_CEREAL         "Compile with cereal support?"         ON)
OPTION(IFEM_USE_ZOLTAN         "Compile with zoltan support?"         OFF)
OPTION(IFEM_AS_SUBMODULE       "Compile IFEM as a submodule of apps?" OFF)
OPTION(IFEM_WHOLE_PROG_OPTIM   "Compile IFEM with link-time optimizations?" OFF)
OPTION(IFEM_TEST_MEMCHECK      "Run tests through valgrind?"          OFF)
OPTION(IFEM_BUILD_TESTING      "Build testing by default?"            OFF)
OPTION(IFEM_SERIAL_TESTS_IN_PARALLEL "Run serial tests in parallel builds?" ON)
OPTION(IFEM_INSTALL_DOXY       "Install documentation?" ON)
OPTION(IFEM_USE_NATIVE         "Enable native tuning in release builds?" ON)
