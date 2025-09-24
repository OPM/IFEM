# -*- mode: cmake -*-
#
#  This file is part of the Feel library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2010-01-22
#
#  Copyright (C) 2010 Université Joseph Fourier
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 3.0 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
set(OLD_PKG $ENV{PKG_CONFIG_PATH})
set(ENV{PKG_CONFIG_PATH} "$ENV{SLEPC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig:$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig")
set(OLD_ALLOW $ENV{PKG_CONFIG_ALLOW_SYSTEM_CFLAGS})
set(ENV{PKG_CONFIG_ALLOW_SYSTEM_CFLAGS} 1)
pkg_check_modules(SLEPC SLEPc>=3.6.3)
set(ENV{PKG_CONFIG_PATH} ${OLD_PKG})
set(ENV{PKG_CONFIG_ALLOW_SYSTEM_CFLAGS} ${OLD_ALLOW})

set(SLEPC_LIBRARIES ${SLEPC_STATIC_LDFLAGS})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SLEPc DEFAULT_MSG
                                  SLEPC_INCLUDE_DIRS SLEPC_LIBRARIES)

if(SLEPc_FOUND)
  add_library(SLEPc::SLEPc INTERFACE IMPORTED)
  set(SLEPC_LIBRARY_DIRS ${SLEPC_LIBRARIES})
  list(FILTER SLEPC_LIBRARY_DIRS INCLUDE REGEX -L.*)
  list(TRANSFORM SLEPC_LIBRARY_DIRS REPLACE "-L" "")
  list(FILTER SLEPC_LIBRARIES EXCLUDE REGEX -L.*)
  target_link_libraries(SLEPc::SLEPc INTERFACE ${SLEPC_LIBRARIES})
  target_include_directories(SLEPc::SLEPc INTERFACE ${SLEPC_INCLUDE_DIRS})
  target_compile_definitions(SLEPc::SLEPc INTERFACE HAS_SLEPC=1)
  target_link_libraries(SLEPc::SLEPc INTERFACE PETSc::PETSc)
  target_link_directories(SLEPc::SLEPc INTERFACE ${SLEPC_LIBRARY_DIRS})
endif()
