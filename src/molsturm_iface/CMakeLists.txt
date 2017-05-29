set(SRC
	hartree_fock.cc
	Parameters.cc
)

#
# Parse build type
#
set(MOLSTURM_IFACE_BUILD_TYPE "Release" CACHE STRING
"The build type to use for the molsturm python bindings. \
Valid are Debug and Release.")

set(module_BUILD_TYPE "DEBUG")
if (MOLSTURM_IFACE_BUILD_TYPE STREQUAL "Release")
	set(module_BUILD_TYPE "RELEASE")
elseif(MOLSTURM_IFACE_BUILD_TYPE  STREQUAL "Debug")
	set(module_BUILD_TYPE "DEBUG")
else()
	message(FATAL_ERROR "The value of MOLSTURM_IFACE_BUILD_TYPE needs to be either \"Debug\" or \"Release\".")
endif()

#
# Find swig and python
#
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

set(Python_ADDITIONAL_VERSIONS 3)
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

#
# Setup cpp flags and includes for swig
# (Note, that swig makes heavy use of old-style casts)
#
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${module_BUILD_TYPE}}")
enable_if_compiles(CMAKE_CXX_FLAGS -Wno-missing-field-initializers)
enable_if_compiles(CMAKE_CXX_FLAGS -Wno-pedantic)
enable_if_compiles(CMAKE_CXX_FLAGS -Wno-old-style-cast)
enable_if_compiles(CMAKE_CXX_FLAGS -Wno-sign-compare)
enable_if_compiles(CMAKE_CXX_FLAGS -Wno-maybe-uninitialized)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#
# molsturm_iface swig module
#
set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(molsturm_iface.i PROPERTIES CPLUSPLUS ON)
swig_add_module(molsturm_iface python molsturm_iface.i ${SRC})

swig_link_libraries(molsturm_iface
	${PYTHON_LIBRARIES}
	${molsturm_${module_BUILD_TYPE}_TARGET}
)
target_compile_definitions(${SWIG_MODULE_molsturm_iface_REAL_NAME} PRIVATE
        ${MOLSTURM_DEFINITIONS} ${MOLSTURM_DEFINITIONS_${molsturm_BUILD_TYPE}}
	${module_BUILD_TYPE}=${module_BUILD_TYPE}
)

