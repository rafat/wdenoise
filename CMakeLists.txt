cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

set(PROJECT_NAME wdenoise)
project(${PROJECT_NAME} CXX C)

# src root path
set(WDENOISE_SRC_ROOT ${PROJECT_SOURCE_DIR} CACHE PATH "wdenoise source root")

# Global settings
set(GLOBAL_OUTPUT_PATH ${PROJECT_BINARY_DIR}/Bin)
set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR}/install)

# Sets global output directory for single configuration (GCC)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${GLOBAL_OUTPUT_PATH})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${GLOBAL_OUTPUT_PATH})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${GLOBAL_OUTPUT_PATH})

# Sets global output directory for sub-configurations (msvc, mingw)
foreach(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${GLOBAL_OUTPUT_PATH})
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${GLOBAL_OUTPUT_PATH})
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${GLOBAL_OUTPUT_PATH})
endforeach(OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES)

# make include globaly visible
set(PROJECT_WIDE_INCLUDE ${WDENOISE_SRC_ROOT}/header)
include_directories(${PROJECT_WIDE_INCLUDE})

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS -std=c++11)
else()
  add_definitions(-DWIN32) # TMP
endif()

option(BUILD_UT "Enable Unit test" ON)

find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND git submodule update --init --recursive wavelib
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/wavelib/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

add_subdirectory(wavelib)

# Global link directories
link_directories(${GLOBAL_OUTPUT_PATH})

# Copies files from build directory into install directory.

file(COPY ${WDENOISE_SRC_ROOT}/wavelib/header/wavelib.h DESTINATION ${WDENOISE_SRC_ROOT}/header)

# Remove executables and libraries that are not needed

file(GLOB_RECURSE MY_FILES
  ${GLOBAL_OUTPUT_PATH}/libwauxlib.a
  ${GLOBAL_OUTPUT_PATH}/wavelibLibTests.exe
  ${CMAKE_SOURCE_DIR}/test
)

#file(REMOVE ${MY_FILES})

if(BUILD_UT)
    include(CTest)
    enable_testing()
    add_subdirectory(unitTests)
endif()

add_subdirectory(src)
add_subdirectory(examples)
