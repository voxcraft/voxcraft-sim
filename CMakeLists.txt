cmake_minimum_required(VERSION 3.12.0)
# set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if (NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE)
endif()

set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)

project(voxcraft-sim LANGUAGES CXX CUDA)
set(CMAKE_CXX_STANDARD 17)

include_directories(src)
include_directories(src/VX3)
include_directories(src/Utils)
include_directories(src/VXA)
include_directories(src/old)

file( GLOB UTILS_SRC_cu src/Utils/*.cu )
file( GLOB VX3_SRC_cu src/VX3/*.cu )
file( GLOB VX3_SRC_cpp src/VX3/*.cpp )
file( GLOB VXA_SRC_cpp src/VXA/*.cpp )
file( GLOB OLD_SRC src/old/*.c* )
file(GLOB_RECURSE COLLISION_DETECT_SRC src/Cu-Collision-Detection/src/* src/Cu-Collision-Detection/include/*)

add_subdirectory(3rdparty/halloc)
include_directories(3rdparty/halloc/src/)

set(Thrust_DIR 3rdparty/thrust/thrust/cmake)
find_package(Thrust REQUIRED CONFIG)
thrust_create_target(Thrust)

add_library(voxelyze.3.0 ${UTILS_SRC_cu} ${VX3_SRC_cu} ${VX3_SRC_cpp} ${VXA_SRC_cpp} ${OLD_SRC} ${COLLISION_DETECT_SRC})
add_executable(voxcraft-sim src/Executables/voxcraft-sim.cpp)
add_executable(vx3_node_worker src/Executables/vx3_node_worker.cu)
add_executable(vx3_node_daemon src/Executables/vx3_node_daemon.cu)
add_executable(vx3_test src/Executables/vx3_test.cu)

target_link_libraries(vx3_node_worker PUBLIC voxelyze.3.0 halloc)
target_link_libraries(vx3_test PUBLIC voxelyze.3.0 halloc)
# for Boost
set(Boost_USE_STATIC_LIBS   ON)
find_package(Boost REQUIRED COMPONENTS filesystem thread chrono program_options)
target_link_libraries(voxelyze.3.0 PUBLIC ${Boost_LIBRARIES} Thrust)
target_link_libraries(voxcraft-sim PUBLIC ${Boost_LIBRARIES})
target_link_libraries(vx3_node_worker PUBLIC ${Boost_LIBRARIES})
target_link_libraries(vx3_node_daemon PUBLIC ${Boost_LIBRARIES})


# CUDA options
find_package(CUDA) # for CUDA_SELECT_NVCC_ARCH_FLAGS function
# if one want to force use certain Gpu Architecture, one can add `CUDA_ARCH` like `cmake -DCUDA_ARCH=arch=compute_61,code=sm_61 ..`
if (NOT DEFINED CUDA_ARCH)
  CUDA_SELECT_NVCC_ARCH_FLAGS(RET_ARCH_FLAGS)
  message("CUDA_SELECT_NVCC_ARCH_FLAGS returns: " ${RET_ARCH_FLAGS})
  list(LENGTH RET_ARCH_FLAGS FLAG_LENGTH)
  message("Flag length: " ${FLAG_LENGTH})
  if (FLAG_LENGTH GREATER 3)
    # CUDA_SELECT_NVCC_ARCH_FLAGS resulted more than 3 items, which is probably GPU not found error.
    message(FATAL_ERROR "Nvidia GPU not found. Abort.")
  endif()
  list(GET RET_ARCH_FLAGS -1 CUDA_ARCH) # use newest architecture
endif()
message("Using GPU Architecture: " ${CUDA_ARCH})
set(ARCH_FLAGS "-gencode" ${CUDA_ARCH})
list(APPEND CMAKE_CUDA_FLAGS -gencode=${CUDA_ARCH})

set_target_properties(voxelyze.3.0 PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(vx3_node_worker PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(vx3_test PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

string(TOUPPER "${CMAKE_BUILD_TYPE}" uppercase_CMAKE_BUILD_TYPE)
if(uppercase_CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  if (NOT DEFINED CUDA_DEBUG)
    set(CUDA_DEBUG "-G or not" ON)
  endif()
endif(uppercase_CMAKE_BUILD_TYPE STREQUAL "DEBUG")
#-G means debug info in GPU, so you can step in GPU kernels
if (CUDA_DEBUG)
  message("Set CUDA_DEBUG to ON.")
  set(CUDA_OPTIONS $<$<COMPILE_LANGUAGE:CUDA>:-G ${ARCH_FLAGS}>)
else (CUDA_DEBUG)
  message("Set CUDA_DEBUG to OFF.")
  set(CUDA_OPTIONS $<$<COMPILE_LANGUAGE:CUDA>: ${ARCH_FLAGS}>)
endif(CUDA_DEBUG)

target_compile_options(voxelyze.3.0     PUBLIC ${CUDA_OPTIONS})
target_compile_options(vx3_node_worker  PUBLIC ${CUDA_OPTIONS})
target_compile_options(vx3_test         PUBLIC ${CUDA_OPTIONS})

#for gtest
################################
# Testing
################################
option(test "Build all tests." OFF) # Makes boolean 'test' available.
if (test)
  if (APPLE)
    add_definitions(-DGTEST_USE_OWN_TR1_TUPLE)
    add_definitions(-D__GLIBCXX__)
  endif (APPLE)

  # This adds another subdirectory, which has 'project(gtest)'.
  add_subdirectory(test/googletest)

  enable_testing()

  # Include the gtest library. gtest_SOURCE_DIR is available due to
  # 'project(gtest)' above.
  include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

  ##############
  # Unit Tests
  ##############
  set(TEST_DIR test/testcase)
  file(GLOB TEST_FILES ${TEST_DIR}/test_*.c*)
  add_executable(runUnitTests ${TEST_FILES})

  # Standard linking to gtest stuff.
  target_link_libraries(runUnitTests gtest gtest_main)

  target_link_libraries(runUnitTests halloc)

  # Extra linking for the project.
  target_link_libraries(runUnitTests voxelyze.3.0)

  # This is so you can do 'make test' to see all your tests run, instead of
  # manually running the executable runUnitTests to see those specific tests.
  add_test(NAME that-test-I-made COMMAND runUnitTests)

  # You can also omit NAME and COMMAND. The second argument could be some other
  # test executable.
  add_test(that-other-test-I-made runUnitTests)

  set_target_properties(runUnitTests PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

  target_compile_options(runUnitTests PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:-G ${ARCH_FLAGS}>) # always with -G
endif()

set(CMAKE_VERBOSE_MAKEFILE OFF) #turn on if we want to debug in compilation.
