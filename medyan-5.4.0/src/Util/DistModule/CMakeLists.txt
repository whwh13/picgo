cmake_minimum_required(VERSION 3.7)
project(simd_dist_moduleV2)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Og -march=native")
#-O3 -D_DEBUG -ffast-math -march=native")
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++ -pthread")

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/doc)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/examples)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/plugins)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/utilities)

set(CMAKE_CXX_STANDARD 14)
set(SRC_DIR #ENTER FULL PATH#)
include_directories(#ENTER FULL PATH TO UMESIMD#)

#file(GLOB_RECURSE simd_dist_module_SRC
    #    ${SRC_DIR}/umesimd/UMESimd.h
    #    "umesimd/*.cpp"
set(simd_dist_module_SRC
        ${SRC_DIR}/dist_avx_aux.h
        ${SRC_DIR}/dist_avx.h
        ${SRC_DIR}/dist_avx_par.cpp
        ${SRC_DIR}/dist_avx_par.h
        ${SRC_DIR}/dist_bench.cpp
        ${SRC_DIR}/dist_common.h
        ${SRC_DIR}/dist_coords.h
        ${SRC_DIR}/dist_driver.h
        ${SRC_DIR}/dist_bench.cpp
        ${SRC_DIR}/dist_example.cpp
        ${SRC_DIR}/dist_main.cpp
        ${SRC_DIR}/dist_mod_vars.cpp
	${SRC_DIR}/dist_out.h
	${SRC_DIR}/dist_serial.h
	${SRC_DIR}/dist_simd.h
	${SRC_DIR}/dist_simd_utils.h
	${SRC_DIR}/dist_test.cpp
        #${SRC_DIR}/umesimd/UMEBasicTypes.h
        #${SRC_DIR}/umesimd/UMESimdInterfaceFunctions.h
        #${SRC_DIR}/umesimd/UMESimdTraits.h
        #${SRC_DIR}/umesimd/UMEInline.h
        #${SRC_DIR}/umesimd/UMESimdInterface.h
        #${SRC_DIR}/umesimd/UMESimdVectorEmulation.h
        #${SRC_DIR}/umesimd/UMEMemory.h
        #${SRC_DIR}/umesimd/UMESimdScalarEmulation.h
        ${SRC_DIR}/umesimd/UMESimd.h
        #${SRC_DIR}/umesimd/UMESimdScalarOperators.h
        )

#install(DIRECTORY doc DESTINATION include/ume/simd/doc)
#install(DIRECTORY examples DESTINATION include/ume/simd/examples)
#install(DIRECTORY plugins DESTINATION include/ume/simd/plugins)
#install(DIRECTORY utilities DESTINATION include/ume/simd/utilities)
#install(FILES ${top_files} DESTINATION include/ume/simd)
#include_directories(${SRC_DIR}/umesimd/)
#include(CPack)
add_executable(simd_dist_moduleV2 ${simd_dist_module_SRC})
