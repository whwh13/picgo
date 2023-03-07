cmake_minimum_required(VERSION 3.7)

project(simd_dist_module)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -D_DEBUG -ffast-math -march=native")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/doc)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/examples)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/plugins)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/utilities)
set(CMAKE_CXX_STANDARD 14)
set(SRC_DIR /home/aravind/Flory_stableversionMar26_18_decoupled_functionalserial/medyan-private/src/)

file(GLOB simd_dist_module
        "umesimd/*.h"
        "umesimd/*.cpp"

        )

install(DIRECTORY doc DESTINATION include/ume/simd/doc)
install(DIRECTORY examples DESTINATION include/ume/simd/examples)
install(DIRECTORY plugins DESTINATION include/ume/simd/plugins)
install(DIRECTORY utilities DESTINATION include/ume/simd/utilities)
install(FILES ${top_files} DESTINATION include/ume/simd)

include(CPack)