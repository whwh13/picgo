cmake_minimum_required(VERSION 3.13)

#######################################
# Configurations
#######################################

# Get git branch.
execute_process(
    COMMAND           git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE   GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
# Get the latest commit hash.
execute_process(
    COMMAND           git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE   GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
# Get latest annotated tag.
execute_process(
    COMMAND           git describe --abbrev=0 --tags
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE   GIT_LATEST_TAG
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# vcpkg triplet config
if(DEFINED ENV{VCPKG_DEFAULT_TRIPLET} AND NOT DEFINED VCPKG_TARGET_TRIPLET)
    set(VCPKG_TARGET_TRIPLET "$ENV{VCPKG_DEFAULT_TRIPLET}" CACHE STRING "")
endif()


project(medyan)

#######################################
# Source files
#######################################

set(src_dir          "src")
set(src_dir_external "external")

# List sources
file(GLOB_RECURSE src_list LIST_DIRECTORIES true CONFIGURE_DEPENDS "${src_dir}/*.cpp")

#######################################
# Compiling configs
#######################################

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Compiler flags
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
    set(CMAKE_CXX_FLAGS "/arch:AVX2 /MP /EHsc /bigobj")
else()
    set(CMAKE_CXX_FLAGS "-Wno-sign-compare -Wno-unused-variable -Wno-reorder -Wno-unused-but-set-variable -Wno-unused-local-typedefs -Wall -ffast-math -fno-finite-math-only -mtune=native -march=native")
    set(CMAKE_CXX_FLAGS_DEBUG "-g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -funroll-loops -flto")
endif()

# Add the executable
add_executable(medyan ${src_list})

# Create the source groups for source tree with root at CMAKE_CURRENT_SOURCE_DIR.
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${src_list})

# Include paths
target_include_directories(medyan PUBLIC
    ${src_dir}
    ${src_dir}/Chemistry
    ${src_dir}/Mechanics
    ${src_dir}/Mechanics/ForceField
    ${src_dir}/Mechanics/ForceField/Boundary
    ${src_dir}/Mechanics/ForceField/Branching
    ${src_dir}/Mechanics/ForceField/Bubble
    ${src_dir}/Mechanics/ForceField/Filament
    ${src_dir}/Mechanics/ForceField/Linker
    ${src_dir}/Mechanics/ForceField/MotorGhost
    ${src_dir}/Mechanics/ForceField/Volume
    ${src_dir}/Mechanics/Minimizer
    ${src_dir}/Structure
    ${src_dir}/Structure/Special
    ${src_dir_external}
)

# Preprocessor macros
target_compile_definitions(medyan PRIVATE
    # Meta: git
    GIT_COMMIT_HASH="${GIT_COMMIT_HASH}"
    GIT_BRANCH="${GIT_BRANCH}"
    GIT_LATEST_TAG="${GIT_LATEST_TAG}"

    # Chemistry algorithm
    TRACK_DEPENDENTS # Needed for Gillespie and NRM algorithm
    # Track zero copy and max copy number.
    # For passivating and activating reactions accordingly
    TRACK_ZERO_COPY_N
    TRACK_UPPER_COPY_N

    # Energy minimization algorithm
    MOVEBEADSLINESEARCH

    # Neighbors list calculations (choose 1 of 3)
    #NLORIGINAL
    #HYBRID_NLSTENCILLIST
    SIMDBINDINGSEARCH

    # Boost memory
    BOOST_MEM_POOL
    BOOL_POOL_NSIZE=65536

    # Debug
    #PLOSFEEDBACK
    CHECKFORCES_INF_NAN
)

if("${MEDYAN_NO_GUI}" STREQUAL "true")
    target_compile_definitions(medyan PRIVATE
        NO_GUI             # Disable GUI
    )
endif()

if(MSVC)
    target_compile_definitions(medyan PRIVATE
        _USE_MATH_DEFINES  # Enable math macros like M_PI
    )
endif()

# Additional link directories
if(MSVC)
else()
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    target_link_libraries(medyan PRIVATE Threads::Threads)

    # GCC <9 requires linking to stdc++fs
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_link_libraries(medyan PRIVATE stdc++fs)
    endif()
endif()

if(MEDYAN_ADDITIONAL_LINK_DIRS)
    target_link_directories(medyan PRIVATE ${MEDYAN_ADDITIONAL_LINK_DIRS})
endif()

# RPath specification
if(MEDYAN_RPATH)
    set_target_properties(medyan PROPERTIES
        BUILD_RPATH ${MEDYAN_RPATH}
    )
endif()

# Build tool specific
if(MSVC)
    set_directory_properties(PROPERTY VS_STARTUP_PROJECT "medyan")
endif()


#######################################
# External dependencies
#######################################

# boost
find_package(Boost REQUIRED COMPONENTS)
target_link_libraries(medyan PRIVATE Boost::boost)

# catch2
find_package(Catch2 CONFIG REQUIRED)
target_link_libraries(medyan PRIVATE Catch2::Catch2)

# eigen3
find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(medyan PRIVATE Eigen3::Eigen)

# fmt
find_package(fmt CONFIG REQUIRED)
target_link_libraries(medyan PRIVATE fmt::fmt)

# highfive
find_package(HighFive CONFIG REQUIRED)
target_link_libraries(medyan PRIVATE libdeps HighFive HighFive_libdeps HighFive_HighFive)

# spectra
find_package(spectra CONFIG REQUIRED)
target_link_libraries(medyan PRIVATE Spectra::Spectra)

# spdlog
find_package(spdlog CONFIG REQUIRED)
target_link_libraries(medyan PRIVATE spdlog::spdlog spdlog::spdlog_header_only)

# xtensor
find_package(xtensor CONFIG REQUIRED)
target_link_libraries(medyan PRIVATE xtensor xtensor::optimize)

# GUI related
if(NOT "${MEDYAN_NO_GUI}" STREQUAL "true")
    # GLFW3
    find_package(glfw3 CONFIG REQUIRED)
    target_link_libraries(medyan PRIVATE glfw)

    # GLAD
    find_package(glad CONFIG REQUIRED)
    target_link_libraries(medyan PRIVATE glad::glad)

    # GLM
    find_package(glm CONFIG REQUIRED)
    target_link_libraries(medyan PRIVATE glm::glm)

    # imgui
    find_package(imgui CONFIG REQUIRED)
    target_link_libraries(medyan PRIVATE imgui::imgui)

    # implot
    find_package(implot CONFIG REQUIRED)
    target_link_libraries(medyan PRIVATE implot::implot)

    # stb
    find_path(STB_INCLUDE_DIRS "stb.h")
    target_include_directories(medyan PRIVATE ${STB_INCLUDE_DIRS})

    # native file dialog
    find_package(unofficial-nativefiledialog CONFIG REQUIRED)
    target_link_libraries(medyan PRIVATE unofficial::nativefiledialog::nfd)

endif()
