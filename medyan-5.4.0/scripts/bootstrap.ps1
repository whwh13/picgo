# This script might change the current directory

$ErrorActionPreference = 'Stop'

if(-NOT $medyan_root_dir) {
    Write-Error "Variable medyan_root_dir is not specified"
    exit 1
}

# Set directories
$medyan_build_dir = "$medyan_root_dir\build"
$script_dir = "$medyan_root_dir\scripts"
$build_dir = "$script_dir\.build"
$vcpkg_dir = "$build_dir\vcpkg"

# Set variables
$medyan_vcpkg_cmake_toolchain = "$medyan_root_dir\scripts\.build\vcpkg\scripts\buildsystems\vcpkg.cmake"

if($MEDYAN_NO_GUI -eq "true") {
    $medyan_cmake_no_gui = "-DMEDYAN_NO_GUI=true"
    $medyan_no_gui = $true
} else {
    $medyan_no_gui = $false
}

if($MEDYAN_ADDITIONAL_LINK_DIRS) {
    $medyan_cmake_additional_link_dirs = "-DMEDYAN_ADDITIONAL_LINK_DIRS=$MEDYAN_ADDITIONAL_LINK_DIRS"
}

if("$MEDYAN_RPATH") {
    $medyan_cmake_rpath = "-DMEDYAN_RPATH=$MEDYAN_RPATH"
}


# Download and setup vcpkg
Function Install-Vcpkg([bool]$required, [bool]$rebuild) {

    if($required) {
        if(-NOT $(Test-Path $vcpkg_dir) -OR $rebuild) {
            Write-Host "Configuring vcpkg..."
            Set-Location $build_dir
            git clone https://github.com/Microsoft/vcpkg.git
            Set-Location $vcpkg_dir
            git checkout c9e786d81a890ef6b3932779925f11e696dc9541
            .\bootstrap-vcpkg.bat
        } else {
            Write-Host "vcpkg is already installed."
        }
    } else {
        Write-Host "Skipping vcpkg installation."
    }

}

# Setup dependencies
Function Install-VcpkgPackages([bool]$no_gui) {

    $Env:VCPKG_DEFAULT_TRIPLET="x64-windows"

    Set-Location $vcpkg_dir
    .\vcpkg install catch2 eigen3 fmt highfive spdlog spectra xtensor boost-pool boost-heap boost-ublas boost-range --binarysource=clear
    if(-Not $no_gui) {
        .\vcpkg install glfw3 glad glm imgui[opengl3-binding,glfw-binding] implot nativefiledialog stb --binarysource=clear
    }
}

# Generate using CMake
Function Use-Cmake() {

    Write-Host "Creating build files..."
    mkdir -Force $medyan_build_dir
    Set-Location $medyan_build_dir
    cmake `
        $medyan_cmake_no_gui `
        $medyan_cmake_additional_link_dirs `
        $medyan_cmake_rpath `
        $MEDYAN_CMAKE_EXTRA_ARGS `
        .. "-DCMAKE_TOOLCHAIN_FILE=$medyan_vcpkg_cmake_toolchain"
}

# Make directories
mkdir -Force $build_dir

# Use vcpkg to resolve dependencies
Install-Vcpkg $true $false
Install-VcpkgPackages $medyan_no_gui

# Use CMake to generate build files
Use-Cmake
