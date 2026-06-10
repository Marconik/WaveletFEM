cmake_minimum_required(VERSION 3.20)
project(WaveletFEM VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---- 编译选项 ----
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    add_compile_options(-Wall -Wextra -O3 -march=native)
endif()

# ---- Eigen ----
# 优先使用 external/ 自带版本
set(EIGEN3_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/external/eigen-3.4.0")
if(NOT EXISTS "${EIGEN3_INCLUDE_DIR}/Eigen/Core")
    # 回退到系统安装的 Eigen
    find_package(Eigen3 QUIET)
    if(Eigen3_FOUND)
        set(EIGEN3_INCLUDE_DIR ${EIGEN3_INCLUDE_DIRS})
        message(STATUS "Using system Eigen: ${EIGEN3_INCLUDE_DIR}")
    else()
        message(FATAL_ERROR "Eigen3 not found. Please install Eigen or place it in external/eigen-3.4.0/")
    endif()
else()
    message(STATUS "Using bundled Eigen: ${EIGEN3_INCLUDE_DIR}")
endif()

# ---- 可执行文件 ----
file(GLOB_RECURSE SOURCES src/*.cpp)
add_executable(wavelet_fem ${SOURCES})
target_include_directories(wavelet_fem PRIVATE src ${EIGEN3_INCLUDE_DIR})

# GNUPlot 为运行时外部调用，无需链接

# ---- 信息输出 ----
message(STATUS "WaveletFEM configured successfully")
message(STATUS "  C++ Standard: 17")
message(STATUS "  Eigen include: ${EIGEN3_INCLUDE_DIR}")
