#!/bin/bash
# ============================================================
# WaveletFEM 依赖检查脚本 (macOS / Linux / Git Bash on Windows)
# 用法: bash scripts/check_deps.sh
# ============================================================

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

PASS=0
FAIL=0
WARN=0

check_ok()   { echo -e "${GREEN}[OK]${NC} $1"; PASS=$((PASS+1)); }
check_fail() { echo -e "${RED}[FAIL]${NC} $1"; FAIL=$((FAIL+1)); }
check_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; WARN=$((WARN+1)); }

echo "=== WaveletFEM Dependency Check ==="
echo ""

# ---- 项目根目录 ----
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_DIR"

# ---- C++ 编译器 ----
echo -n "Checking C++ compiler... "
if command -v g++ &>/dev/null; then
    CXX=g++
elif command -v clang++ &>/dev/null; then
    CXX=clang++
else
    check_fail "No C++ compiler found (g++ or clang++)"
    CXX=""
fi

if [ -n "$CXX" ]; then
    CXX_VER=$("$CXX" --version 2>&1 | head -1)
    # Check C++17 support
    if echo '#include <iostream>' | "$CXX" -std=c++17 -x c++ -c - -o /dev/null 2>/dev/null; then
        check_ok "C++ compiler: $CXX_VER (C++17 supported)"
    else
        check_fail "C++ compiler: $CXX_VER (C++17 NOT supported)"
    fi
fi

# ---- CMake ----
echo -n "Checking CMake... "
if command -v cmake &>/dev/null; then
    CMAKE_VER=$(cmake --version 2>&1 | head -1)
    CMAKE_MAJOR=$(echo "$CMAKE_VER" | grep -oE '[0-9]+\.[0-9]+' | head -1 | cut -d. -f1)
    CMAKE_MINOR=$(echo "$CMAKE_VER" | grep -oE '[0-9]+\.[0-9]+' | head -1 | cut -d. -f2)
    if [ "$CMAKE_MAJOR" -gt 3 ] || ([ "$CMAKE_MAJOR" -eq 3 ] && [ "$CMAKE_MINOR" -ge 20 ]); then
        check_ok "CMake $CMAKE_VER"
    else
        check_fail "CMake $CMAKE_VER (need >= 3.20)"
    fi
else
    check_fail "CMake not found. Install from https://cmake.org/download/"
fi

# ---- GNUPlot ----
echo -n "Checking GNUPlot... "
if command -v gnuplot &>/dev/null; then
    GNUPLOT_VER=$(gnuplot --version 2>&1)
    GNUPLOT_MAJOR=$(echo "$GNUPLOT_VER" | grep -oE '[0-9]+' | head -1)
    if [ "$GNUPLOT_MAJOR" -ge 5 ]; then
        check_ok "GNUPlot $GNUPLOT_VER"
    else
        check_warn "GNUPlot $GNUPLOT_VER (>= 5.0 recommended)"
    fi
else
    check_warn "GNUPlot not found. GIF generation will fail."
    echo "         Install: brew install gnuplot / apt install gnuplot"
fi

# ---- Eigen 3 ----
echo -n "Checking Eigen... "
EIGEN_CORE="external/eigen-3.4.0/Eigen/Core"
if [ -f "$EIGEN_CORE" ]; then
    check_ok "Eigen: $EIGEN_CORE"
else
    check_fail "Eigen Core missing: $EIGEN_CORE"
    echo "         Run: git add -f $EIGEN_CORE"
    echo "         Or download from: https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip"
fi

# ---- 编译测试 ----
echo -n "Running Eigen compile test... "
cat > /tmp/waveletfem_eigen_test.cpp << 'TESTEOF'
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
int main() {
    Eigen::MatrixXd m(2,2); m.setZero();
    Eigen::SparseMatrix<double> s(2,2);
    return 0;
}
TESTEOF

if "$CXX" -std=c++17 -I external/eigen-3.4.0 /tmp/waveletfem_eigen_test.cpp \
   -o /tmp/waveletfem_eigen_test 2>/dev/null; then
    check_ok "Eigen compiles and links correctly"
    rm -f /tmp/waveletfem_eigen_test /tmp/waveletfem_eigen_test.cpp
else
    check_fail "Eigen compile test failed"
    rm -f /tmp/waveletfem_eigen_test.cpp
fi

# ---- 汇总 ----
echo ""
echo "=== Summary ==="
echo -e "${GREEN}Passed:${NC} $PASS"
if [ $WARN -gt 0 ]; then echo -e "${YELLOW}Warnings:${NC} $WARN"; fi
if [ $FAIL -gt 0 ]; then echo -e "${RED}Failed:${NC} $FAIL"; fi
echo ""

if [ $FAIL -eq 0 ]; then
    echo -e "${GREEN}All dependencies satisfied.${NC}"
    echo "Run: cmake -B build -S . -DCMAKE_BUILD_TYPE=Release"
    echo "     cmake --build build"
    exit 0
else
    echo -e "${RED}Please fix the failed checks above before building.${NC}"
    exit 1
fi
