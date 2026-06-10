# ============================================================
# WaveletFEM 依赖检查脚本 (Windows PowerShell)
# 用法: .\scripts\check_deps.ps1
# ============================================================

$ErrorActionPreference = "Continue"
$Pass = 0
$Fail = 0
$Warn = 0

function Write-OK($msg) {
    Write-Host "[OK] " -NoNewline -ForegroundColor Green
    Write-Host $msg
    $script:Pass++
}

function Write-FAIL($msg) {
    Write-Host "[FAIL] " -NoNewline -ForegroundColor Red
    Write-Host $msg
    $script:Fail++
}

function Write-WARN($msg) {
    Write-Host "[WARN] " -NoNewline -ForegroundColor Yellow
    Write-Host $msg
    $script:Warn++
}

Write-Host "=== WaveletFEM Dependency Check ==="
Write-Host ""

# ---- Project root ----
$ProjectDir = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $ProjectDir

# ---- C++ Compiler ----
Write-Host -NoNewline "Checking C++ compiler... "
$clPath = Get-Command cl.exe -ErrorAction SilentlyContinue
$gppPath = Get-Command g++.exe -ErrorAction SilentlyContinue
$clangPath = Get-Command clang++.exe -ErrorAction SilentlyContinue

if ($clPath) {
    $clVer = & cl.exe 2>&1 | Select-Object -First 1
    Write-OK "MSVC compiler found: $clVer"
} elseif ($gppPath) {
    $gppVer = & g++.exe --version 2>&1 | Select-Object -First 1
    Write-OK "GCC compiler found: $gppVer"
} elseif ($clangPath) {
    $clangVer = & clang++.exe --version 2>&1 | Select-Object -First 1
    Write-OK "Clang compiler found: $clangVer"
} else {
    Write-FAIL "No C++ compiler found. Install Visual Studio 2022 or MinGW."
}

# ---- CMake ----
Write-Host -NoNewline "Checking CMake... "
$cmakePath = Get-Command cmake.exe -ErrorAction SilentlyContinue
if ($cmakePath) {
    $cmakeVer = & cmake.exe --version 2>&1 | Select-Object -First 1
    if ($cmakeVer -match '(\d+)\.(\d+)') {
        $major = [int]$Matches[1]
        $minor = [int]$Matches[2]
        if ($major -gt 3 -or ($major -eq 3 -and $minor -ge 20)) {
            Write-OK "CMake $cmakeVer"
        } else {
            Write-FAIL "CMake $cmakeVer (need >= 3.20)"
        }
    } else {
        Write-OK "CMake $cmakeVer"
    }
} else {
    Write-FAIL "CMake not found. Install from https://cmake.org/download/"
}

# ---- GNUPlot ----
Write-Host -NoNewline "Checking GNUPlot... "
$gpPath = Get-Command gnuplot.exe -ErrorAction SilentlyContinue
if ($gpPath) {
    $gpVer = & gnuplot.exe --version 2>&1
    Write-OK "GNUPlot $gpVer"
} else {
    Write-WARN "GNUPlot not found. GIF generation will fail."
    Write-Host "         Install: winget install gnuplot.gnuplot"
    Write-Host "         Or: https://sourceforge.net/projects/gnuplot/"
}

# ---- Eigen 3 ----
Write-Host -NoNewline "Checking Eigen... "
$eigenCore = "external\eigen-3.4.0\Eigen\Core"
if (Test-Path $eigenCore) {
    Write-OK "Eigen: $eigenCore"
} else {
    Write-FAIL "Eigen Core missing: $eigenCore"
    Write-Host "         Download from: https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip"
    Write-Host "         Extract to: external\eigen-3.4.0\"
}

# ---- Summary ----
Write-Host ""
Write-Host "=== Summary ==="
Write-Host "Passed:  $Pass" -ForegroundColor Green
if ($Warn -gt 0) { Write-Host "Warnings: $Warn" -ForegroundColor Yellow }
if ($Fail -gt 0) { Write-Host "Failed:  $Fail" -ForegroundColor Red }
Write-Host ""

if ($Fail -eq 0) {
    Write-Host "All dependencies satisfied." -ForegroundColor Green
    Write-Host "Run: cmake -B build -S . -DCMAKE_BUILD_TYPE=Release"
    Write-Host "     cmake --build build --config Release"
} else {
    Write-Host "Please fix the failed checks above before building." -ForegroundColor Red
}

# Keep window open if double-clicked
if ($Host.Name -eq "ConsoleHost") {
    Write-Host "`nPress any key to exit..."
    $null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")
}
