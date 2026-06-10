#ifndef WAVELET_FEM_EXAMPLES_H
#define WAVELET_FEM_EXAMPLES_H

#include "types.h"

// ============================================================
// 内置算例
// ============================================================

/// 算例 1: Robin 边界条件
/// 精确解: u = e^{-2x} sin(πy/2) * (t+1)/(t+2)
/// 边界类型: ∂u/∂n + u = g (Robin)
struct Example1 {
    static ProblemConfig config();

    // PDE: ∂²u/∂x² + ∂²u/∂y² + f = ∂u/∂t
    static double source(double x, double y, double t);

    static double initial(double x, double y);
    static double exact(double x, double y, double t);

    // 边界函数 g_i
    static double g_left(double x, double y, double t);
    static double g_right(double x, double y, double t);
    static double g_bottom(double x, double y, double t);
    static double g_top(double x, double y, double t);
};

/// 算例 2: Neumann 边界条件
/// 精确解: u = exp(-((x+1)+(y+1))/(2t+1))
/// 边界类型: ∂u/∂n = g (Neumann)
struct Example2 {
    static ProblemConfig config();

    static double source(double x, double y, double t);
    static double initial(double x, double y);
    static double exact(double x, double y, double t);

    static double g_left(double x, double y, double t);
    static double g_right(double x, double y, double t);
    static double g_bottom(double x, double y, double t);
    static double g_top(double x, double y, double t);
};

#endif // WAVELET_FEM_EXAMPLES_H
