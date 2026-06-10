#include "examples.h"
#include <cmath>

// ============================================================
// 算例 1: Robin 边界条件
// u(x,y,t) = e^{-2x} * sin(πy/2) * (t+1)/(t+2)
// ∂u/∂n + u = g on boundary
// a₁ = a₂ = 1
// ============================================================

ProblemConfig Example1::config() {
    ProblemConfig cfg;
    cfg.ax = 0.0; cfg.bx = 2.0;
    cfg.ay = 0.0; cfg.by = 2.0;
    cfg.a1 = 1.0; cfg.a2 = 1.0;

    // Robin BC: ∂u/∂n + u = g → b_i = 1, d_i = 1
    cfg.b = {1.0, 1.0, 1.0, 1.0};
    cfg.d = {1.0, 1.0, 1.0, 1.0};

    cfg.N = 6;
    cfg.j = 4;
    cfg.T = 1.0;
    cfg.dt0 = 0.005;
    cfg.variable_dt = true;

    return cfg;
}

double Example1::source(double x, double y, double t) {
    // f = ∂u/∂t - ∂²u/∂x² - ∂²u/∂y²
    // u = e^{-2x} * sin(πy/2) * (t+1)/(t+2)
    //
    // ∂u/∂t = e^{-2x} * sin(πy/2) * 1/(t+2)²
    // ∂²u/∂x² = 4 * e^{-2x} * sin(πy/2) * (t+1)/(t+2)
    // ∂²u/∂y² = -(π²/4) * e^{-2x} * sin(πy/2) * (t+1)/(t+2)

    double factor = std::exp(-2.0 * x) * std::sin(M_PI * y / 2.0);
    double time_term = (t + 1.0) / (t + 2.0);
    double dt_term = 1.0 / ((t + 2.0) * (t + 2.0));

    return factor * (dt_term - (4.0 - M_PI * M_PI / 4.0) * time_term);
}

double Example1::initial(double x, double y) {
    return exact(x, y, 0.0);
}

double Example1::exact(double x, double y, double t) {
    return std::exp(-2.0 * x) * std::sin(M_PI * y / 2.0)
           * (t + 1.0) / (t + 2.0);
}

double Example1::g_left(double x, double y, double t) {
    // x = 0, ∂u/∂n = -∂u/∂x
    // g = -∂u/∂x + u (evaluated at x=0)
    // ∂u/∂x = -2 * e^{-2x} * sin(πy/2) * (t+1)/(t+2)
    // At x=0: ∂u/∂x = -2 * sin(πy/2) * (t+1)/(t+2)
    // -∂u/∂x + u = 2*sin(πy/2)*(t+1)/(t+2) + sin(πy/2)*(t+1)/(t+2)
    // = 3*sin(πy/2)*(t+1)/(t+2)
    (void)x;
    return 3.0 * std::sin(M_PI * y / 2.0) * (t + 1.0) / (t + 2.0);
}

double Example1::g_right(double x, double y, double t) {
    // x = 2, ∂u/∂n = ∂u/∂x
    // g = ∂u/∂x + u (evaluated at x=2)
    // ∂u/∂x = -2 * e^{-4} * sin(πy/2) * (t+1)/(t+2)
    // u = e^{-4} * sin(πy/2) * (t+1)/(t+2)
    // g = (-2*e^{-4} + e^{-4}) * sin(πy/2) * (t+1)/(t+2)
    //   = -e^{-4} * sin(πy/2) * (t+1)/(t+2)
    (void)x;
    return -std::exp(-4.0) * std::sin(M_PI * y / 2.0) * (t + 1.0) / (t + 2.0);
}

double Example1::g_bottom(double x, double y, double t) {
    // y = 0, ∂u/∂n = -∂u/∂y
    // g = -∂u/∂y + u (evaluated at y=0)
    // u(x,0,t) = e^{-2x} * 0 * (t+1)/(t+2) = 0
    // ∂u/∂y = e^{-2x} * (π/2)*cos(πy/2) * (t+1)/(t+2)
    // At y=0: ∂u/∂y = (π/2) * e^{-2x} * (t+1)/(t+2)
    // g = -(π/2) * e^{-2x} * (t+1)/(t+2)
    (void)y;
    return -M_PI / 2.0 * std::exp(-2.0 * x) * (t + 1.0) / (t + 2.0);
}

double Example1::g_top(double x, double y, double t) {
    // y = 2, ∂u/∂n = ∂u/∂y
    // u(x,2,t) = e^{-2x} * sin(π) * (t+1)/(t+2) = 0
    // ∂u/∂y = e^{-2x} * (π/2)*cos(π) * (t+1)/(t+2)
    // = -e^{-2x} * (π/2) * (t+1)/(t+2)
    // g = ∂u/∂y + u = -e^{-2x} * (π/2) * (t+1)/(t+2)
    (void)y;
    return -M_PI / 2.0 * std::exp(-2.0 * x) * (t + 1.0) / (t + 2.0);
}

// ============================================================
// 算例 2: Neumann 边界条件
// u(x,y,t) = exp(-((x+1)+(y+1))/(2t+1))
// ∂u/∂n = g on boundary
// a₁ = a₂ = 1
// ============================================================

ProblemConfig Example2::config() {
    ProblemConfig cfg;
    cfg.ax = 0.0; cfg.bx = 2.0;
    cfg.ay = 0.0; cfg.by = 2.0;
    cfg.a1 = 1.0; cfg.a2 = 1.0;

    // Neumann BC: ∂u/∂n = g → b_i = 1, d_i = 0
    cfg.b = {1.0, 1.0, 1.0, 1.0};
    cfg.d = {0.0, 0.0, 0.0, 0.0};

    cfg.N = 6;
    cfg.j = 4;
    cfg.T = 1.0;
    cfg.dt0 = 0.005;
    cfg.variable_dt = true;

    return cfg;
}

double Example2::source(double x, double y, double t) {
    // u = exp(-((x+1)+(y+1))/(2t+1))
    // Let s = 2t+1, r = (x+1)+(y+1)
    // u = exp(-r/s)
    //
    // ∂u/∂t = u * r/(s²) * 2 = 2r/s² * u
    // ∂²u/∂x² = u * (-1/s)² = u/s²
    // ∂²u/∂y² = u/s²
    //
    // f = ∂u/∂t - ∂²u/∂x² - ∂²u/∂y²
    //   = u * (2r/s² - 2/s²) = u * 2(r-1)/s²

    double s = 2.0 * t + 1.0;
    double r = (x + 1.0) + (y + 1.0);
    double u = std::exp(-r / s);

    return u * 2.0 * (r - 1.0) / (s * s);
}

double Example2::initial(double x, double y) {
    return exact(x, y, 0.0);
}

double Example2::exact(double x, double y, double t) {
    double r = (x + 1.0) + (y + 1.0);
    double s = 2.0 * t + 1.0;
    return std::exp(-r / s);
}

double Example2::g_left(double x, double y, double t) {
    // x = 0, ∂u/∂n = -∂u/∂x
    // ∂u/∂x = -u/s
    // g = -∂u/∂x = u/s
    double s = 2.0 * t + 1.0;
    double u = exact(0.0, y, t);
    (void)x;
    return u / s;
}

double Example2::g_right(double x, double y, double t) {
    // x = 2, ∂u/∂n = ∂u/∂x = -u/s
    double s = 2.0 * t + 1.0;
    double u = exact(2.0, y, t);
    (void)x;
    return -u / s;
}

double Example2::g_bottom(double x, double y, double t) {
    // y = 0, ∂u/∂n = -∂u/∂y
    // ∂u/∂y = -u/s, g = -∂u/∂y = u/s
    double s = 2.0 * t + 1.0;
    double u = exact(x, 0.0, t);
    (void)y;
    return u / s;
}

double Example2::g_top(double x, double y, double t) {
    // y = 2, ∂u/∂n = ∂u/∂y = -u/s
    double s = 2.0 * t + 1.0;
    double u = exact(x, 2.0, t);
    (void)y;
    return -u / s;
}
