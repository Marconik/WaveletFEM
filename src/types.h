#ifndef WAVELET_FEM_TYPES_H
#define WAVELET_FEM_TYPES_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <array>
#include <vector>
#include <functional>
#include <string>

// ============================================================
// 全局类型别名
// ============================================================

using VectorXd = Eigen::VectorXd;
using MatrixXd = Eigen::MatrixXd;
using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

// 标量函数类型
// ScalarFunction2D: f(x, y)
// ScalarFunction3D: f(x, y, t)
using ScalarFunction2D = std::function<double(double, double)>;
using ScalarFunction3D = std::function<double(double, double, double)>;

// ============================================================
// 问题配置结构体
// ============================================================

struct ProblemConfig {
    // 求解区域 Ω = [ax, bx] × [ay, by]
    double ax = 0.0, bx = 2.0;
    double ay = 0.0, by = 2.0;

    // PDE 系数
    double a1 = 1.0;   // ∂²u/∂x² 的系数
    double a2 = 1.0;   // ∂²u/∂y² 的系数

    // 边界条件系数: b_i * ∂u/∂n + d_i * u = g_i
    // 边顺序: 0=左(x=ax), 1=右(x=bx), 2=下(y=ay), 3=上(y=by)
    std::array<double, 4> b = {1.0, 1.0, 1.0, 1.0};
    std::array<double, 4> d = {1.0, 1.0, 1.0, 1.0};

    // Daubechies 小波参数
    int N = 6;          // 阶数（建议 4, 6, 8, 10）
    int j = 4;          // 分辨率级别

    // 时间参数
    double T = 1.0;     // 终止时间
    double dt0 = 0.005; // 初始时间步长
    bool variable_dt = true;  // 是否使用变步长策略

    // 可视化参数
    int nviz = 101;     // 可视化网格点数（每方向）
};

// ============================================================
// 时间步记录
// ============================================================

struct TimeStepRecord {
    double t;                    // 当前时间
    double dt;                   // 本步使用的步长
    MatrixXd U;                  // n × n 系数矩阵
    double l2_error = -1.0;      // L² 相对误差（-1 表示未计算）
};

#endif // WAVELET_FEM_TYPES_H
