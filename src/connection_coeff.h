#ifndef WAVELET_FEM_CONNECTION_COEFF_H
#define WAVELET_FEM_CONNECTION_COEFF_H

#include <vector>
#include <Eigen/Dense>

// ============================================================
// 连接系数预计算
// 通过级联算法 + 高精度数值积分计算质量和刚度矩阵所需的
// 两-term 连接系数 C_k^n(x) 和矩连接系数 M_k^r(x)
// ============================================================

/// 预计算的连接系数表
/// 在级联网格上存储 phi(x) 和 phi''(x) 的值，
/// 支持快速查值和积分
struct ConnectionTable {
    int N;                     // Daubechies 阶数
    int cascade_level;         // 级联算法级别
    int support_len;           // 支撑集整数区间数 = 2N-1
    int n_samples;             // 总采样点数

    double dx;                 // 网格间距 = 2^{-cascade_level}
    std::vector<double> phi;   // phi(i*dx), i=0..n_samples-1
    std::vector<double> d2phi; // phi''(i*dx), i=0..n_samples-1

    /// 构造函数：预计算所有需要的值
    /// @param N Daubechies 阶数
    /// @param cascade_level 级联级别（建议 ≥ 12）
    ConnectionTable(int N, int cascade_level);

    /// 在任意 x 处求值 phi(x)（线性插值）
    double eval_phi(double x) const;

    /// 在任意 x 处求值 phi''(x)（线性插值）
    double eval_d2phi(double x) const;

    /// 在任意 x 处求值连接系数 C_k^n(x) = ∫_0^x phi^{(n)}(y-k) phi(y) dy
    /// @param k 平移指标 k ∈ [-(2N-2), 2N-2]
    /// @param n 导数阶数: 0=质量型, 2=刚度型
    /// @param x 积分上限
    double eval_C(int k, int n, double x) const;

    /// 在任意 x 处求值矩连接系数 M_k^r(x) = ∫_0^x y^r phi(y-k) dy
    /// @param k 平移指标
    /// @param r 矩的阶数
    /// @param x 积分上限
    double eval_M(int k, int r, double x) const;

private:
    /// 使用 Simpson 法则计算积分 ∫_a^b f(i*dx) dx
    double integrate_simpson(const std::vector<double>& f,
                             int i_start, int i_end) const;
};

#endif // WAVELET_FEM_CONNECTION_COEFF_H
