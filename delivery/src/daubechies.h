#ifndef WAVELET_FEM_DAUBECHIES_H
#define WAVELET_FEM_DAUBECHIES_H

#include <vector>

// ============================================================
// Daubechies 小波滤波器系数与尺度函数求值
// ============================================================

/// 生成 N 阶 Daubechies 小波的滤波器系数 {h_k}, k = 0..2N-1
/// N 必须是偶数，建议取 2, 4, 6, 8, 10
std::vector<double> daubechies_filter_coefficients(int N);

/// 利用双尺度方程递推计算尺度函数 phi(x) 在整数点的值
/// phi_ints[k] = phi(k),  k = 0..2N-1
/// dphi_ints[k] = phi'(k), k = 0..2N-1（使用有限差分近似）
void evaluate_scaling_function(int N,
                               std::vector<double>& phi_ints,
                               std::vector<double>& dphi_ints);

/// 通过级联算法在细网格上计算 phi(x) 的值
/// 返回在间距 h = 2^{-level} 上的采样值
/// 结果数组大小 = (2N-1) * 2^level + 1
/// phi_samples[i] ≈ phi(i * 2^{-level})
std::vector<double> cascade_scaling_function(int N, int level);

/// 在级联网格上计算 phi''(x) 的值（四阶中心差分）
/// 结果数组大小与 cascade_scaling_function 相同
std::vector<double> cascade_scaling_function_d2(int N, int level);

/// 用 Latto–Resnikoff–Tenenbaum (LRT) 方法精确计算全直线二阶连接系数
///   Γ_d = ∫_{-∞}^{∞} φ'(x) φ'(x - d) dx,  d = -(2N-2) .. (2N-2)
/// 通过双尺度方程导出特征方程 Γ = T Γ（特征值 1），再用矩条件
///   Σ_d d² Γ_d = -2 唯一定标。完全不依赖 φ 的逐点采样或数值积分。
/// 返回长度 4N-3 的向量，下标 i 对应 d = i - (2N-2)。
std::vector<double> connection_coeff_d1d1(int N);

#endif // WAVELET_FEM_DAUBECHIES_H
