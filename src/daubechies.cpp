#include "daubechies.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <complex>

using namespace Eigen;

// ============================================================
// 已知的 Daubechies 滤波器系数（高精度）
// ============================================================

namespace {

struct DaubCoeffs {
    int N;
    std::vector<double> h;
};

const DaubCoeffs KNOWN_COEFFS[] = {
    {2, {
        0.4829629131445341,
        0.8365163037378079,
        0.2241438680420134,
        -0.1294095225512604
    }},
    {3, {
        0.3326705529500825,
        0.8068915093110924,
        0.4598775021184914,
        -0.1350110200102546,
        -0.0854412738820267,
        0.0352262918857095
    }},
    {4, {
        0.2303778133088964,
        0.7148465705529154,
        0.6308807679398587,
        -0.0279837694168599,
        -0.1870348117190931,
        0.0308413818355607,
        0.0328830116668852,
        -0.0105974017850690
    }},
    {5, {
        0.1601023979741929,
        0.6038292697971895,
        0.7243085284377726,
        0.1384281459013203,
        -0.2422948870663823,
        -0.0322448695846381,
        0.0775714938400459,
        -0.0062414902127983,
        -0.0125807519990820,
        0.0033357252854738
    }},
    {6, {
        0.1115407433501095,
        0.4946238903984533,
        0.7511339080210959,
        0.3152503517091982,
        -0.2262646939654400,
        -0.1297668675672625,
        0.0975016055873225,
        0.0275228655303053,
        -0.0315820393174862,
        0.0005538422011614,
        0.0047772575109455,
        -0.0010773010853085
    }}
};

const int NUM_KNOWN = sizeof(KNOWN_COEFFS) / sizeof(KNOWN_COEFFS[0]);

// ============================================================
// 辅助函数：计算二项式系数
// ============================================================

double binomial(int n, int k) {
    if (k < 0 || k > n) return 0.0;
    if (k == 0 || k == n) return 1.0;
    double result = 1.0;
    for (int i = 1; i <= k; ++i) {
        result = result * (n - k + i) / i;
    }
    return result;
}

// ============================================================
// 通过谱分解计算任意 N 的 Daubechies 滤波器系数
// 参考：Daubechies (1988) "Orthonormal bases of compactly
// supported wavelets"
// ============================================================

std::vector<double> compute_daubechies_coeffs(int N) {
    // Step 1: 构造多项式 P(y) = Σ_{k=0}^{N-1} C(N-1+k, k) y^k
    // 其中 y = sin²(ω/2) = (1-cos ω)/2
    // P(y) 是 |m₀(ω)|² 的多项式部分

    // 构造伴随矩阵求 P(y) 的根
    if (N == 1) {
        return {std::sqrt(2.0) / 2.0, std::sqrt(2.0) / 2.0};
    }

    int deg = N - 1;
    MatrixXd companion(deg, deg);
    companion.setZero();

    if (deg == 1) {
        // P(y) = C(N-1,0) + C(N,1)*y ... for N=2: C(1,0)+C(2,1)*y = 1+2y
        // Root: y = -1/2
        companion(0, 0) = -binomial(N - 1 + 0, 0) / binomial(N - 1 + 1, 1);
    } else {
        double leading = binomial(N - 1 + deg, deg);
        for (int i = 0; i < deg; ++i) {
            if (i < deg - 1) companion(i, i + 1) = 1.0;
            companion(deg - 1, i) = -binomial(N - 1 + i, i) / leading;
        }
    }

    EigenSolver<MatrixXd> solver(companion);
    VectorXcd roots = solver.eigenvalues();

    // Step 2: 将 y-平面根映射到 z-平面
    // y = -(z-1)²/(4z)  →  z² + (4y-2)z + 1 = 0
    std::vector<std::complex<double>> z_roots;
    for (int i = 0; i < deg; ++i) {
        std::complex<double> y = roots(i);
        // 处理数值误差导致的虚部
        if (std::abs(y.imag()) < 1e-12) y = std::complex<double>(y.real(), 0.0);

        // z² + (4y-2)z + 1 = 0
        std::complex<double> b = 4.0 * y - 2.0;
        std::complex<double> disc = std::sqrt(b * b - 4.0);
        std::complex<double> z1 = (-b + disc) / 2.0;
        std::complex<double> z2 = (-b - disc) / 2.0;

        // 选择单位圆内的根（或单位圆上）
        if (std::abs(z1) <= 1.0 + 1e-10) {
            z_roots.push_back(z1);
        } else {
            z_roots.push_back(z2);
        }
    }

    // Step 3: 构造多项式 Q(z) = Π (z - z_k)
    // m₀(z) = ((1+z)/2)^N * Q(z) / Q(1) (归一化)

    // 计算 Q(1) = Π (1 - z_k)
    std::complex<double> Q1(1.0, 0.0);
    for (auto& zk : z_roots) {
        Q1 *= (1.0 - zk);
    }

    // 展开 ((1+z)/2)^N * Q(z)
    // 使用多项式乘法
    // 先展开 ((1+z)/2)^N
    std::vector<std::complex<double>> factor1(N + 1, 0.0);
    for (int k = 0; k <= N; ++k) {
        factor1[k] = binomial(N, k) / std::pow(2.0, N);
    }

    // Q(z) = Π (z - z_k)，从常数项开始
    std::vector<std::complex<double>> Q_poly = {1.0};
    for (auto& zk : z_roots) {
        std::vector<std::complex<double>> new_poly(Q_poly.size() + 1, 0.0);
        for (size_t j = 0; j < Q_poly.size(); ++j) {
            new_poly[j] += Q_poly[j] * (-zk);        // 常数项 × (-zk)
            new_poly[j + 1] += Q_poly[j];             // z × 系数
        }
        Q_poly = new_poly;
    }

    // m₀(z) 的系数 = factor1 * Q_poly / Q(1)
    // 多项式乘法
    int m_deg = N + (int)Q_poly.size() - 1;  // should be 2N-1
    std::vector<std::complex<double>> m_poly(m_deg + 1, 0.0);

    for (int i = 0; i <= N; ++i) {
        for (size_t j = 0; j < Q_poly.size(); ++j) {
            m_poly[i + j] += factor1[i] * Q_poly[j];
        }
    }

    // 除以 Q(1)
    for (auto& c : m_poly) c /= Q1;

    // Step 4: 取实部并乘以 √2 得到滤波器系数
    std::vector<double> h(2 * N);
    for (int k = 0; k < 2 * N; ++k) {
        h[k] = std::sqrt(2.0) * m_poly[k].real();
        // 抑制数值噪声
        if (std::abs(h[k]) < 1e-15) h[k] = 0.0;
    }

    return h;
}

} // anonymous namespace

// ============================================================
// 公共接口
// ============================================================

std::vector<double> daubechies_filter_coefficients(int N) {
    if (N < 1 || N > 20) {
        throw std::invalid_argument("Daubechies order N must be between 1 and 20");
    }

    // 查表获取已知系数
    for (int i = 0; i < NUM_KNOWN; ++i) {
        if (KNOWN_COEFFS[i].N == N) {
            return KNOWN_COEFFS[i].h;
        }
    }

    // 计算未知 N 的系数
    return compute_daubechies_coeffs(N);
}

void evaluate_scaling_function(int N,
                               std::vector<double>& phi_ints,
                               std::vector<double>& dphi_ints) {
    // ============================================================
    // 计算 phi(x) 在整数点的值
    // 通过求解 phi(k) = √2 Σ h_i phi(2k-i) 的特征值问题
    // phi(0) = phi(2N-1) = 0
    // ============================================================

    auto h = daubechies_filter_coefficients(N);
    int M = 2 * N - 2;  // 内部整数点数量: 1, 2, ..., 2N-2

    phi_ints.assign(2 * N, 0.0);   // 包括 phi(0)..phi(2N-1), phi(2N-1)=0
    dphi_ints.assign(2 * N, 0.0);

    if (M <= 0) return;  // N=1 的情况

    // 构造矩阵 A_{k,j} = √2 * h_{2k+1 - (j+1)}
    // 即 A_{k,j} = √2 * h_{2k - j + 1}
    // 其中 k, j = 0..M-1 对应 phi(1)..phi(2N-2)
    MatrixXd A(M, M);
    A.setZero();

    for (int k = 0; k < M; ++k) {
        int xk = k + 1;  // 实际整数坐标
        for (int j = 0; j < M; ++j) {
            int xj = j + 1;  // 实际整数坐标
            int idx = 2 * xk - xj;
            if (idx >= 0 && idx < 2 * N) {
                A(k, j) = std::sqrt(2.0) * h[idx];
            }
        }
    }

    // 求解特征向量（特征值 ≈ 1）
    EigenSolver<MatrixXd> solver(A);
    VectorXcd evals = solver.eigenvalues();
    MatrixXcd evecs = solver.eigenvectors();

    // 找到最接近 1 的特征值
    int best_idx = 0;
    double best_dist = std::abs(evals(0) - 1.0);
    for (int i = 1; i < M; ++i) {
        double dist = std::abs(evals(i) - 1.0);
        if (dist < best_dist) {
            best_dist = dist;
            best_idx = i;
        }
    }

    // 提取实部特征向量并归一化 Σ φ(k) = 1
    VectorXd phi_vec(M);
    for (int i = 0; i < M; ++i) {
        phi_vec(i) = evecs(i, best_idx).real();
    }

    // 确保非负（特征向量可能整体翻转符号）
    double sum_phi = phi_vec.sum();
    if (sum_phi < 0) phi_vec = -phi_vec;
    sum_phi = phi_vec.sum();
    if (std::abs(sum_phi) > 1e-12) {
        phi_vec /= sum_phi;
    }

    // 填充 phi_ints: [phi(0), phi(1), ..., phi(2N-2)]
    for (int i = 0; i < M; ++i) {
        phi_ints[i + 1] = phi_vec(i);
    }

    // ============================================================
    // 计算 phi'(x) 在整数点的值
    // 通过求解类似的线性系统，使用额外的矩条件
    // phi'(x) = 2√2 Σ h_i phi'(2x-i)
    // 矩条件: Σ k * phi'(k) = -1 (由 ∫ φ'(x) dx 在支撑集上的积分)
    // ============================================================

    // 同样构造矩阵，但这次是对 phi' 的齐次系统
    // 特征值应为 1/2（因为 phi'(x) = 2√2 Σ h_i phi'(2x-i)，在整数点
    // phi'(k) = 2√2 Σ h_i phi'(2k-i)，所以特征值是 1/2）

    MatrixXd B(M, M);
    B.setZero();
    for (int k = 0; k < M; ++k) {
        int xk = k + 1;
        for (int j = 0; j < M; ++j) {
            int xj = j + 1;
            int idx = 2 * xk - xj;
            if (idx >= 0 && idx < 2 * N) {
                B(k, j) = 2.0 * std::sqrt(2.0) * h[idx];
            }
        }
    }

    EigenSolver<MatrixXd> solverB(B);
    VectorXcd evalsB = solverB.eigenvalues();
    MatrixXcd evecsB = solverB.eigenvectors();

    // 找最接近 1 的特征值（phi' 的特征值也应为 1？让我检查...
    // phi'(k) = 2√2 Σ h_i phi'(2k-i)
    // 这与 phi 相同，所以特征值也是 1）
    best_idx = 0;
    best_dist = std::abs(evalsB(0) - 1.0);
    for (int i = 1; i < M; ++i) {
        double dist = std::abs(evalsB(i) - 1.0);
        if (dist < best_dist) {
            best_dist = dist;
            best_idx = i;
        }
    }

    VectorXd dphi_vec(M);
    for (int i = 0; i < M; ++i) {
        dphi_vec(i) = evecsB(i, best_idx).real();
    }

    // 归一化: 使用矩条件 Σ k * phi'(k) = -1
    double sum_k_dphi = 0.0;
    for (int i = 0; i < M; ++i) {
        sum_k_dphi += (i + 1) * dphi_vec(i);
    }
    if (std::abs(sum_k_dphi) > 1e-12) {
        dphi_vec /= -sum_k_dphi;
    }

    for (int i = 0; i < M; ++i) {
        dphi_ints[i + 1] = dphi_vec(i);
    }
}

std::vector<double> cascade_scaling_function(int N, int level) {
    // ============================================================
    // 标准级联算法：从 Kronecker delta 出发迭代应用胀缩方程
    // φ_{m+1}(x) = √2 Σ_k h_k φ_m(2x - k)
    // 经过足够多次迭代后收敛到真实的 φ
    // ============================================================

    auto h = daubechies_filter_coefficients(N);
    int L = 2 * N;
    double sqrt2 = std::sqrt(2.0);
    int support_len = 2 * N - 1;

    // 从 Kronecker delta 开始（初始化为 δ 函数在原点）
    std::vector<double> cur = {1.0};

    // 级联迭代次数：需足够多以收敛
    int cascade_iters = level + 12;  // 额外 12 次确保收敛

    for (int iter = 0; iter < cascade_iters; ++iter) {
        int nc = (int)cur.size();
        // 上采样（插入零）
        int n_up = 2 * nc;
        std::vector<double> up(n_up, 0.0);
        for (int i = 0; i < nc; ++i) up[2 * i] = cur[i];

        // 卷积滤波
        int n_next = n_up + L - 1;
        std::vector<double> next(n_next, 0.0);
        for (int n = 0; n < n_next; ++n) {
            double sum = 0.0;
            for (int k = 0; k < L; ++k) {
                int src = n - k;
                if (src >= 0 && src < n_up)
                    sum += h[k] * up[src];
            }
            next[n] = sqrt2 * sum;
        }
        cur = next;
    }

    // 级联收敛后，phi 的中心位于约 (L-1)/2 处
    // 提取支撑 [0, 2N-1] 对应的数值
    // 经过 cascade_iters 次迭代，间距 = 2^{-cascade_iters}
    // 最终目标间距 = 2^{-level}
    int downsample = 1 << (cascade_iters - level);
    int n_out = support_len * (1 << level) + 1;

    std::vector<double> result(n_out, 0.0);

    // 寻找 phi 的起始位置：在级联输出中，phi 的中心对应级联起始的 delta 位置
    // 经过 convergent 迭代，delta 扩散到覆盖整个支撑
    // phi 的支撑 [0, 2N-1] 在级联输出中对应某个偏移范围
    // 使用 phi(0) ≈ 0 和 phi(2N-1) ≈ 0 来确定偏移量

    // 计算级联输出的总支撑对应的物理区间
    // 查找最大值位置以确定 phi 的峰值位置
    double max_val = 0.0;
    int max_idx = 0;
    for (int i = 0; i < (int)cur.size(); ++i) {
        if (std::abs(cur[i]) > max_val) {
            max_val = std::abs(cur[i]);
            max_idx = i;
        }
    }

    // phi 的峰值大约在 N - 0.5 附近（对于 Daubechies N）
    // 级联中峰值在 max_idx，对应物理位置约为 max_idx / 2^{cascade_iters}
    // phi 的支撑起始物理位置 = 峰值位置 - (峰值在支撑内的偏移)
    // 峰值在支撑内的偏移 ≈ N - 0.5（实际略有不同）
    // 支撑起始物理位置 ≈ max_idx/2^{cascade_iters} - (N - 0.5)

    // 简化处理：利用 phi(0)=0，找到从 max_idx 向左的第一个零交叉
    // 但这对于光滑的 phi 不太可靠

    // 实际上级联从 delta(0) 出发，经过 iter 次迭代后
    // φ(x) 近似位于 x ∈ [shift, shift + 2N-1]
    // 其中 shift 由滤波器系数决定
    // 对于 Daubechies 正交小波，shift = 0（精确地）
    // 所以级联输出的位置 0 对应于 φ(0)

    // 级联中第 k 个采样点对应物理位置 k / 2^{cascade_iters}
    // 要获取 φ(i / 2^{level})，需取级联输出中的索引 i * 2^{cascade_iters - level}
    // 即 i * downsample

    for (int i = 0; i < n_out; ++i) {
        int src_idx = i * downsample;
        if (src_idx >= 0 && src_idx < (int)cur.size())
            result[i] = cur[src_idx];
    }

    return result;
}

std::vector<double> cascade_scaling_function_d2(int N, int level) {
    // ============================================================
    // 使用级联算法 + 四阶中心差分计算 phi''(x)
    // ============================================================

    auto phi = cascade_scaling_function(N, level);
    int n = (int)phi.size();
    std::vector<double> d2phi(n, 0.0);

    double h_step = 1.0 / (1 << level);  // 网格间距

    // 四阶中心差分系数
    // f''(x) ≈ (-f(x+2h) + 16f(x+h) - 30f(x) + 16f(x-h) - f(x-2h)) / (12h²)
    double denom = 12.0 * h_step * h_step;

    for (int i = 0; i < n; ++i) {
        double val = 0.0;

        // f(x+2h)
        int i2 = i + 2;
        if (i2 >= 0 && i2 < n) val -= phi[i2];

        // f(x+h)
        int i1 = i + 1;
        if (i1 >= 0 && i1 < n) val += 16.0 * phi[i1];

        // f(x)
        val -= 30.0 * phi[i];

        // f(x-h)
        int im1 = i - 1;
        if (im1 >= 0 && im1 < n) val += 16.0 * phi[im1];

        // f(x-2h)
        int im2 = i - 2;
        if (im2 >= 0 && im2 < n) val -= phi[im2];

        d2phi[i] = val / denom;
    }

    return d2phi;
}
