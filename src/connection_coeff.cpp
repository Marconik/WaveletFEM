#include "connection_coeff.h"
#include "daubechies.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

ConnectionTable::ConnectionTable(int N, int cascade_level)
    : N(N)
    , cascade_level(cascade_level)
    , support_len(2 * N - 1)
{
    // 预计算 phi 和 phi'' 在细网格上的值
    phi = cascade_scaling_function(N, cascade_level);
    d2phi = cascade_scaling_function_d2(N, cascade_level);

    n_samples = (int)phi.size();
    dx = 1.0 / (1 << cascade_level);
}

double ConnectionTable::eval_phi(double x) const {
    if (x < 0 || x > support_len) return 0.0;

    double idx_f = x / dx;
    int i = (int)idx_f;
    double frac = idx_f - i;

    if (i < 0) return phi[0];
    if (i >= n_samples - 1) return phi[n_samples - 1];

    // 线性插值
    return (1.0 - frac) * phi[i] + frac * phi[i + 1];
}

double ConnectionTable::eval_d2phi(double x) const {
    if (x < 0 || x > support_len) return 0.0;

    double idx_f = x / dx;
    int i = (int)idx_f;
    double frac = idx_f - i;

    if (i < 0) return d2phi[0];
    if (i >= n_samples - 1) return d2phi[n_samples - 1];

    return (1.0 - frac) * d2phi[i] + frac * d2phi[i + 1];
}

double ConnectionTable::integrate_simpson(
    const std::vector<double>& f, int i_start, int i_end) const {

    if (i_start >= i_end) return 0.0;

    int n_intervals = i_end - i_start;
    if (n_intervals < 1) return 0.0;

    double sum = 0.0;

    // 确保区间数为偶数以使用 Simpson 1/3 法则
    if (n_intervals % 2 == 1) {
        // 用梯形法则处理最后一个区间
        double trap = (f[i_end - 1] + f[i_end]) * dx / 2.0;
        sum += trap;
        i_end--;
        n_intervals--;
    }

    // Simpson 1/3 法则: h/3 * (f0 + 4f1 + 2f2 + 4f3 + ... + fn)
    double simpson = f[i_start] + f[i_end];
    for (int i = i_start + 1; i < i_end; ++i) {
        simpson += (i - i_start) % 2 == 1 ? 4.0 * f[i] : 2.0 * f[i];
    }
    sum += simpson * dx / 3.0;

    return sum;
}

double ConnectionTable::eval_C(int k, int n, double x) const {
    // C_k^n(x) = ∫_0^x phi^{(n)}(y - k) * phi(y) dy
    // 支撑集条件: phi(y) ≠ 0 only if y ∈ [0, 2N-1]
    //             phi^{(n)}(y-k) ≠ 0 only if y-k ∈ [0, 2N-1]
    //                              i.e., y ∈ [k, k+2N-1]

    if (x <= 0) return 0.0;

    // 确定积分区间
    double a = std::max(0.0, (double)k);
    double b = std::min(x, (double)(k + support_len));
    // phi(y) 也在 [0, 2N-1] 上非零
    b = std::min(b, (double)support_len);

    if (a >= b) return 0.0;

    // 在级联网格上做数值积分
    int i_a = (int)std::ceil(a / dx);
    int i_b = (int)std::floor(b / dx);
    if (i_a > i_b) return 0.0;

    // 构建被积函数在网格点上的值
    int n_pts = i_b - i_a + 1;
    std::vector<double> integrand(n_pts);

    const std::vector<double>* kernel_ptr;
    if (n == 0) kernel_ptr = &phi;
    else if (n == 2) kernel_ptr = &d2phi;
    else throw std::invalid_argument("Only n=0 (mass) and n=2 (stiffness) are supported");

    const auto& kernel = *kernel_ptr;

    for (int i = i_a; i <= i_b; ++i) {
        double y = i * dx;
        double phi_y = phi[i];  // phi(y)

        // phi^{(n)}(y - k): 在网格点上求值
        double arg = y - k;
        double kernel_val = 0.0;
        if (arg >= 0 && arg <= support_len) {
            int j = (int)std::round(arg / dx);
            if (j >= 0 && j < (int)kernel.size()) {
                kernel_val = kernel[j];
            }
        }

        integrand[i - i_a] = phi_y * kernel_val;
    }

    return integrate_simpson(integrand, 0, n_pts - 1);
}

double ConnectionTable::eval_M(int k, int r, double x) const {
    // M_k^r(x) = ∫_0^x y^r * phi(y - k) dy

    if (x <= 0) return 0.0;

    double a = std::max(0.0, (double)k);
    double b = std::min(x, (double)(k + support_len));

    if (a >= b) return 0.0;

    int i_a = (int)std::ceil(a / dx);
    int i_b = (int)std::floor(b / dx);
    if (i_a > i_b) return 0.0;

    int n_pts = i_b - i_a + 1;
    std::vector<double> integrand(n_pts);

    for (int i = i_a; i <= i_b; ++i) {
        double y = i * dx;

        // phi(y - k)
        double arg = y - k;
        double phi_val = 0.0;
        if (arg >= 0 && arg <= support_len) {
            int j = (int)std::round(arg / dx);
            if (j >= 0 && j < (int)phi.size()) {
                phi_val = phi[j];
            }
        }

        integrand[i - i_a] = std::pow(y, r) * phi_val;
    }

    return integrate_simpson(integrand, 0, n_pts - 1);
}
