#include "wavelet_fem.h"
#include "daubechies.h"
#include <Eigen/SparseLU>
#include <cmath>
#include <iostream>
#include <algorithm>

// ============================================================
// 构造函数
// ============================================================

WaveletFEMSolver::WaveletFEMSolver(const ProblemConfig& config)
    : config_(config)
{
    // 预计算连接系数（级联级别 14 以获得高精度）
    conn_ = std::make_unique<ConnectionTable>(config.N, 14);

    // 确定基函数指标范围
    determine_index_range();

    // 组装一维矩阵
    assemble_1d_matrices();
}

// ============================================================
// 设置算例函数
// ============================================================

void WaveletFEMSolver::set_initial_condition(ScalarFunction2D s) { s_ = s; }
void WaveletFEMSolver::set_source_term(ScalarFunction3D f) { f_ = f; }
void WaveletFEMSolver::set_boundary_functions(std::array<ScalarFunction3D, 4> g) { g_ = g; }

void WaveletFEMSolver::set_exact_solution(ScalarFunction3D u_exact) {
    u_exact_ = u_exact;
    has_exact_ = true;
}

// ============================================================
// 向量/矩阵转换
// ============================================================

VectorXd WaveletFEMSolver::mat_to_vec(const MatrixXd& M) {
    // 列主序拉直: vec(M)[i + j*n] = M(i,j) (Eigen 默认列主序)
    VectorXd v(M.size());
    int idx = 0;
    for (int j = 0; j < M.cols(); ++j)
        for (int i = 0; i < M.rows(); ++i)
            v(idx++) = M(i, j);
    return v;
}

MatrixXd WaveletFEMSolver::vec_to_mat(const VectorXd& v, int n) {
    MatrixXd M(n, n);
    int idx = 0;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            M(i, j) = v(idx++);
    return M;
}

// ============================================================
// 指标范围
// ============================================================

void WaveletFEMSolver::determine_index_range() {
    // I = { -N+2, -N+3, ..., 2^{j+1} - 1 }
    k_min_ = -config_.N + 2;
    int k_max = (1 << (config_.j + 1)) - 1;
    n_ = k_max - k_min_ + 1;
}

// ============================================================
// 组装一维矩阵
// ============================================================

void WaveletFEMSolver::assemble_1d_matrices() {
    // 使用公式:
    // ã_{l,k} = C_{k-l}^0(2^{j+1} - l) - C_{k-l}^0(-l)
    // c̃_{l,k} = C_{k-l}^2(2^{j+1} - l) - C_{k-l}^2(-l)
    //
    // 其中 l, k ∈ I, 且 C_k^n(x) 已对无量纲变量定义
    // 注意：在 section 2 的公式中，积分区间是 [0,2] 和
    // 缩放后的基函数，需要相应的变量转换

    A_ = MatrixXd::Zero(n_, n_);
    C_ = MatrixXd::Zero(n_, n_);

    double L = config_.bx - config_.ax;  // 区间长度 = 2
    int pow2_j = 1 << config_.j;          // 2^j
    double ymax = pow2_j * L;             // 缩放后 y 区间右端 = 2^{j+1}
    int support_len = 2 * config_.N - 1;  // φ 支撑整数长度
    int Dmax = 2 * config_.N - 2;         // |k-l| 的最大值
    double scale = std::pow(2.0, 2 * config_.j);

    // LRT 精确二阶连接系数 Γ_d = ∫ φ'(x)φ'(x-d) dx,  d ∈ [-Dmax, Dmax]
    // 内部刚度： ∫ φ''(y-k)φ(y-l) dy = -Γ_{k-l}（分部积分），再乘 2^{2j}
    std::vector<double> Gamma = connection_coeff_d1d1(config_.N);

    for (int l_idx = 0; l_idx < n_; ++l_idx) {
        int l = k_min_ + l_idx;
        for (int k_idx = 0; k_idx < n_; ++k_idx) {
            int k = k_min_ + k_idx;
            int d = k - l;

            // 支撑不重叠 → 元素为 0
            if (std::abs(d) > Dmax) {
                A_(l_idx, k_idx) = 0.0;
                C_(l_idx, k_idx) = 0.0;
                continue;
            }

            // 判断两个基函数的重叠区是否完全落在区间 [0, ymax] 内
            //   φ(y-k) 支撑 [k, k+support_len]，φ(y-l) 支撑 [l, l+support_len]
            //   重叠区 = [max(k,l), min(k,l)+support_len]
            int hi = std::max(k, l);
            int lo = std::min(k, l);
            bool interior = (hi >= 0) && (lo + support_len <= ymax);

            if (interior) {
                // 内部：质量矩阵由正交归一性精确给出 = δ；刚度由 LRT 精确给出
                A_(l_idx, k_idx) = (d == 0) ? 1.0 : 0.0;
                C_(l_idx, k_idx) = -scale * Gamma[d + Dmax];
            } else {
                // 边界函数：支撑被区间截断，暂用原数值积分（截断连接系数）
                double upper = ymax - l;
                double lower = -l;
                A_(l_idx, k_idx) = conn_->eval_C(d, 0, upper)
                                 - conn_->eval_C(d, 0, lower);
                C_(l_idx, k_idx) = scale * (
                    conn_->eval_C(d, 2, upper) - conn_->eval_C(d, 2, lower)
                );
            }
        }
    }

    // 确保对称性（数值误差可能导致微小不对称）
    A_ = 0.5 * (A_ + A_.transpose());
    C_ = 0.5 * (C_ + C_.transpose());
}

// ============================================================
// Kronecker 积系统
// ============================================================

void WaveletFEMSolver::build_kronecker_system(double dt,
                                               SparseMatrix& L_sp,
                                               SparseMatrix& R_sp) const {
    // Ā = A ⊗ A (n² × n²)
    // C̄ = a₁(C ⊗ A) + a₂(A ⊗ C)
    // L = Ā - Δt/2 · C̄
    // R = Ā + Δt/2 · C̄

    int n2 = n_ * n_;
    double half_dt = dt / 2.0;

    // 使用 triplet 列表构建稀疏矩阵
    std::vector<Triplet> L_trip, R_trip;
    L_trip.reserve(n_ * n_ * n_ * n_ / 4);  // 粗略估计
    R_trip.reserve(n_ * n_ * n_ * n_ / 4);

    for (int i = 0; i < n_; ++i) {
        for (int j = 0; j < n_; ++j) {
            double A_ij = A_(i, j);
            double C_ij = C_(i, j);

            if (std::abs(A_ij) < 1e-16 && std::abs(C_ij) < 1e-16) continue;

            for (int p = 0; p < n_; ++p) {
                for (int q = 0; q < n_; ++q) {
                    double A_pq = A_(p, q);
                    double C_pq = C_(p, q);

                    if (std::abs(A_pq) < 1e-16 && std::abs(C_pq) < 1e-16) continue;

                    // Kronecker 积索引: (i,j)⊗(p,q) → 全局行 i*n+p, 列 j*n+q
                    int row = i * n_ + p;
                    int col = j * n_ + q;

                    double A_bar = A_ij * A_pq;
                    double C_bar = config_.a1 * C_ij * A_pq
                                 + config_.a2 * A_ij * C_pq;

                    double L_val = A_bar - half_dt * C_bar;
                    double R_val = A_bar + half_dt * C_bar;

                    if (std::abs(L_val) > 1e-16)
                        L_trip.emplace_back(row, col, L_val);
                    if (std::abs(R_val) > 1e-16)
                        R_trip.emplace_back(row, col, R_val);
                }
            }
        }
    }

    L_sp.resize(n2, n2);
    R_sp.resize(n2, n2);
    L_sp.setFromTriplets(L_trip.begin(), L_trip.end());
    R_sp.setFromTriplets(R_trip.begin(), R_trip.end());

    L_sp.makeCompressed();
    R_sp.makeCompressed();
}

// ============================================================
// 函数投影
// ============================================================

MatrixXd WaveletFEMSolver::project_2d(ScalarFunction2D func) const {
    // 投影二维函数到基函数空间
    // S_{l,m} = ∫∫ s(x,y) φ_{j,l}(x) φ_{j,m}(y) dx dy
    //
    // 使用二维 Gauss-Legendre 求积在细化网格上计算

    MatrixXd S = MatrixXd::Zero(n_, n_);

    double Lx = config_.bx - config_.ax;
    double Ly = config_.by - config_.ay;

    // 积分网格：在区间上均匀取点
    int nq = 200;  // 每方向求积点数
    double hx = Lx / nq;
    double hy = Ly / nq;
    double dA = hx * hy;

    // 预计算所有求积点上 phi_{j,k}(x) 的值
    int pow2_j = 1 << config_.j;
    double sqrt_pow2_j = std::sqrt((double)pow2_j);  // 2^{j/2}

    // 对每个求积点，计算 phi 值
    std::vector<std::vector<double>> phi_x(n_, std::vector<double>(nq + 1));
    std::vector<std::vector<double>> phi_y(n_, std::vector<double>(nq + 1));

    for (int idx = 0; idx < n_; ++idx) {
        int k = k_min_ + idx;
        for (int iq = 0; iq <= nq; ++iq) {
            double x = config_.ax + iq * hx;
            double arg = pow2_j * x - k;
            phi_x[idx][iq] = sqrt_pow2_j * conn_->eval_phi(arg);
        }
        for (int iq = 0; iq <= nq; ++iq) {
            double y = config_.ay + iq * hy;
            double arg = pow2_j * y - k;
            phi_y[idx][iq] = sqrt_pow2_j * conn_->eval_phi(arg);
        }
    }

    // 二维求积
    for (int ix = 0; ix < nq; ++ix) {
        double x_mid = config_.ax + (ix + 0.5) * hx;
        for (int iy = 0; iy < nq; ++iy) {
            double y_mid = config_.ay + (iy + 0.5) * hy;
            double f_val = func(x_mid, y_mid);

            if (std::abs(f_val) < 1e-16) continue;

            for (int l = 0; l < n_; ++l) {
                double px = 0.5 * (phi_x[l][ix] + phi_x[l][ix + 1]);
                for (int m = 0; m < n_; ++m) {
                    double py = 0.5 * (phi_y[m][iy] + phi_y[m][iy + 1]);
                    S(l, m) += f_val * px * py * dA;
                }
            }
        }
    }

    return S;
}

MatrixXd WaveletFEMSolver::project_3d(ScalarFunction3D func, double t) const {
    // 将 func(x,y,t) 投影为 n×n 矩阵
    auto func_t = [&func, t](double x, double y) -> double {
        return func(x, y, t);
    };
    return project_2d(func_t);
}

// ============================================================
// 施加边界条件
// ============================================================

void WaveletFEMSolver::apply_boundary_conditions(
    SparseMatrix& L, VectorXd& rhs, double t) const {

    int n2 = n_ * n_;
    double pow2_j = 1 << config_.j;
    double sqrt_pow2_j = std::sqrt(pow2_j);
    double Lx = config_.bx - config_.ax;
    double Ly = config_.by - config_.ay;

    // ------------------------------------------------------------
    // 基函数在边界点的求值（对每个基函数指标 a 都计算）
    //   φ_a(x) = 2^{j/2} φ(2^j x - (k_min + a))
    //   φ_a'(x) = d/dx，用物理坐标中心差分（自动含 2^j 因子）
    // ------------------------------------------------------------
    auto phi_at = [&](int a, double xx) -> double {
        return sqrt_pow2_j * conn_->eval_phi(pow2_j * xx - (k_min_ + a));
    };
    auto dphi_at = [&](int a, double xx) -> double {
        double ex = 1e-5;
        return (phi_at(a, xx + ex) - phi_at(a, xx - ex)) / (2.0 * ex);
    };

    // 边界行系数因子 B_a = ±b φ_a'(x_b) + d φ_a(x_b)
    // 投影后的边界方程：  Σ_{a,b} U_{a,b} B_a A_{b,m} = g_m
    // （x 边界对 x 指标 a 作用、对 y 用质量矩阵投影；y 边界对称处理）

    // ---- 边 0: x = ax (左边界)  -b ∂u/∂x + d u = g_0 ----
    {
        double b_bc = config_.b[0], d_bc = config_.d[0];
        double x0 = config_.ax;
        std::vector<double> B(n_);
        for (int a = 0; a < n_; ++a)
            B[a] = -b_bc * dphi_at(a, x0) + d_bc * phi_at(a, x0);

        for (int m = 0; m < n_; ++m) {
            int row = 0 * n_ + m;          // x 指标 = 0，y 指标 = m
            int k_m = k_min_ + m;
            for (int col = 0; col < n2; ++col) L.coeffRef(row, col) = 0.0;

            double g_int = 0.0;
            int nq = 100; double hy = Ly / nq;
            for (int iy = 0; iy < nq; ++iy) {
                double y_mid = config_.ay + (iy + 0.5) * hy;
                double g_val = g_[0](config_.ax, y_mid, t);
                g_int += g_val * sqrt_pow2_j * conn_->eval_phi(pow2_j * y_mid - k_m) * hy;
            }

            for (int a = 0; a < n_; ++a) {
                if (std::abs(B[a]) < 1e-15) continue;
                for (int b = 0; b < n_; ++b) {
                    double coeff = B[a] * A_(b, m);
                    if (std::abs(coeff) < 1e-15) continue;
                    L.coeffRef(row, a * n_ + b) = coeff;
                }
            }
            rhs(row) = g_int;
        }
    }

    // ---- 边 1: x = bx (右边界)  +b ∂u/∂x + d u = g_1 ----
    {
        double b_bc = config_.b[1], d_bc = config_.d[1];
        double x1 = config_.bx;
        std::vector<double> B(n_);
        for (int a = 0; a < n_; ++a)
            B[a] = b_bc * dphi_at(a, x1) + d_bc * phi_at(a, x1);

        for (int m = 0; m < n_; ++m) {
            int row = (n_ - 1) * n_ + m;   // x 指标 = n-1，y 指标 = m
            int k_m = k_min_ + m;
            for (int col = 0; col < n2; ++col) L.coeffRef(row, col) = 0.0;

            double g_int = 0.0;
            int nq = 100; double hy = Ly / nq;
            for (int iy = 0; iy < nq; ++iy) {
                double y_mid = config_.ay + (iy + 0.5) * hy;
                double g_val = g_[1](config_.bx, y_mid, t);
                g_int += g_val * sqrt_pow2_j * conn_->eval_phi(pow2_j * y_mid - k_m) * hy;
            }

            for (int a = 0; a < n_; ++a) {
                if (std::abs(B[a]) < 1e-15) continue;
                for (int b = 0; b < n_; ++b) {
                    double coeff = B[a] * A_(b, m);
                    if (std::abs(coeff) < 1e-15) continue;
                    L.coeffRef(row, a * n_ + b) = coeff;
                }
            }
            rhs(row) = g_int;
        }
    }

    // ---- 边 2: y = ay (下边界)  -b ∂u/∂y + d u = g_2 ----
    {
        double b_bc = config_.b[2], d_bc = config_.d[2];
        double y0 = config_.ay;
        std::vector<double> B(n_);   // 依赖 y 指标 b
        for (int b = 0; b < n_; ++b)
            B[b] = -b_bc * dphi_at(b, y0) + d_bc * phi_at(b, y0);

        for (int l = 0; l < n_; ++l) {
            int row = l * n_ + 0;          // x 指标 = l，y 指标 = 0
            int k_l = k_min_ + l;
            for (int col = 0; col < n2; ++col) L.coeffRef(row, col) = 0.0;

            double g_int = 0.0;
            int nq = 100; double hx = Lx / nq;
            for (int ix = 0; ix < nq; ++ix) {
                double x_mid = config_.ax + (ix + 0.5) * hx;
                double g_val = g_[2](x_mid, config_.ay, t);
                g_int += g_val * sqrt_pow2_j * conn_->eval_phi(pow2_j * x_mid - k_l) * hx;
            }

            for (int a = 0; a < n_; ++a) {       // x 指标
                double Aal = A_(a, l);
                if (std::abs(Aal) < 1e-15) continue;
                for (int b = 0; b < n_; ++b) {   // y 指标
                    double coeff = B[b] * Aal;
                    if (std::abs(coeff) < 1e-15) continue;
                    L.coeffRef(row, a * n_ + b) = coeff;
                }
            }
            rhs(row) = g_int;
        }
    }

    // ---- 边 3: y = by (上边界)  +b ∂u/∂y + d u = g_3 ----
    {
        double b_bc = config_.b[3], d_bc = config_.d[3];
        double y1 = config_.by;
        std::vector<double> B(n_);
        for (int b = 0; b < n_; ++b)
            B[b] = b_bc * dphi_at(b, y1) + d_bc * phi_at(b, y1);

        for (int l = 0; l < n_; ++l) {
            int row = l * n_ + (n_ - 1);   // x 指标 = l，y 指标 = n-1
            int k_l = k_min_ + l;
            for (int col = 0; col < n2; ++col) L.coeffRef(row, col) = 0.0;

            double g_int = 0.0;
            int nq = 100; double hx = Lx / nq;
            for (int ix = 0; ix < nq; ++ix) {
                double x_mid = config_.ax + (ix + 0.5) * hx;
                double g_val = g_[3](x_mid, config_.by, t);
                g_int += g_val * sqrt_pow2_j * conn_->eval_phi(pow2_j * x_mid - k_l) * hx;
            }

            for (int a = 0; a < n_; ++a) {
                double Aal = A_(a, l);
                if (std::abs(Aal) < 1e-15) continue;
                for (int b = 0; b < n_; ++b) {
                    double coeff = B[b] * Aal;
                    if (std::abs(coeff) < 1e-15) continue;
                    L.coeffRef(row, a * n_ + b) = coeff;
                }
            }
            rhs(row) = g_int;
        }
    }

    L.prune(0.0);  // 清除显式零元
}

// ============================================================
// Crank-Nicolson 时间步
// ============================================================

VectorXd WaveletFEMSolver::time_step(
    const SparseMatrix& L, const SparseMatrix& R,
    const VectorXd& U_curr_vec, const VectorXd& F_curr_vec,
    const VectorXd& F_next_vec) const {

    // RHS = R · U^n + Δt/2 (F^{n+1} + F^n)
    // The F vectors are the pure Galerkin projections of f(x,y,t)
    // The Δt/2 factor is applied here

    // NOTE: This function receives dt implicitly through L and R
    // which already contain the Δt factor. The F terms need Δt/2 scaling.
    // However, this function is not currently used from solve();
    // the solve() method applies the scaling directly.

    VectorXd rhs = R * U_curr_vec + 0.5 * (F_next_vec + F_curr_vec);

    // 求解 L · U_next = rhs
    Eigen::SparseLU<SparseMatrix> solver;
    solver.compute(L);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("SparseLU decomposition failed");
    }

    VectorXd U_next_vec = solver.solve(rhs);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("SparseLU solve failed");
    }

    return U_next_vec;
}

// ============================================================
// 重构近似解
// ============================================================

double WaveletFEMSolver::reconstruct(double x, double y, int step) const {
    if (step < 0 || step >= (int)history_.size()) return 0.0;

    const MatrixXd& U = history_[step].U;
    double pow2_j = 1 << config_.j;
    double sqrt_pow2_j = std::sqrt(pow2_j);

    double sum = 0.0;
    for (int l = 0; l < n_; ++l) {
        int k_l = k_min_ + l;
        double arg_x = pow2_j * x - k_l;
        double phi_x = sqrt_pow2_j * conn_->eval_phi(arg_x);

        for (int m = 0; m < n_; ++m) {
            int k_m = k_min_ + m;
            double arg_y = pow2_j * y - k_m;
            double phi_y = sqrt_pow2_j * conn_->eval_phi(arg_y);

            sum += U(l, m) * phi_x * phi_y;
        }
    }

    return sum;
}

MatrixXd WaveletFEMSolver::reconstruct_on_grid(int step) const {
    int nv = config_.nviz;
    MatrixXd grid(nv, nv);

    double hx = (config_.bx - config_.ax) / (nv - 1);
    double hy = (config_.by - config_.ay) / (nv - 1);

    for (int i = 0; i < nv; ++i) {
        double x = config_.ax + i * hx;
        for (int j = 0; j < nv; ++j) {
            double y = config_.ay + j * hy;
            grid(i, j) = reconstruct(x, y, step);
        }
    }

    return grid;
}

// ============================================================
// L² 误差
// ============================================================

double WaveletFEMSolver::compute_l2_error(const MatrixXd& U, double t) const {
    if (!has_exact_) return -1.0;

    double Lx = config_.bx - config_.ax;
    double Ly = config_.by - config_.ay;
    double dx = Lx / 150;
    double dy = Ly / 150;
    double dA = dx * dy;

    double error_sq = 0.0;
    double norm_sq = 0.0;

    double pow2_j = 1 << config_.j;
    double sqrt_pow2_j = std::sqrt(pow2_j);

    for (int ix = 0; ix < 150; ++ix) {
        double x = config_.ax + (ix + 0.5) * dx;
        for (int iy = 0; iy < 150; ++iy) {
            double y = config_.ay + (iy + 0.5) * dy;

            // 重构近似解
            double u_approx = 0.0;
            for (int l = 0; l < n_; ++l) {
                int k_l = k_min_ + l;
                double arg_x = pow2_j * x - k_l;
                double phi_x = sqrt_pow2_j * conn_->eval_phi(arg_x);

                for (int m = 0; m < n_; ++m) {
                    int k_m = k_min_ + m;
                    double arg_y = pow2_j * y - k_m;
                    double phi_y = sqrt_pow2_j * conn_->eval_phi(arg_y);
                    u_approx += U(l, m) * phi_x * phi_y;
                }
            }

            double u_ex = u_exact_(x, y, t);
            double diff = u_approx - u_ex;
            error_sq += diff * diff * dA;
            norm_sq += u_ex * u_ex * dA;
        }
    }

    if (norm_sq < 1e-16) return -1.0;
    return std::sqrt(error_sq / norm_sq);
}

// ============================================================
// 主求解流程
// ============================================================

std::vector<TimeStepRecord> WaveletFEMSolver::solve() {
    history_.clear();

    std::cout << "=== WaveletFEM Solver ===" << std::endl;
    std::cout << "Daubechies N = " << config_.N
              << ", resolution j = " << config_.j
              << ", basis count n = " << n_ << std::endl;
    std::cout << "System size: " << n_ * n_ << " × " << n_ * n_ << std::endl;

    // ---- 初始条件投影 ----
    std::cout << "Projecting initial condition..." << std::endl;
    MatrixXd S_mat = project_2d(s_);

    // 解 Ā · Ū⁰ = vec(S)
    SparseMatrix A_bar, dummy;
    build_kronecker_system(0.0, A_bar, dummy);  // Δt=0 时 L=R=Ā
    VectorXd S_vec = mat_to_vec(S_mat);

    Eigen::SparseLU<SparseMatrix> solver_init;
    solver_init.compute(A_bar);
    VectorXd U0_vec = solver_init.solve(S_vec);
    MatrixXd U0 = vec_to_mat(U0_vec, n_);

    double t = 0.0;
    double dt = config_.dt0;

    TimeStepRecord rec0;
    rec0.t = t;
    rec0.dt = dt;
    rec0.U = U0;
    rec0.l2_error = compute_l2_error(U0, t);
    history_.push_back(rec0);

    std::cout << "t = " << t << ", dt = " << dt
              << (has_exact_ ? ", L2 error = " + std::to_string(rec0.l2_error) : "")
              << std::endl;

    // ---- 时间推进 ----
    int step = 0;

    while (t < config_.T - 1e-12) {
        // 调整步长不超过剩余时间
        if (t + dt > config_.T) {
            dt = config_.T - t;
        }

        // 构造系统矩阵
        SparseMatrix L, R;
        build_kronecker_system(dt, L, R);

        // 计算非齐次项投影
        MatrixXd F_curr_mat = project_3d(f_, t);
        MatrixXd F_next_mat = project_3d(f_, t + dt);

        VectorXd F_curr_vec = mat_to_vec(F_curr_mat);
        VectorXd F_next_vec = mat_to_vec(F_next_mat);

        // 构建完整 RHS = R · U^n + Δt/2 (F^{n+1} + F^n)
        VectorXd U_curr_vec = mat_to_vec(history_.back().U);
        VectorXd RHS = R * U_curr_vec + 0.5 * dt * (F_next_vec + F_curr_vec);

        // 施加边界条件：修改 L 的边界行和 RHS 的边界行
        apply_boundary_conditions(L, RHS, t + dt);

        // 求解
        Eigen::SparseLU<SparseMatrix> solver_lu;
        solver_lu.compute(L);

        if (solver_lu.info() != Eigen::Success) {
            throw std::runtime_error("SparseLU decomposition failed at t="
                                     + std::to_string(t));
        }

        VectorXd U_next_vec = solver_lu.solve(RHS);

        if (solver_lu.info() != Eigen::Success) {
            throw std::runtime_error("SparseLU solve failed at t="
                                     + std::to_string(t));
        }

        MatrixXd U_next = vec_to_mat(U_next_vec, n_);
        t += dt;

        TimeStepRecord rec;
        rec.t = t;
        rec.dt = dt;
        rec.U = U_next;
        rec.l2_error = compute_l2_error(U_next, t);
        history_.push_back(rec);

        std::cout << "t = " << t << ", dt = " << dt
                  << (has_exact_ ? ", L2 error = " + std::to_string(rec.l2_error) : "")
                  << std::endl;

        // ---- 变步长调整 ----
        if (config_.variable_dt) {
            // 单调不减步长策略：每 10 步增大一次
            if ((step + 1) % 10 == 0) {
                double dt_new = dt * 1.5;
                dt = std::min(dt_new, config_.T / 20.0);  // 上限
            }
        }

        ++step;

        // 安全检查
        if (step > 100000) {
            std::cerr << "Warning: too many time steps, stopping." << std::endl;
            break;
        }
    }

    std::cout << "Solved in " << history_.size() << " time steps." << std::endl;
    return history_;
}
