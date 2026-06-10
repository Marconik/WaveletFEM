#ifndef WAVELET_FEM_SOLVER_H
#define WAVELET_FEM_SOLVER_H

#include "types.h"
#include "connection_coeff.h"
#include <memory>

// ============================================================
// 小波-Galerkin 有限元求解器
// ============================================================

class WaveletFEMSolver {
public:
    /// 构造函数，接受问题配置
    explicit WaveletFEMSolver(const ProblemConfig& config);

    /// 设置算例函数
    void set_initial_condition(ScalarFunction2D s);
    void set_source_term(ScalarFunction3D f);
    void set_boundary_functions(std::array<ScalarFunction3D, 4> g);

    /// 设置精确解（用于误差评估，可选）
    void set_exact_solution(ScalarFunction3D u_exact);

    /// 执行求解，返回各时间步的记录
    std::vector<TimeStepRecord> solve();

    /// 重构近似解 u(x, y, t_n) 在给定点上的值
    double reconstruct(double x, double y, int step) const;

    // ---- 查询接口 ----

    int num_basis() const { return n_; }
    int nviz() const { return config_.nviz; }
    int k_min() const { return k_min_; }
    int resolution_j() const { return config_.j; }
    double ax() const { return config_.ax; }
    double bx() const { return config_.bx; }
    double ay() const { return config_.ay; }
    double by() const { return config_.by; }

    /// 在可视化网格上重构解（供 GNUPlot 输出使用）
    MatrixXd reconstruct_on_grid(int step) const;

    /// 获取连接系数表（外部使用）
    const ConnectionTable& conn_table() const { return *conn_; }

private:
    ProblemConfig config_;

    // 基函数指标范围
    int n_;                        // 每方向的基函数数量
    int k_min_;                    // 最小平移指标

    // 一维矩阵
    MatrixXd A_;                   // 质量矩阵 (n × n, 带状对称)
    MatrixXd C_;                   // 刚度矩阵 (n × n, 带状对称)

    // 连接系数
    std::unique_ptr<ConnectionTable> conn_;

    // 算例函数
    ScalarFunction2D s_;           // 初始条件 u(x,y,0)
    ScalarFunction3D f_;           // 非齐次源项 f(x,y,t)
    std::array<ScalarFunction3D, 4> g_; // 边界函数 g_i(y,t) 或 g_i(x,t)
    ScalarFunction3D u_exact_;    // 精确解（可选）

    // 边界条件标记
    bool has_exact_ = false;

    // 求解历史
    std::vector<TimeStepRecord> history_;

    // ---- 内部方法 ----

    /// 确定基函数指标范围
    void determine_index_range();

    /// 组装一维质量矩阵 A 和刚度矩阵 C
    void assemble_1d_matrices();

    /// 构建 Kronecker 积系统矩阵
    /// L = Ā - Δt/2 * C̄,  R = Ā + Δt/2 * C̄
    void build_kronecker_system(double dt,
                                SparseMatrix& L,
                                SparseMatrix& R) const;

    /// 将二维函数投影到基函数空间：返回 n×n 矩阵
    MatrixXd project_2d(ScalarFunction2D func) const;

    /// 将三维函数在时间 t 投影到基函数空间：返回 n×n 矩阵
    MatrixXd project_3d(ScalarFunction3D func, double t) const;

    /// 施加边界条件到系统矩阵和右端向量
    void apply_boundary_conditions(SparseMatrix& L,
                                   VectorXd& rhs,
                                   double t) const;

    /// 单个 Crank-Nicolson 时间步
    VectorXd time_step(const SparseMatrix& L,
                       const SparseMatrix& R,
                       const VectorXd& U_curr_vec,
                       const VectorXd& F_curr_vec,
                       const VectorXd& F_next_vec) const;

    /// 计算 L² 相对误差
    double compute_l2_error(const MatrixXd& U, double t) const;

    /// 将 n×n 矩阵拉直为 n² 向量（列主序）
    static VectorXd mat_to_vec(const MatrixXd& M);

    /// 将 n² 向量重整为 n×n 矩阵（列主序）
    static MatrixXd vec_to_mat(const VectorXd& v, int n);
};

#endif // WAVELET_FEM_SOLVER_H
