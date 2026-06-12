#include "types.h"
#include "wavelet_fem.h"
#include "examples.h"
#include "gnuplot_output.h"
#include <iostream>
#include <string>
#include <cstring>
#include <cmath>

void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " [options]\n"
              << "Options:\n"
              << "  --example N     Example number (1 or 2, default: 1)\n"
              << "  --j LEVEL       Resolution level (default: 4)\n"
              << "  --N ORDER       Daubechies order (default: 6)\n"
              << "  --dt DT         Initial time step (default: 0.005)\n"
              << "  --T TIME        Final time (default: 1.0)\n"
              << "  --variable-dt   Use variable time step (default: on)\n"
              << "  --fixed-dt      Use fixed time step\n"
              << "  --help          Show this help\n";
}

int main(int argc, char* argv[]) {
    // 默认参数
    int example_num = 1;
    int j = 4;
    int N = 10;
    double dt = 0.001;
    double T = 1.0;
    bool variable_dt = true;

    // 解析命令行参数
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "--example" && i + 1 < argc) {
            example_num = std::stoi(argv[++i]);
        } else if (arg == "--j" && i + 1 < argc) {
            j = std::stoi(argv[++i]);
        } else if (arg == "--N" && i + 1 < argc) {
            N = std::stoi(argv[++i]);
        } else if (arg == "--dt" && i + 1 < argc) {
            dt = std::stod(argv[++i]);
        } else if (arg == "--T" && i + 1 < argc) {
            T = std::stod(argv[++i]);
        } else if (arg == "--variable-dt") {
            variable_dt = true;
        } else if (arg == "--fixed-dt") {
            variable_dt = false;
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    // 验证参数
    if (example_num < 1 || example_num > 2) {
        std::cerr << "Error: example must be 1 or 2\n";
        return 1;
    }
    if (N < 2 || N > 10) {
        std::cerr << "Error: N must be between 2 and 10\n";
        return 1;
    }
    if (j < 1 || j > 7) {
        std::cerr << "Error: j must be between 1 and 7\n";
        return 1;
    }

    std::cout << "========================================\n";
    std::cout << "  WaveletFEM - Wavelet-Galerkin FEM\n";
    std::cout << "  2D Heat Equation Solver\n";
    std::cout << "========================================\n\n";

    // 获取算例配置
    ProblemConfig config;
    ScalarFunction2D initial;
    ScalarFunction3D source, exact;
    std::array<ScalarFunction3D, 4> g;

    if (example_num == 1) {
        std::cout << "Example 1: Robin BC\n";
        std::cout << "  u(x,y,t) = exp(-2x) * sin(pi*y/2) * (t+1)/(t+2)\n\n";
        config = Example1::config();
        initial = Example1::initial;
        source = Example1::source;
        exact = Example1::exact;
        g = {Example1::g_left, Example1::g_right,
             Example1::g_bottom, Example1::g_top};
    } else {
        std::cout << "Example 2: Neumann BC\n";
        std::cout << "  u(x,y,t) = exp(-((x+1)+(y+1))/(2t+1))\n\n";
        config = Example2::config();
        initial = Example2::initial;
        source = Example2::source;
        exact = Example2::exact;
        g = {Example2::g_left, Example2::g_right,
             Example2::g_bottom, Example2::g_top};
    }

    // 覆盖用户参数
    config.j = j;
    config.N = N;
    config.dt0 = dt;
    config.T = T;
    config.variable_dt = variable_dt;

    std::cout << "Configuration:\n";
    std::cout << "  Domain: [" << config.ax << ", " << config.bx << "] × ["
              << config.ay << ", " << config.by << "]\n";
    std::cout << "  Daubechies N = " << config.N << ", resolution j = "
              << config.j << "\n";
    std::cout << "  Time: T = " << config.T << ", dt0 = " << config.dt0
              << (config.variable_dt ? " (variable)" : " (fixed)") << "\n";
    std::cout << "  BC: b = [" << config.b[0] << "," << config.b[1] << ","
              << config.b[2] << "," << config.b[3] << "]"
              << ", d = [" << config.d[0] << "," << config.d[1] << ","
              << config.d[2] << "," << config.d[3] << "]\n\n";

    try {
        // 创建求解器
        std::cout << "Initializing solver..." << std::endl;
        WaveletFEMSolver solver(config);

        // 设置算例函数
        solver.set_initial_condition(initial);
        solver.set_source_term(source);
        solver.set_boundary_functions(g);
        solver.set_exact_solution(exact);

        std::cout << "Solver initialized. Basis count: "
                  << solver.num_basis() << " per direction.\n\n";

        // 求解
        auto history = solver.solve();

        // 输出最终误差
        if (!history.empty() && history.back().l2_error >= 0) {
            std::cout << "\nFinal L2 relative error: "
                      << history.back().l2_error << "\n";
        }

        // 输出结果文件
        std::cout << "\nExporting results..." << std::endl;
        export_all_steps("results", solver, history, config);

        std::cout << "\nDone! Results saved in results/ directory.\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
