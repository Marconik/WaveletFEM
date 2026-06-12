#ifndef WAVELET_FEM_GNUPLOT_OUTPUT_H
#define WAVELET_FEM_GNUPLOT_OUTPUT_H

#include "types.h"
#include "wavelet_fem.h"
#include <string>
#include <vector>

// ============================================================
// GNUPlot 数据输出与 GIF 生成
// ============================================================

/// 将单步重构解写入 .dat 纯文本文件
/// 格式: 矩阵形式，空格分隔列，行对应 x，列对应 y
void write_timestep_dat(const std::string& filename,
                        const WaveletFEMSolver& solver,
                        const MatrixXd& U,
                        double t,
                        int nx = 100, int ny = 100);

/// 生成 GNUPlot 脚本并调用 gnuplot 生成动画 GIF
void generate_gif(const std::string& output_gif,
                  const std::vector<std::string>& dat_files,
                  const std::vector<double>& timestamps,
                  const std::string& title = "2D Heat Equation Evolution",
                  const ProblemConfig& config = ProblemConfig(),
                  int delay_ms = 10);

/// 将所有时间步写入 .dat 文件并生成 GIF
void export_all_steps(const std::string& output_dir,
                      const WaveletFEMSolver& solver,
                      const std::vector<TimeStepRecord>& history,
                      const ProblemConfig& config);

#endif // WAVELET_FEM_GNUPLOT_OUTPUT_H
