#include "fluid.h"

// ==================== 计时工具 ====================
using Clock     = std::chrono::high_resolution_clock;
using TimePoint = std::chrono::time_point<Clock>;

inline TimePoint now() { return Clock::now(); }
inline double elapsed_ms(TimePoint s, TimePoint e) {
    return std::chrono::duration<double, std::milli>(e - s).count();
}

// ==================== DG CFL 安全上限 ====================
//
// SSP-RK3 + DG P 阶的理论稳定条件（Cockburn-Shu 1989）：
//
//   Δt ≤ C / (2P+1) · h / λ_max
//
//   其中 C 对 SSP-RK3 约为 1.0（1D 分析值）
//
// 推广到 2D 张量积网格（x/y 方向独立施加 CFL）：
//   每个方向的谱半径均为 λ_max，两方向叠加后将有效 CFL 减半：
//
//   Δt_safe = 1 / (2·(2P+1)) · h / λ_max
//
//   P=1 → Δt_safe ≈ h / (6  · λ_max)    CFL_lim ≈ 0.167
//   P=2 → Δt_safe ≈ h / (10 · λ_max)    CFL_lim ≈ 0.100
//   P=3 → Δt_safe ≈ h / (14 · λ_max)    CFL_lim ≈ 0.071
//
// 注：此处 CFL 定义为  ν = λ_max · Δt / h（无量纲），
//     与 FV 的 CFL ≤ 0.5 **不可直接比较**。
// =============================================================================
static constexpr double DG_CFL_SAFETY = 0.9;   // 留 10% 安全余量

inline double dgCFLLimit() {
    return DG_CFL_SAFETY / (2.0 * (2 * DG_P + 1));
}

/// 由当前场的最大波速计算 DG 安全时间步
inline double dgSafeDt(double lambda_max, double h) {
    const double eps = 1e-7;
    lambda_max = std::max(lambda_max, eps);
    return dgCFLLimit() * h / lambda_max;
}

// =============================================================================
// [修复 1] computeMaxSpeedDG
//
// 原版 computeMaxSpeed 只遍历单元平均值（模式 0），低估了 DG 多项式在界面处
// 的实际波速。DG 数值通量 HLLC 使用的是界面处的重构值：
//
//   U(ξ=±1, η=±1) = Σ_m û_m · φ_m(±1, ±1)
//
// 对 P=2，L_2(±1) = 1，界面处修正量 ≈ û_{2,0} + û_{0,2}，量级不可忽略。
//
// 此函数对每个内部单元的 4 个角点求波速，取全局最大值。
// 角点覆盖了所有 4 条界面的两端，已足以估计最坏情况。
// =============================================================================
static double computeMaxSpeedDG(const Mesh& mesh)
{
    const double gam = mesh.gamma;
    constexpr double eps = 1e-12;

    // 4 个界面角点的参考坐标
    const double XI [4] = {-1.0, +1.0, -1.0, +1.0};
    const double ETA[4] = {-1.0, -1.0, +1.0, +1.0};

    double umax = 0.0;

    for (int i = 0; i < mesh.ny; ++i) {
        for (int j = 0; j < mesh.nx; ++j) {
            if (mesh.bctype(i, j) != 0) continue;

            for (int k = 0; k < 4; ++k) {
                double xi = XI[k], eta = ETA[k];

                // DG 多项式重构：U(xi,eta) = Σ_m û_m φ_m(xi,eta)
                double U[4] = {0.0, 0.0, 0.0, 0.0};
                for (int m = 0; m < DG_NM; ++m) {
                    double b = phi2(m, xi, eta);
                    U[0] += mesh.dof[0][m](i, j) * b;
                    U[1] += mesh.dof[1][m](i, j) * b;
                    U[2] += mesh.dof[2][m](i, j) * b;
                    U[3] += mesh.dof[3][m](i, j) * b;
                }

                double rho = std::max(U[0], eps);
                double u   = U[1] / rho;
                double v   = U[2] / rho;
                double p   = std::max((gam - 1.0) * (U[3] - 0.5 * rho * (u*u + v*v)), eps);
                double a   = std::sqrt(gam * p / rho);
                double spd = std::sqrt(u*u + v*v) + a;

                umax = std::max(umax, spd);
            }
        }
    }

    return umax;
}

// ==================== 性能报告 ====================
void printPerfReport(int rank, int num_procs,
                     double wall_ms,
                     double compute_ms,
                     double io_ms,
                     int    timesteps,
                     const Mesh& local_mesh)
{
    std::vector<double> all_compute(num_procs);
    std::vector<int>    all_cells(num_procs);
    int local_cells = local_mesh.nx * local_mesh.ny;

    MPI_Gather(&compute_ms,  1, MPI_DOUBLE, all_compute.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&local_cells, 1, MPI_INT,    all_cells.data(),   1, MPI_INT,    0, MPI_COMM_WORLD);

    double global_io_max = 0.0;
    MPI_Reduce(&io_ms, &global_io_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank != 0) return;

    double avg_step_ms = (timesteps > 0) ? wall_ms / timesteps : 0.0;

    std::cout << "\n╔══════════════════════════════════════╗\n"
              <<   "║           性能统计报告               ║\n"
              <<   "╚══════════════════════════════════════╝\n"
              << std::fixed << std::setprecision(3);
    std::cout << "  总运行时间   : " << std::setw(10) << wall_ms       << " ms"
              << "  (" << wall_ms / 1000.0 << " s)\n"
              << "  每步平均耗时 : " << std::setw(10) << avg_step_ms   << " ms/step\n"
              << "  最大I/O耗时  : " << std::setw(10) << global_io_max << " ms"
              << "  (占总时间 "
              << std::setprecision(1) << global_io_max / wall_ms * 100.0 << "%)\n";

    double max_t     = *std::max_element(all_compute.begin(), all_compute.end());
    double min_t     = *std::min_element(all_compute.begin(), all_compute.end());
    double avg_t     = std::accumulate(all_compute.begin(), all_compute.end(), 0.0) / num_procs;
    double imbalance = (avg_t > 0) ? max_t / avg_t : 1.0;

    std::cout << "\n  ── 负载均衡 ──────────────────────────\n"
              << std::setprecision(3)
              << std::setw(6)  << "进程"
              << std::setw(14) << "计算时间(ms)"
              << std::setw(10) << "单元数"
              << std::setw(8)  << "占比\n"
              << "  " << std::string(36, '-') << "\n";

    for (int i = 0; i < num_procs; ++i) {
        double pct     = (max_t > 0) ? all_compute[i] / max_t * 100.0 : 0.0;
        int    bar_len = (int)(pct / 10.0);
        std::cout << std::setw(6)  << i
                  << std::setw(14) << all_compute[i]
                  << std::setw(10) << all_cells[i]
                  << "  [" << std::string(bar_len, '|')
                           << std::string(10 - bar_len, ' ') << "] "
                  << std::setprecision(1) << pct << "%\n";
    }

    std::cout << "  " << std::string(36, '-') << "\n"
              << std::setprecision(3)
              << "  最大: " << max_t << " ms  最小: " << min_t
              << " ms  均值: " << avg_t << " ms\n"
              << "  不均衡系数: " << std::setprecision(3) << imbalance << "\n"
              << "  并行效率:   " << std::setprecision(1) << 100.0 / imbalance << "%\n"
              << "════════════════════════════════════════\n";
}

// ==================== 主函数 ====================
int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    auto wall_start = now();

    // ── 参数解析 ──────────────────────────────────────────────────────────
    // 用法 A（指定 dt）：  solver <网格> <dt>  <步数>
    // 用法 B（自动 CFL）： solver <网格> auto  <步数>
    // ─────────────────────────────────────────────────────────────────────
    if (argc != 4) {
        if (rank == 0) {
            std::cerr << "用法 A（指定dt）  : mpirun -np N ./solver_DG <网格> <dt>   <步数>\n"
                      << "用法 B（自动CFL） : mpirun -np N ./solver_DG <网格> auto  <步数>\n"
                      << "例如              : mpirun -np 4 ./solver_DG ./mesh auto  2000\n";
        }
        MPI_Finalize();
        return 1;
    }

    std::string mesh_folder = argv[1];
    std::string dt_arg      = argv[2];
    int         timesteps   = std::stoi(argv[3]);
    bool        auto_dt     = (dt_arg == "auto");
    double      dt_user     = auto_dt ? 0.0 : std::stod(dt_arg);

    // ── 加载 & 分割网格 ───────────────────────────────────────────────────
    if (rank == 0) std::cout << "正在加载网格...\n";
    Mesh full_mesh(mesh_folder);
    if (rank == 0) {
        std::cout << "网格大小 : " << full_mesh.ny << " x " << full_mesh.nx << "\n"
                  << "Gamma    : " << full_mesh.gamma << "\n"
                  << "DG 阶数  : P = " << DG_P
                  << "  模式数 = " << DG_NM << " / 格\n";
    }

    std::vector<Mesh> sub_meshes = splitMeshVertically(full_mesh, num_procs);
    Mesh local_mesh = sub_meshes[rank];

    // ── CFL 检查与 dt 确定 ────────────────────────────────────────────────
    //
    // [修复 1] 使用 computeMaxSpeedDG 代替 computeMaxSpeed：
    //   后者只读单元平均值（模式 0），会低估 DG 界面处的实际波速。
    //   computeMaxSpeedDG 对每格 4 个角点做多项式重构后求波速，
    //   正确反映 HLLC 数值通量实际遇到的最大特征速度。
    // ─────────────────────────────────────────────────────────────────────
    double local_umax  = computeMaxSpeedDG(local_mesh);
    double global_umax = 0.0;
    MPI_Allreduce(&local_umax, &global_umax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    const double h         = local_mesh.da;
    const double cfl_limit = dgCFLLimit();
    const double dt_safe   = dgSafeDt(global_umax, h);
    double physical_time = 0.0;
    double dt;
    if (auto_dt) {
        dt = dt_safe;
        if (rank == 0)
            std::cout << "自动 dt 模式：dt = " << dt
                      << "  (CFL_lim=" << cfl_limit << ")\n";
    } else {
        dt = dt_user;
    }

    double cfl_actual = global_umax * dt / h;

    if (rank == 0) {
        std::cout << "================================\n"
                  << "二维欧拉方程求解器 - DG 并行版本\n"
                  << "================================\n"
                  << "网格文件夹 : " << mesh_folder  << "\n"
                  << "时间步长   : " << dt           << "\n"
                  << "时间步数   : " << timesteps    << "\n"
                  << "MPI进程数  : " << num_procs    << "\n"
                  << "格子尺寸 h : " << h            << "\n"
                  << "最大波速(DG界面重构) : " << global_umax << "\n"
                  << "--------------------------------\n"
                  << "DG CFL 信息（P=" << DG_P << ", 2D, SSP-RK3）\n"
                  << "  安全上限  ν_lim = 1/(2·(2P+1)) · safety\n"
                  << "           = 1/(2·" << (2*DG_P+1) << ") · "
                  << DG_CFL_SAFETY << " ≈ " << cfl_limit << "\n"
                  << "  安全 dt          = " << dt_safe   << "\n"
                  << "  实际 ν = λΔt/h   = " << cfl_actual << "\n";

        if (cfl_actual > cfl_limit) {
            std::cerr << "\n[错误] CFL 数 " << cfl_actual
                      << " 超过 DG P=" << DG_P
                      << " 稳定上限 " << cfl_limit << "\n"
                      << "       建议 dt ≤ " << dt_safe << "\n"
                      << "       或使用 'auto' 自动设定 dt\n";
        } else if (cfl_actual > cfl_limit * 0.8) {
            std::cerr << "\n[警告] CFL 数 " << cfl_actual
                      << " 接近 DG 稳定上限 " << cfl_limit
                      << "，建议适当减小 dt\n";
        } else {
            std::cout << "  CFL 检查通过 ✓\n";
        }
        std::cout << "================================\n\n开始时间步进...\n";
    }

    // [说明] cfl_actual 由 global_umax（MPI_Allreduce）和相同的 dt/h 计算得出，
    // 各进程值完全一致，无需 Bcast。直接用本地值判断即可。
    if (cfl_actual > cfl_limit) {
        MPI_Finalize();
        return 1;
    }

    // ── 主时间步循环 ──────────────────────────────────────────────────────
    double compute_ms = 0.0;
    double io_ms      = 0.0;
    const int print_interval        = std::max(1, timesteps / 1000);
    const int output_interval       = 100;
    const int cfl_recheck_interval  = 50;

    for (int step = 0; step < timesteps; ++step)
    {
        // ── 定期重新计算最大波速（无论 auto_dt 与否均执行）─────────────
        //
        // [修复 2] 原代码只在 auto_dt 模式下重算波速，固定 dt 模式完全
        // 不检查，流场激波增强后可能无声发散。
        //
        // 修复策略：
        //   - auto_dt 模式：重算后更新 dt（行为不变）
        //   - 固定 dt 模式：仅检查并输出警告，dt 不变
        // ─────────────────────────────────────────────────────────────────
        if (step > 0 && step % cfl_recheck_interval == 0) {
            // [修复 1] 动态重算同样使用 computeMaxSpeedDG
            double lmax_local  = computeMaxSpeedDG(local_mesh);
            double lmax_global = 0.0;
            MPI_Allreduce(&lmax_local, &lmax_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

            if (auto_dt) {
                dt = dgSafeDt(lmax_global, h);
            } else {
                // [修复 2] 固定 dt 模式：检查 CFL 是否超限，超限时输出警告
                double cfl_now = lmax_global * dt / h;
                if (rank == 0 && cfl_now > cfl_limit) {
                    std::cerr << std::fixed << std::setprecision(6)
                              << "[警告] 第 " << step << " 步：CFL = " << cfl_now
                              << " 超过稳定上限 " << cfl_limit
                              << "  (λ_max=" << lmax_global
                              << ", dt=" << dt << ")\n"
                              << "         固定dt模式无法自动调整，"
                              << "建议减小dt或改用 'auto' 模式\n";
                }
            }
        }

        // ── 计算 ─────────────────────────────────────────────────────────
        auto t0 = now();
        updateMesh(local_mesh, dt, rank, num_procs);
        compute_ms += elapsed_ms(t0, now());
        physical_time += dt;

        // ── 进度输出 ──────────────────────────────────────────────────────
        if (rank == 0 && (step + 1) % print_interval == 0) {
            double used = elapsed_ms(wall_start, now());
            double eta  = used / (step + 1) * (timesteps - step - 1);
            std::cout << std::fixed << std::setprecision(2)
                      << std::setw(6) << (step + 1) << " / " << timesteps                      
                      << "  dt=" << std::setprecision(12) << dt
                      << "  t="  << std::setprecision(6) << physical_time

                      << "  已用: " << std::setprecision(2) << used / 1000.0 << " s"
                      << "  预计剩余: " << eta / 1000.0 << " s\n";
        }

        // ── 定期保存 ──────────────────────────────────────────────────────
        if ((step + 1) % output_interval == 0) {
            std::string folder = "result/" + std::to_string(step + 1);
            if (rank == 0) fs::create_directories(folder);
            MPI_Barrier(MPI_COMM_WORLD);

            auto t_io = now();
            saveMeshData(local_mesh, rank, folder);
            io_ms += elapsed_ms(t_io, now());
        }
    }

    double wall_ms = elapsed_ms(wall_start, now());
    MPI_Barrier(MPI_COMM_WORLD);

    printPerfReport(rank, num_procs, wall_ms, compute_ms, io_ms,
                    timesteps, local_mesh);

    MPI_Finalize();
    return 0;
}