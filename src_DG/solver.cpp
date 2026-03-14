#include "fluid.h"

// ==================== 计时工具 ====================
using Clock     = std::chrono::high_resolution_clock;
using TimePoint = std::chrono::time_point<Clock>;

inline TimePoint now() { return Clock::now(); }
inline double elapsed_ms(TimePoint s, TimePoint e) {
    return std::chrono::duration<double, std::milli>(e - s).count();
}

// ==================== DG CFL 安全上限 ====================
static constexpr double DG_CFL_SAFETY = 0.9;

inline double dgCFLLimit() {
    return DG_CFL_SAFETY / (2.0 * (2 * DG_P + 1));
}

inline double dgSafeDt(double lambda_max, double h) {
    const double eps = 1e-7;
    lambda_max = std::max(lambda_max, eps);
    return dgCFLLimit() * h / lambda_max;
}

// ==================== 优化后的 DG 最大波速 ====================
static double computeMaxSpeedDG_omp(const Mesh& mesh)
{
    const double gam = mesh.gamma;
    constexpr double eps = 1e-12;
    const double XI [4] = {-1.0, +1.0, -1.0, +1.0};
    const double ETA[4] = {-1.0, -1.0, +1.0, +1.0};

    // 预计算 phi2(m, xi, eta)  -> DG_NM x 4
    static double PHI[DG_NM][4];
    for (int m = 0; m < DG_NM; ++m)
        for (int k = 0; k < 4; ++k)
            PHI[m][k] = phi2(m, XI[k], ETA[k]);

    double umax = 0.0;

    #pragma omp parallel for collapse(2) reduction(max:umax)
    for (int i = 0; i < mesh.ny; ++i) {
        for (int j = 0; j < mesh.nx; ++j) {
            if (mesh.bctype(i,j) != 0) continue;
            for (int k = 0; k < 4; ++k) {
                double U[4] = {0.0,0.0,0.0,0.0};
                for (int m = 0; m < DG_NM; ++m) {
                    double b = PHI[m][k];
                    U[0] += mesh.dof[0][m](i,j) * b;
                    U[1] += mesh.dof[1][m](i,j) * b;
                    U[2] += mesh.dof[2][m](i,j) * b;
                    U[3] += mesh.dof[3][m](i,j) * b;
                }
                double rho = std::max(U[0], eps);
                double u   = U[1]/rho;
                double v   = U[2]/rho;
                double p   = std::max((gam-1.0)*(U[3]-0.5*rho*(u*u+v*v)), eps);
                double a   = std::sqrt(gam*p/rho);
                umax = std::max(umax, std::sqrt(u*u+v*v) + a);
            }
        }
    }
    return umax;
}

// ==================== 自适应 dt 更新 ====================
inline double adaptiveDt(double& dt, double h,
                         double local_umax,
                         bool auto_dt, bool time_mode,
                         double physical_time=0.0, double t_end=0.0)
{
    double global_umax = 0.0;
    const double dt_min = 1e-6;  // 最小时间步限制
    MPI_Allreduce(&local_umax, &global_umax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    double dt_safe = dgSafeDt(global_umax, h);

    if (auto_dt || time_mode) {
        dt = dt_safe;
        if (time_mode && physical_time + dt > t_end)
            dt = t_end - physical_time; // 防止 overshoot
    }
    dt = std::max(dt, dt_min);
    return dt;
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
    MPI_Gather(&local_cells, 1, MPI_INT,   all_cells.data(),   1, MPI_INT,   0, MPI_COMM_WORLD);

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

    // ── 参数解析 ──────────────────────────────────────────
    if (argc != 4) {
        if (rank == 0) {
            std::cerr << "用法 A（指定dt）  : mpirun -np N ./solver_DG <网格> <dt>   <步数>\n"
                      << "用法 B（自动CFL） : mpirun -np N ./solver_DG <网格> auto  <步数>\n"
                      << "用法 C（到指定时间）: mpirun -np N ./solver_DG <网格> time  <终止时间>\n"
                      << "例如              : mpirun -np 4 ./solver_DG ./mesh auto  2000\n";
        }
        MPI_Finalize();
        return 1;
    }

    std::string mesh_folder = argv[1];
    std::string dt_arg      = argv[2];
    std::string third_arg   = argv[3];

    bool auto_dt = (dt_arg == "auto");
    bool time_mode = (dt_arg == "time");

    int timesteps = 0;
    double t_end = 0.0;
    double dt_user = 0.0;

    if (time_mode) {
        t_end = std::stod(third_arg);
    } else {
        timesteps = std::stoi(third_arg);
        if (!auto_dt) dt_user = std::stod(dt_arg);
    }

    // ── 加载 & 分割网格 ─────────────────────────────────────
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

    // ── CFL 检查 & dt 确定 ─────────────────────────────────
    double local_umax  = computeMaxSpeedDG_omp(local_mesh);
    double global_umax = 0.0;
    MPI_Allreduce(&local_umax, &global_umax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    const double h         = local_mesh.da;
    const double cfl_limit = dgCFLLimit();
    const double dt_safe   = dgSafeDt(global_umax, h);

    double physical_time = 0.0;
    double dt;
    if (auto_dt || time_mode) {
        dt = dt_safe;
        if (rank == 0) {
            std::cout << (auto_dt ? "自动 dt 模式" : "time 模式")
                      << "：初始 dt = " << dt
                      << "  (CFL_lim=" << cfl_limit << ")\n";
        }
    } else {
        dt = dt_user;
    }

    double cfl_actual = global_umax * dt / h;

    if (rank == 0) {
        std::cout << "================================\n"
                  << "二维欧拉方程求解器 - DG 并行版本\n"
                  << "================================\n"
                  << "网格文件夹 : " << mesh_folder  << "\n"
                  << "时间步长   : " << dt              << "\n"
                  << "MPI进程数  : " << num_procs       << "\n"
                  << "格子尺寸 h : " << h               << "\n"
                  << "最大波速(DG界面重构) : " << global_umax << "\n"
                  << "--------------------------------\n"
                  << "DG CFL 信息（P=" << DG_P << ", 2D, SSP-RK3）\n"
                  << "  安全上限  ν_lim = 1/(2·(2P+1)) · safety"
                  << " = " << cfl_limit << "\n"
                  << "  安全 dt      = " << dt_safe << "\n"
                  << "  实际 ν = λΔt/h = " << cfl_actual << "\n";

        if (cfl_actual > cfl_limit) {
            std::cerr << "\n[错误] CFL 数 " << cfl_actual
                      << " 超过 DG P=" << DG_P
                      << " 稳定上限 " << cfl_limit
                      << "\n建议 dt ≤ " << dt_safe
                      << " 或使用 'auto' 模式\n";
        } else {
            std::cout << "  CFL 检查通过 ✓\n";
        }
        std::cout << "================================\n\n开始时间步进...\n";
    }

    if (cfl_actual > cfl_limit) {
        MPI_Finalize();
        return 1;
    }

    // ── 主时间步循环 ─────────────────────────────────────────
    double compute_ms = 0.0;
    double io_ms      = 0.0;
    const int print_interval       = std::max(1, (time_mode ? 10 : timesteps / 10000));
    const int output_interval      = 100;
    const int cfl_recheck_interval = 10;

    int step = 0;
    while ((!time_mode && step < timesteps) || (time_mode && physical_time < t_end)) {

        // 重新计算最大波速并调整 dt
    if (step > 0 && step % cfl_recheck_interval == 0) {
    double lmax_local  = computeMaxSpeedDG_omp(local_mesh);

    if (auto_dt || time_mode) {
        dt = adaptiveDt(dt, h, lmax_local, auto_dt, time_mode, physical_time, t_end);
    } else {
        double cfl_now = lmax_local * dt / h;
        if (rank == 0 && cfl_now > cfl_limit) {
            std::cerr << std::fixed << std::setprecision(6)
                      << "[警告] 第 " << step << " 步：CFL = " << cfl_now
                      << " 超过稳定上限 " << cfl_limit
                      << "  (λ_max=" << lmax_local
                      << ", dt=" << dt << ")\n"
                      << "  固定dt模式无法自动调整，建议减小dt或改用 'auto' 模式\n";
        }
    }
    }

        // 计算
        auto t0 = now();
        updateMesh(local_mesh, dt, rank, num_procs);
        compute_ms += elapsed_ms(t0, now());
        physical_time += dt;

        // 进度输出
        if (rank == 0 && step % print_interval == 0) {
            double used = elapsed_ms(wall_start, now());
            double eta  = time_mode ? (used / (physical_time) * (t_end - physical_time))
                                    : (used / (step+1) * ((double)timesteps - step - 1));
            std::cout << std::fixed << std::setprecision(2)
                      << std::setw(6) << (step+1) << " / "
                      << (time_mode ? "-" : std::to_string(timesteps))
                      << "  dt=" << std::setprecision(12) << dt
                      << "  t="  << std::setprecision(6) << physical_time
                      << "  已用: " << std::setprecision(2) << used / 1000.0 << " s"
                      << "  预计剩余: " << eta / 1000.0 << " s\n";
        }

        // 定期保存
        if ((step + 1) % output_interval == 0|| (step + 1) == timesteps|| (time_mode && physical_time >= t_end)) {
            std::string folder = "result/" + std::to_string(step+1);
            if (rank == 0) fs::create_directories(folder);
            MPI_Barrier(MPI_COMM_WORLD);

            auto t_io = now();
            saveMeshData(local_mesh, rank, folder);
            io_ms += elapsed_ms(t_io, now());
        }

        ++step;
    }

    double wall_ms = elapsed_ms(wall_start, now());
    MPI_Barrier(MPI_COMM_WORLD);
    printPerfReport(rank, num_procs, wall_ms, compute_ms, io_ms,
                    (!time_mode ? timesteps : step), local_mesh);

    MPI_Finalize();
    return 0;
}