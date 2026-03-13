#include "fluid.h"
#include <chrono>
#include <iomanip>
#include <numeric>
#include <filesystem>

namespace fs = std::filesystem;

// ==================== 计时工具 ====================
using Clock     = std::chrono::high_resolution_clock;
using TimePoint = std::chrono::time_point<Clock>;

inline TimePoint now() { return Clock::now(); }
inline double elapsed_ms(TimePoint s, TimePoint e) {
    return std::chrono::duration<double, std::milli>(e - s).count();
}

// ==================== 性能报告 ====================
void printPerfReport(int rank, int num_procs,
                     double wall_ms,
                     double compute_ms,   // 本进程纯计算累计
                     double io_ms,        // 本进程I/O累计
                     const std::string& step_info,
                     const Mesh& local_mesh)
{
    std::vector<double> all_compute(num_procs);
    std::vector<int>    all_cells(num_procs);
    int local_cells = local_mesh.nx * local_mesh.ny;

    MPI_Gather(&compute_ms,   1, MPI_DOUBLE, all_compute.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&local_cells,  1, MPI_INT,    all_cells.data(),   1, MPI_INT,    0, MPI_COMM_WORLD);

    double global_io_max = 0.0;
    MPI_Reduce(&io_ms, &global_io_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(rank != 0) return;

    double avg_step_ms = (all_compute.size() > 0) ? wall_ms / std::max(1, (int)all_compute.size()) : 0.0;

    std::cout << "\n╔══════════════════════════════════════╗\n"
              <<   "║           性能统计报告               ║\n"
              <<   "╚══════════════════════════════════════╝\n"
              << std::fixed << std::setprecision(3);
    std::cout << "  总运行时间   : " << std::setw(10) << wall_ms       << " ms"
              << "  (" << wall_ms / 1000.0 << " s)\n";
    std::cout << "  每步平均耗时 : " << std::setw(10) << avg_step_ms   << " ms/step\n";
    std::cout << "  最大I/O耗时  : " << std::setw(10) << global_io_max << " ms"
              << "  (占总时间 "
              << std::setprecision(1) << global_io_max / wall_ms * 100.0 << "%)\n";

    double max_t = *std::max_element(all_compute.begin(), all_compute.end());
    double min_t = *std::min_element(all_compute.begin(), all_compute.end());
    double avg_t = std::accumulate(all_compute.begin(), all_compute.end(), 0.0) / num_procs;
    double imbalance = (avg_t > 0) ? max_t / avg_t : 1.0;

    std::cout << "\n  ── 负载均衡 ──────────────────────────\n"
              << std::setprecision(3)
              << std::setw(6)  << "进程"
              << std::setw(14) << "计算时间(ms)"
              << std::setw(10) << "单元数"
              << std::setw(8)  << "占比\n"
              << "  " << std::string(36, '-') << "\n";

    for(int i = 0; i < num_procs; ++i) {
        double pct = (max_t > 0) ? all_compute[i] / max_t * 100.0 : 0.0;
        int bar_len = (int)(pct / 10.0);
        std::cout << std::setw(6)  << i
                  << std::setw(14) << all_compute[i]
                  << std::setw(10) << all_cells[i]
                  << "  [" << std::string(bar_len, '|') << std::string(10 - bar_len, ' ') << "] "
                  << std::setprecision(1) << pct << "%\n";
    }

    std::cout << "  " << std::string(36, '-') << "\n";
    std::cout << std::setprecision(3);
    std::cout << "  最大: " << max_t << " ms  最小: " << min_t
              << " ms  均值: " << avg_t << " ms\n"
              << "  不均衡系数: " << imbalance
              << (imbalance < 1.05 ? "  ✓ 优秀"
                 : imbalance < 1.20 ? "  ○ 良好"
                                     : "  ✗ 需优化") << "\n"
              << "  并行效率:   " << std::setprecision(1)
              << 100.0 / imbalance << "%\n"
              << "  时间信息:   " << step_info << "\n"
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

    if(argc != 4) {
        if(rank == 0) {
            std::cerr << "用法: mpirun -np <进程数> ./solver <网格> <dt/auto/time> <步数或结束时间>\n";
            std::cerr << "示例1: mpirun -np 4 ./solver ./mesh 0.0001 1000  # 指定 dt+步数\n";
            std::cerr << "示例2: mpirun -np 4 ./solver ./mesh auto 2000    # 自动 CFL+步数\n";
            std::cerr << "示例3: mpirun -np 4 ./solver ./mesh time 0.1     # 指定 dt 到 t_end=0.1\n";
        }
        MPI_Finalize();
        return 1;
    }

    std::string mesh_folder = argv[1];
    std::string dt_mode     = argv[2];
    std::string step_arg    = argv[3];

    bool auto_dt   = (dt_mode == "auto");
    bool time_mode = (dt_mode == "time");
    double dt_user = 0.0;
    double t_end   = 0.0;
    int timesteps  = 0;

    if(time_mode) {
        t_end   = std::stod(step_arg);
        dt_user = 0.0; // 初始 dt 可用 CFL 决定
    } else if(!auto_dt) {
        dt_user  = std::stod(dt_mode);
        timesteps = std::stoi(step_arg);
    } else {
        timesteps = std::stoi(step_arg);
    }

    if(rank == 0) {
        std::cout << "================================\n"
                  << "二维欧拉方程 WENO 求解器 - MPI 并行\n"
                  << "================================\n"
                  << "网格文件夹: " << mesh_folder << "\n"
                  << "模式:      " << dt_mode << "\n"
                  << "参数:      " << step_arg << "\n"
                  << "MPI进程数: " << num_procs << "\n"
                  << "================================\n";
    }

    // ── 加载 & 分割网格 ──
    if(rank == 0) std::cout << "正在加载网格...\n";
    Mesh full_mesh(mesh_folder);
    std::vector<Mesh> sub_meshes = splitMeshVertically(full_mesh, num_procs);
    Mesh local_mesh = sub_meshes[rank];

    // ── 初始 CFL 和 dt ──
    double global_umax = computeMaxSpeed(local_mesh, local_mesh.gamma);
    MPI_Allreduce(MPI_IN_PLACE, &global_umax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    double dt = auto_dt || time_mode ? global_umax * local_mesh.da : dt_user;
    double physical_time = 0.0;

    // ── 性能计数器 ──
    double compute_ms = 0.0;
    double io_ms      = 0.0;
    const int output_interval = 100;
    int step = 0;

    while(true) {

        if(time_mode && physical_time >= t_end) break;
        if(!time_mode && step >= timesteps) break;

        // 重算 CFL
        if(auto_dt || time_mode) {
            double lmax_local  = computeMaxSpeed(local_mesh, local_mesh.gamma);
            double lmax_global = 0.0;
            MPI_Allreduce(&lmax_local, &lmax_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            dt = lmax_global * local_mesh.da; // 自动 CFL
            if(time_mode && physical_time + dt > t_end) dt = t_end - physical_time;
        }

        // 计算
        auto t0 = now();
        updateMesh(local_mesh, dt, rank, num_procs);
        compute_ms += elapsed_ms(t0, now());
        physical_time += dt;
        step++;

        // 输出进度
        if(rank == 0 && step % std::max(1, timesteps/1000) == 0) {
            double used = elapsed_ms(wall_start, now());
            double eta  = time_mode ? (t_end - physical_time)/dt * used/step
                                    : used / step * (timesteps - step);
            std::cout << std::fixed << std::setprecision(2)
                      << "Step " << step
                      << "  t=" << std::setprecision(6) << physical_time
                      << "  已用: " << used/1000.0 << " s"
                      << "  预计剩余: " << eta/1000.0 << " s\n";
        }

        // 保存数据
       if ((step ) % output_interval == 0|| (step ) == timesteps|| (time_mode && physical_time >= t_end)) {
            std::string folder = "result/" + std::to_string(step);
            if(rank == 0) fs::create_directories(folder);
            MPI_Barrier(MPI_COMM_WORLD);

            auto t_io = now();
            saveMeshData(local_mesh, rank, folder);
            io_ms += elapsed_ms(t_io, now());
        }
    }

    double wall_ms = elapsed_ms(wall_start, now());
    MPI_Barrier(MPI_COMM_WORLD);

    std::string step_info = time_mode ? "t_end=" + std::to_string(t_end)
                                      : "steps=" + std::to_string(timesteps);
    printPerfReport(rank, num_procs, wall_ms, compute_ms, io_ms, step_info, local_mesh);

    MPI_Finalize();
    return 0;
}