#include "fluid.h"

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
                     int    timesteps,
                     const Mesh& local_mesh)
{
    // ---------- 负载均衡：收集各进程计算时间 ----------
    std::vector<double> all_compute(num_procs);
    std::vector<int>    all_cells(num_procs);
    int local_cells = local_mesh.nx * local_mesh.ny;

    MPI_Gather(&compute_ms,   1, MPI_DOUBLE, all_compute.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&local_cells,  1, MPI_INT,    all_cells.data(),   1, MPI_INT,    0, MPI_COMM_WORLD);

    // ---------- I/O：取各进程最大值作为代表 ----------
    double global_io_max = 0.0;
    MPI_Reduce(&io_ms, &global_io_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(rank != 0) return;

    // ── 基础时间统计 ──
    double avg_step_ms = (timesteps > 0) ? wall_ms / timesteps : 0.0;

    std::cout << "\n╔══════════════════════════════════════╗" << std::endl;
    std::cout <<   "║           性能统计报告               ║" << std::endl;
    std::cout <<   "╚══════════════════════════════════════╝" << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  总运行时间   : " << std::setw(10) << wall_ms       << " ms"
              << "  (" << wall_ms / 1000.0 << " s)\n";
    std::cout << "  每步平均耗时 : " << std::setw(10) << avg_step_ms   << " ms/step\n";
    std::cout << "  最大I/O耗时  : " << std::setw(10) << global_io_max << " ms"
              << "  (占总时间 "
              << std::setprecision(1) << global_io_max / wall_ms * 100.0 << "%)\n";

    // ── 负载均衡 ──
    double max_t = *std::max_element(all_compute.begin(), all_compute.end());
    double min_t = *std::min_element(all_compute.begin(), all_compute.end());
    double avg_t = std::accumulate(all_compute.begin(), all_compute.end(), 0.0) / num_procs;
    double imbalance = (avg_t > 0) ? max_t / avg_t : 1.0;

    std::cout << "\n  ── 负载均衡 ──────────────────────────\n";
    std::cout << std::setprecision(3);
    std::cout << std::setw(6)  << "进程"
              << std::setw(14) << "计算时间(ms)"
              << std::setw(10) << "单元数"
              << std::setw(8)  << "占比\n";
    std::cout << "  " << std::string(36, '-') << "\n";

    for(int i = 0; i < num_procs; ++i) {
        double pct = (max_t > 0) ? all_compute[i] / max_t * 100.0 : 0.0;
        // 简单ASCII柱状条（宽度10）
        int bar_len = (int)(pct / 10.0);
        std::cout << std::setw(6)  << i
                  << std::setw(14) << all_compute[i]
                  << std::setw(10) << all_cells[i]
                  << "  [" << std::string(bar_len, '|')
                           << std::string(10 - bar_len, ' ') << "] "
                  << std::setprecision(1) << pct << "%\n";
    }

    std::cout << "  " << std::string(36, '-') << "\n";
    std::cout << std::setprecision(3);
    std::cout << "  最大: " << max_t << " ms  最小: " << min_t
              << " ms  均值: " << avg_t << " ms\n";
    std::cout << "  不均衡系数: " << std::setprecision(3) << imbalance
              << (imbalance < 1.05 ? "  ✓ 优秀"
                : imbalance < 1.20 ? "  ○ 良好"
                                   : "  ✗ 需优化") << "\n";
    std::cout << "  并行效率:   " << std::setprecision(1)
              << 100.0 / imbalance << "%\n";
    std::cout << "════════════════════════════════════════\n";
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
            std::cerr << "用法: mpirun -np <进程数> ./solver <网格文件夹> <dt> <时间步数>" << std::endl;
            std::cerr << "例如: mpirun -np 4 ./solver ./mesh_data 0.0001 1000" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    std::string mesh_folder = argv[1];
    double dt               = std::stod(argv[2]);
    int    timesteps        = std::stoi(argv[3]);

    if(rank == 0) {
        std::cout << "================================\n"
                  << "二维欧拉方程求解器 - MPI并行版本\n"
                  << "================================\n"
                  << "网格文件夹: " << mesh_folder << "\n"
                  << "时间步长:   " << dt          << "\n"
                  << "时间步数:   " << timesteps   << "\n"
                  << "MPI进程数:  " << num_procs   << "\n"
                  << "================================\n";
    }

    // ── 加载 & 分割网格 ──
    if(rank == 0) std::cout << "正在加载网格...\n";
    Mesh full_mesh(mesh_folder);
    if(rank == 0) {
        std::cout << "网格大小: " << full_mesh.ny << " x " << full_mesh.nx << "\n"
                  << "Gamma值:  " << full_mesh.gamma << "\n";
    }

    std::vector<Mesh> sub_meshes = splitMeshVertically(full_mesh, num_procs);
    Mesh local_mesh = sub_meshes[rank];

    // ── CFL 检查 ──
    double global_umax = computeMaxSpeed(local_mesh, local_mesh.gamma);
    MPI_Allreduce(MPI_IN_PLACE, &global_umax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    double cfl = global_umax * dt / local_mesh.da;
    if(rank == 0) {
        std::cout << "CFL数 = " << cfl << " (应 ≤ 0.5)\n";
        if(cfl > 0.5) std::cerr << "警告: CFL数过大，可能不稳定！\n";
        std::cout << "\n开始时间步进...\n";
    }

    // ── 性能计数器（仅三项）──
    double compute_ms = 0.0;
    double io_ms      = 0.0;
    const int print_interval  = std::max(1, timesteps / 10);
    const int output_interval = 100;

    for(int step = 0; step < timesteps; ++step) {

        // 计算计时
        auto t0 = now();
        updateMesh(local_mesh, dt, rank, num_procs);
        compute_ms += elapsed_ms(t0, now());

        // 进度输出
        if(rank == 0 && (step + 1) % print_interval == 0) {
            double used = elapsed_ms(wall_start, now());
            double eta  = used / (step + 1) * (timesteps - step - 1);
            std::cout << std::fixed << std::setprecision(2)
                      << std::setw(6) << (step + 1) << " / " << timesteps
                      << "  已用: " << used / 1000.0 << " s"
                      << "  预计剩余: " << eta / 1000.0 << " s\n";
        }

        // 定期输出
        if((step + 1) % output_interval == 0) {
            recoverPrimitives(local_mesh);
            std::string folder = "result/" + std::to_string(step + 1);
            if(rank == 0) fs::create_directories(folder);
            MPI_Barrier(MPI_COMM_WORLD);

            // I/O 计时
            auto t_io = now();
            saveMeshData(local_mesh, rank, folder);
            io_ms += elapsed_ms(t_io, now());
        }
    }

    double wall_ms = elapsed_ms(wall_start, now());
    MPI_Barrier(MPI_COMM_WORLD);

    printPerfReport(rank, num_procs, wall_ms, compute_ms, io_ms, timesteps, local_mesh);

    MPI_Finalize();
    return 0;
}