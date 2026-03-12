#include "fluid.h"
// ==================== 主函数 ====================
int main(int argc, char* argv[]) 
{    
    // 初始化MPI环境
    MPI_Init(&argc, &argv);
    
    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    
    // ========================================
    // 参数检查和解析
    // ========================================
    if(argc != 4) {
        if(rank == 0) {
            std::cerr << "用法: mpirun -np <进程数> ./solver <网格文件夹> <dt> <时间步数>" << std::endl;
            std::cerr << "例如: mpirun -np 4 ./solver ./mesh_data 0.0001 1000" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }
    
    std::string mesh_folder = argv[1];
    double dt = std::stod(argv[2]);
    int timesteps = std::stoi(argv[3]);
    int n_splits = num_procs;  // 按MPI进程数分割
    
    if(rank == 0) {
        std::cout << "================================" << std::endl;
        std::cout << "二维欧拉方程求解器 - MPI并行版本" << std::endl;
        std::cout << "================================" << std::endl;
        std::cout << "网格文件夹: " << mesh_folder << std::endl;
        std::cout << "时间步长:   " << dt << std::endl;
        std::cout << "时间步数:   " << timesteps << std::endl;
        std::cout << "MPI进程数:  " << num_procs << std::endl;
        std::cout << "================================" << std::endl;
    }
    
    
    // ========================================
    // 加载完整网格
    // ========================================
    if(rank == 0) {
            std::cout << "正在加载网格..." << std::endl;
    }
    Mesh full_mesh(mesh_folder);
        
    if(rank == 0) {
            std::cout << "网格大小: " << full_mesh.ny << " x " << full_mesh.nx << std::endl;
            std::cout << "网格间距: " << full_mesh.da << std::endl;
            std::cout << "Gamma值:  " << full_mesh.gamma << std::endl;
    }
        
    // ========================================
    // 垂直分割网格给各进程
    // ========================================
    std::vector<Mesh> sub_meshes = splitMeshVertically(full_mesh, n_splits);
    Mesh local_mesh = sub_meshes[rank];
        
    if(rank == 0) {
        std::cout << "网格分割完成: " << n_splits << " 个子域" << std::endl;
    }
    
    // 在时间循环前添加
    double global_umax = computeMaxSpeed(local_mesh, local_mesh.gamma);
    MPI_Allreduce(MPI_IN_PLACE, &global_umax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    double cfl = global_umax * dt / local_mesh.da;
    if(rank == 0) {
    std::cout << "CFL数 = " << cfl << " (应 ≤ 0.5 确保稳定)" << std::endl;
    if(cfl > 0.5) {
        std::cerr << "警告: CFL数过大，可能不稳定！" << std::endl;
    }
    }
        int local_nx = local_mesh.nx;
        int local_ny = local_mesh.ny;
        
        // ========================================
        // 时间循环
        // ========================================
        if(rank == 0) {
            std::cout << "开始时间步进..." << std::endl;
        }
        
        for(int step = 0; step < timesteps; ++step) {
            
            // 更新本地网格
            updateMesh(local_mesh, dt, rank, num_procs);
            
            // 周期性输出
            if(rank == 0) {
                std::cout << "当前时间步: " << step + 1 << " / " << timesteps << std::endl;
            }

        // 修改输出逻辑
        if((step + 1) % 100 == 0) {
        recoverPrimitives(local_mesh);
        std::string timestep_folder = "result/" + std::to_string(step+1);
        if(rank == 0) {
        fs::create_directories(timestep_folder);
         }
         MPI_Barrier(MPI_COMM_WORLD);
         saveMeshData(local_mesh, rank, timestep_folder);
         }
            
        }
        
        
    
    MPI_Finalize();
    return 0;
}