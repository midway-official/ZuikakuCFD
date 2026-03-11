#ifndef FLUID_H
#define FLUID_H

#include <iostream>
#include <iomanip>
#include <eigen3/Eigen/Sparse>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <vector>
#include <filesystem>
#include <mpi.h>
#include <omp.h>
#include <chrono>
using namespace Eigen;
using namespace std;
namespace fs = std::filesystem;


// ============================================================================
// Mesh 类 —— 网格数据容器
// ============================================================================


class Mesh {
public:
    //守恒变量
    MatrixXd U0; // rho     
    MatrixXd U1; // rho*u
    MatrixXd U2; // rho*v
    MatrixXd U3; // E
    //物理量
    MatrixXd rho;      ///< 密度矩阵，ny×nx
    MatrixXd u;        ///< x 方向速度矩阵，ny×nx
    MatrixXd v;        ///< y 方向速度矩阵，ny×nx
    MatrixXd p;        ///< 压力矩阵，ny×nx
    // ── 网格拓扑与边界 ────────────────────────────────────────────────────
    MatrixXi bctype;  ///< 边界类型标记，ny×nx（编码见文件头注释）

    int nx;           ///< x 方向单元数
    int ny;           ///< y 方向单元数
    double da;        //网格边长
    double gamma;        //比热比
    // ── 构造函数 ──────────────────────────────────────────────────────────

    /** @brief 默认构造函数（允许延迟初始化） */
    Mesh() = default;

    /**
     * @brief 按尺寸构造空网格，所有矩阵分配内存但不初始化
     * @param n_y  y 方向单元数
     * @param n_x  x 方向单元数
     * @param da   网格边长
     * @param gamma  比热比
     */
    Mesh(int n_y, int n_x, double da, double gamma);

    Mesh(const std::string& folderPath);

    
};


// ============================================================================
// 并行计算相关函数
// ============================================================================

/**
 * @brief 将完整网格沿 x 方向分割为 n 个子网格（MPI 域分解）
 *
 * @details
 * 分割策略：
 * 1. 尽量均匀分配：widths[k] ≈ nx / n
 * 2. 在接口处各添加 2 列 ghost 层（bctype=-3），用于通信边界插值
 * 3. 子网格完整继承原始网格的 zoneu/zonev 边界速度配置
 * 4. 各子网格独立调用 initGeometry() 和 createInterId()
 *
 * 子网格列数：real_w + left_ghost(0或2) + right_ghost(0或2)
 *
 * @param original  原始完整网格（只读）
 * @param n         分割数（通常等于 MPI 进程数）
 * @return          长度为 n 的子网格向量，sub_meshes[rank] 分配给对应进程
 */
vector<Mesh> splitMeshVertically(const Mesh& original, int n);

// ============================================================================
// 数据通信函数
// ============================================================================

/**
 * @brief 在相邻 MPI 进程间交换 ghost 列数据（列方向域分解专用）
 *
 * @details
 * 通信模式（双向 Sendrecv）：
 * - 向左邻进程发送本进程第 [2,3] 列，接收并写入第 [0,1] 列（左 ghost）
 * - 向右邻进程发送本进程第 [nx-4,nx-3] 列，接收并写入第 [nx-2,nx-1] 列（右 ghost）
 *
 * 最左/最右进程使用 MPI_PROC_NULL 跳过不存在方向的通信。
 *
 * @param matrix     待交换的场变量矩阵（行优先，ny×nx），原地更新 ghost 列
 * @param rank       当前进程编号
 * @param num_procs  总进程数
 *
 * @note 矩阵须为 Eigen ColMajor 格式（默认），列内数据连续，可直接映射为发送缓冲区
 * @note 使用固定 tag（0 和 1），多字段并发通信时需注意 tag 冲突
 */
void exchangeColumns(MatrixXd& matrix, int rank, int num_procs);

void muscl_reconstruct(
    double UL2, double UL1, double UR1, double UR2,
    int limiter,
    double& UL, double& UR);

vector<double> updateCenterCell(
    const std::vector<double>& Up,
    const std::vector<double>& Ur1,
    const std::vector<double>& Ur2,
    const std::vector<double>& Ul1,
    const std::vector<double>& Ul2,
    const std::vector<double>& Uu1,
    const std::vector<double>& Uu2,
    const std::vector<double>& Ud1,
    const std::vector<double>& Ud2,
    double gamma,
    double dt,
    double dx,
    double dy,
    int limiter= 2);
void hllcFlux(
    double UL0, double UL1, double UL2, double UL3,
    double UR0, double UR1, double UR2, double UR3,
    double gamma,
    double& F0, double& F1, double& F2, double& F3);    
void exchangeConservativeColumns(Mesh& mesh, int rank, int num_procs);

void updateMesh(Mesh& mesh, double dt,int rank, int num_procs);      
  
void recoverPrimitives(Mesh& mesh);
void saveMeshData(
    const Mesh& mesh,
    int rank,
    const std::string& timestep_folder);
double computeMaxSpeed(const Mesh& mesh, double gamma);    
#endif // FLUID_H