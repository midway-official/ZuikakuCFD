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
#include <iomanip>
#include <numeric>      // std::accumulate
// ============================================================================
// 命名空间 / 类型别名
// ============================================================================
using namespace std;
using namespace Eigen;
namespace fs = std::filesystem;
// ============================================================================
// 网格类
// ============================================================================
class Mesh
{
public:
    // -----------------------------------------------------------------------
    // 网格参数
    // -----------------------------------------------------------------------
    int    nx;       ///< x 方向格子数（含 ghost）
    int    ny;       ///< y 方向格子数
    double da;       ///< 均匀格子尺寸（dx = dy = da）
    double gamma;    ///< 比热比
 
    // -----------------------------------------------------------------------
    // 原始变量
    // -----------------------------------------------------------------------
    MatrixXd rho;    ///< 密度
    MatrixXd u;      ///< x 方向速度
    MatrixXd v;      ///< y 方向速度
    MatrixXd p;      ///< 压力
 
    // -----------------------------------------------------------------------
    // 守恒变量   U = [rho, rho*u, rho*v, E]^T
    // -----------------------------------------------------------------------
    MatrixXd U0;     ///< 密度          rho
    MatrixXd U1;     ///< x 动量        rho*u
    MatrixXd U2;     ///< y 动量        rho*v
    MatrixXd U3;     ///< 总能量        E
 
    // -----------------------------------------------------------------------
    // 边界类型标记矩阵
    //   0  : 内部流体单元（参与计算）
    //  -1  : 固壁 ghost cell（复制法则）
    //  -3  : MPI 进程间 ghost cell（通信填充）
    // -----------------------------------------------------------------------
    MatrixXi bctype;
 
    // -----------------------------------------------------------------------
    // 构造函数
    // -----------------------------------------------------------------------
 
    /**
     * @brief 尺寸初始化：所有物理量置零
     * @param n_y    y 方向格子数
     * @param n_x    x 方向格子数
     * @param da_    均匀格子尺寸
     * @param gamma_ 比热比
     */
    Mesh(int n_y, int n_x, double da_, double gamma_);
 
    /**
     * @brief 从文件夹读取网格（params.txt + 各物理量 .dat 文件）
     * @param folderPath 存放 params.txt / rho.dat / u.dat / v.dat /
     *                   p.dat / bctype.dat 的目录路径
     */
    explicit Mesh(const std::string& folderPath);
};
 
// ============================================================================
// 网格分割（MPI 域分解）
// ============================================================================
 
/**
 * @brief 将全局网格沿 x 方向切分为 n 个子网格（各含 3 层 ghost）
 * @param original 全局网格
 * @param n        进程数（切分份数）
 * @return 子网格向量，长度为 n
 */
vector<Mesh> splitMeshVertically(const Mesh& original, int n);
 
// ============================================================================
// MPI 通信：守恒变量边界列交换
// ============================================================================
 
/**
 * @brief 交换单个矩阵的左右各 3 列 ghost 数据
 * @param matrix     需要交换的矩阵（行主序：行=y，列=x）
 * @param rank       当前进程号
 * @param num_procs  总进程数
 */
void exchangeColumns(MatrixXd& matrix, int rank, int num_procs);
 
/**
 * @brief 对 mesh 的四个守恒变量矩阵分别调用 exchangeColumns
 */
void exchangeConservativeColumns(Mesh& mesh, int rank, int num_procs);
 
// ============================================================================
// 重构格式
// ============================================================================
 
/**
 * @brief MUSCL 二阶重构（配合斜率限制器）
 *
 * 重构 [UL1 | UP] 界面（即 i-1/2 到 i+1/2 之间）处的左右状态。
 * 模板：UL2  UL1  UR1  UR2
 *        i-2  i-1  i+1  i+2
 *
 * @param UL2,UL1,UR1,UR2  四点模板值
 * @param limiter           0=minmod, 1=van Leer, 2=superbee
 * @param UL                输出：界面左状态
 * @param UR                输出：界面右状态
 */
void muscl_reconstruct(
    double UL2, double UL1, double UR1, double UR2,
    int limiter,
    double& UL, double& UR);
 
/**
 * @brief WENO5 五阶重构（Jiang-Shu 光滑指示子）
 *
 * 重构 UP 单元右界面 (i+1/2) 处的左右状态。
 * 模板：UL2   UL1   UP   UR1   UR2   UR3
 *        i-2   i-1   i   i+1   i+2   i+3
 *                         ^
 *                    重构此界面
 *
 * @param UL2~UR3  六点模板值
 * @param UL       输出：界面左状态（左偏重构）
 * @param UR       输出：界面右状态（右偏重构）
 */
void weno5_reconstruct(
    double UL2, double UL1, double UP,
    double UR1, double UR2, double UR3,
    double& UL, double& UR);
 
// ============================================================================
// 数值通量
// ============================================================================
 
/**
 * @brief HLLC 黎曼求解器（x 方向，法向速度为 U1/U0）
 *
 * 输入守恒变量左右状态 [rho, rho*u_n, rho*u_t, E]，
 * 输出 x 方向数值通量 [F0, F1, F2, F3]。
 *
 * @note y 方向通量：调用前需将 U1/U2（法向/切向）互换后传入，
 *       得到的 F1/F2 同样需要互换后写回。
 */
void hllcFlux(
    double UL0, double UL1, double UL2, double UL3,
    double UR0, double UR1, double UR2, double UR3,
    double gamma,
    double& F0, double& F1, double& F2, double& F3);
 
// ============================================================================
// 单元更新（WENO5 + HLLC，有限体积格式）
// ============================================================================
 
/**
 * @brief 对中心单元 Up 做一步 FVM 更新（返回新守恒变量）
 *
 * 各邻格命名约定（以 Up 所在单元为原点）：
 *   Ur1/Ur2/Ur3 : x 方向右侧 +1/+2/+3 格
 *   Ul1/Ul2/Ul3 : x 方向左侧 -1/-2/-3 格
 *   Uu1/Uu2/Uu3 : y 方向上侧 -1/-2/-3 格（行索引减小方向）
 *   Ud1/Ud2/Ud3 : y 方向下侧 +1/+2/+3 格（行索引增大方向）
 *
 * @return 更新后的守恒变量向量 [U0, U1, U2, U3]
 */
vector<double> updateCenterCell(
    const std::vector<double>& Up,
    const std::vector<double>& Ur1,
    const std::vector<double>& Ur2,
    const std::vector<double>& Ur3,
    const std::vector<double>& Ul1,
    const std::vector<double>& Ul2,
    const std::vector<double>& Ul3,
    const std::vector<double>& Uu1,
    const std::vector<double>& Uu2,
    const std::vector<double>& Uu3,
    const std::vector<double>& Ud1,
    const std::vector<double>& Ud2,
    const std::vector<double>& Ud3,
    double gamma,
    double dt,
    double dx,
    double dy);
 
// ============================================================================
// 时间推进（SSP-RK2）
// ============================================================================
 
/**
 * @brief 对网格做一步 SSP-RK2 时间推进（含 MPI ghost 列交换）
 *
 * Stage 1: U*      = U^n + dt·L(U^n)
 * Stage 2: U^{n+1} = 1/2·U^n + 1/2·(U* + dt·L(U*))
 *
 * @param mesh      当前网格（原地修改）
 * @param dt        时间步长
 * @param rank      当前 MPI 进程号
 * @param num_procs 总进程数
 */
void updateMesh(Mesh& mesh, double dt, int rank, int num_procs);
void computeRHS(
    Mesh& mesh,
    double dt,
    MatrixXd& dU0, MatrixXd& dU1, MatrixXd& dU2, MatrixXd& dU3);
// ============================================================================
// 后处理
// ============================================================================
 
/**
 * @brief 从守恒变量 U0~U3 恢复原始变量 rho/u/v/p（全场）
 */
void recoverPrimitives(Mesh& mesh);
 
/**
 * @brief 将网格守恒变量以文本格式保存到指定目录
 *
 * 输出文件：U0_<rank>.dat  U1_<rank>.dat  U2_<rank>.dat  U3_<rank>.dat
 *
 * @param mesh              网格对象
 * @param rank              当前 MPI 进程号（用于文件命名）
 * @param timestep_folder   目标目录路径（空字符串则写当前目录）
 */
void saveMeshData(
    const Mesh& mesh,
    int rank,
    const std::string& timestep_folder);
 
// ============================================================================
// CFL 辅助
// ============================================================================
 
/**
 * @brief 计算网格内最大特征速度（|u|+|v|+a 的最大值）
 * @return 最大波速（用于 CFL 限制时间步）
 */
double computeMaxSpeed(const Mesh& mesh, double gamma);

#endif // FLUID_H