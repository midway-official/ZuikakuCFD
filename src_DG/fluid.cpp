// =============================================================================
// fluid_dg.cpp  —  间断伽辽金（DG）版本
//
// 核心改动（对比 WENO5 有限体积版本）：
//
//  ┌──────────────────────────────────────────────────────────────────────┐
//  │  WENO5 FV            │  DG（本文件）                                │
//  ├──────────────────────┼───────────────────────────────────────────────┤
//  │  每格存单元平均值    │  每格存 DG_NM 个模态系数（模式 0 = 平均值）  │
//  │  WENO5 重构界面值    │  多项式直接在界面处求值（无需外部模板）       │
//  │  有限体积 RHS        │  弱形式：体积分 + 面积分（HLLC 数值通量）    │
//  │  无限制器（TVD 隐含）│  Cockburn-Shu 矩限制器（消除 Gibbs 振荡）    │
//  │  ghost 深度 3        │  ghost 深度 1 即够（DG 无宽模板依赖）         │
//  └──────────────────────┴───────────────────────────────────────────────┘
//
// DG 方程（对每个单元 K 和测试函数 φ_k）：
//
//   M_k · dû_k/dt = ∫∫_K (F ∂φ_k/∂ξ + G ∂φ_k/∂η) dξ dη    ← 体积分（通量 ∇·分部积分）
//                 - ∮_∂K F*·n φ_k dγ                         ← 面积分（HLLC 数值通量）
//
//   M_k = (h/2)² / ((2px+1)(2py+1))  （Legendre 正交归一质量矩阵，对角）
//
// 参考单元映射：ξ = 2(x - x_c)/h，η = 2(y - y_c)/h，ξ,η ∈ [-1,1]
// =============================================================================

#include "fluid.h"
#include <omp.h>

// ============================================================================
// DG 辅助：Legendre 多项式基函数
// ============================================================================

// L_n(x)，n ≤ DG_P
static inline double legP(int n, double x) noexcept {
    switch (n) {
        case 0: return 1.0;
        case 1: return x;
        case 2: return 0.5 * (3.0*x*x - 1.0);
        case 3: return 0.5 * (5.0*x*x*x - 3.0*x);
        default: return 0.0;
    }
}

// dL_n/dx
static inline double legD(int n, double x) noexcept {
    switch (n) {
        case 0: return 0.0;
        case 1: return 1.0;
        case 2: return 3.0 * x;
        case 3: return 0.5 * (15.0*x*x - 3.0);
        default: return 0.0;
    }
}

// ============================================================================
// DG 辅助：张量积基函数索引与求值
//
// 模式 m = px*(DG_P+1) + py
// φ_m(ξ,η) = L_{px}(ξ) · L_{py}(η)
// ============================================================================
static inline int midx(int px, int py) noexcept { return px*(DG_P+1) + py; }
static inline int mpx (int m)          noexcept { return m / (DG_P+1); }
static inline int mpy (int m)          noexcept { return m % (DG_P+1); }

 double phi2(int m, double xi, double eta) noexcept {
    return legP(mpx(m), xi) * legP(mpy(m), eta);
}

// 参考坐标梯度 ∂φ_m/∂ξ 和 ∂φ_m/∂η
static inline void dphi2(int m, double xi, double eta,
                          double& dxi, double& det) noexcept {
    dxi = legD(mpx(m), xi) * legP(mpy(m), eta);
    det = legP(mpx(m), xi) * legD(mpy(m), eta);
}

// 逆质量系数：(2px+1)*(2py+1)
// 推导：M_mm = (h/2)² * 1/((2px+1)(2py+1))  →  invM_m = (2px+1)(2py+1)
static inline double invMass(int m) noexcept {
    return (double)((2*mpx(m)+1) * (2*mpy(m)+1));
}

// ============================================================================
// DG 辅助：3 点 Gauss-Legendre 求积（精确到 5 次多项式）
//
// 用于体积分（NQ×NQ = 9 个积分点）和面积分（NQ = 3 个积分点）
// P=2 时通量最高次约为 2P=4，3 点 Gauss 满足精度要求
// ============================================================================
static constexpr int NQ = DG_P + 1;   // = 3 for P=2

struct GaussQuadrature {
    // 存储 2阶, 3阶, 4阶 的节点
    static constexpr double POINTS[5][4] = {
        {}, // 0阶 (不使用)
        {}, // 1阶 (不使用)
        {-0.5773502691896257, 0.5773502691896257, 0, 0}, // NQ=2
        {-0.7745966692414834, 0.0, 0.7745966692414834, 0}, // NQ=3
        {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526} // NQ=4
    };

    static constexpr double WEIGHTS[5][4] = {
        {}, 
        {},
        {1.0, 1.0, 0, 0}, // NQ=2
        {0.5555555555555556, 0.8888888888888888, 0.5555555555555556, 0}, // NQ=3
        {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538} // NQ=4
    };
};

// 使用示例：
static const double* GL_PT = GaussQuadrature::POINTS[NQ];
static const double* GL_WT = GaussQuadrature::WEIGHTS[NQ];
// ============================================================================
// DG 辅助：Euler 通量函数
//
// 守恒变量 U = [ρ, ρu, ρv, E]
// x 通量 F = [ρu, ρu²+p, ρuv, (E+p)u]
// y 通量 G = [ρv, ρuv, ρv²+p, (E+p)v]
// ============================================================================
static inline void eulerFluxX(const double U[4], double F[4], double gam) noexcept {
    constexpr double eps = 1e-3;
    double rho = std::max(U[0], eps);
    double u   = U[1] / rho,  v = U[2] / rho;
    double p   = std::max((gam-1.0)*(U[3] - 0.5*rho*(u*u+v*v)), eps);
    F[0] = rho*u;
    F[1] = rho*u*u + p;
    F[2] = rho*u*v;
    F[3] = (U[3]+p)*u;
}

static inline void eulerFluxY(const double U[4], double G[4], double gam) noexcept {
    constexpr double eps = 1e-3;
    double rho = std::max(U[0], eps);
    double u   = U[1] / rho,  v = U[2] / rho;
    double p   = std::max((gam-1.0)*(U[3] - 0.5*rho*(u*u+v*v)), eps);
    G[0] = rho*v;
    G[1] = rho*u*v;
    G[2] = rho*v*v + p;
    G[3] = (U[3]+p)*v;
}

// ============================================================================
// 工具函数（与原版相同）
// ============================================================================
static void readMatrix(const std::string& filename, MatrixXd& mat)
{
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Cannot open file: " + filename);
    for (int i = 0; i < mat.rows(); ++i)
        for (int j = 0; j < mat.cols(); ++j)
            file >> mat(i, j);
}

static void readMatrixInt(const std::string& filename, MatrixXi& mat)
{
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Cannot open file: " + filename);
    for (int i = 0; i < mat.rows(); ++i)
        for (int j = 0; j < mat.cols(); ++j)
            file >> mat(i, j);
}

// ============================================================================
// Mesh::initDGModes
//
// 从当前单元平均值（U0-U3）初始化 DG 自由度：
//   - 模式 0（单元平均）= 现有 U0-U3
//   - 模式 1..DG_NM-1  = 零（分段常数初始条件，高阶模式由时间推进建立）
// ============================================================================
void Mesh::initDGModes()
{
    for (int c = 0; c < 4; ++c) {
        dof[c].assign(DG_NM, MatrixXd::Zero(ny, nx));
    }
    dof[0][0] = U0;
    dof[1][0] = U1;
    dof[2][0] = U2;
    dof[3][0] = U3;
}

// ============================================================================
// Mesh::syncCellAverages
//
// 将 DG 模式 0（单元平均值）同步回 U0-U3
// 必须在 recoverPrimitives、computeMaxSpeed、saveMeshData 之前调用
// ============================================================================
void Mesh::syncCellAverages()
{
    U0 = dof[0][0];
    U1 = dof[1][0];
    U2 = dof[2][0];
    U3 = dof[3][0];
}

// ============================================================================
// 构造函数
// ============================================================================
Mesh::Mesh(int n_y, int n_x, double da_, double gamma_)
    : ny(n_y), nx(n_x), da(da_), gamma(gamma_)
{
    rho    = MatrixXd::Zero(ny, nx);
    u      = MatrixXd::Zero(ny, nx);
    v      = MatrixXd::Zero(ny, nx);
    p      = MatrixXd::Zero(ny, nx);
    U0     = MatrixXd::Zero(ny, nx);
    U1     = MatrixXd::Zero(ny, nx);
    U2     = MatrixXd::Zero(ny, nx);
    U3     = MatrixXd::Zero(ny, nx);
    bctype = MatrixXi::Zero(ny, nx);
    initDGModes();
}

Mesh::Mesh(const std::string& folderPath)
{
    std::ifstream paramFile(folderPath + "/params.txt");
    if (!paramFile) throw std::runtime_error("Cannot open params.txt");
    paramFile >> nx >> ny >> da >> gamma;

    rho    = MatrixXd(ny, nx); u = MatrixXd(ny, nx);
    v      = MatrixXd(ny, nx); p = MatrixXd(ny, nx);
    U0     = MatrixXd(ny, nx); U1 = MatrixXd(ny, nx);
    U2     = MatrixXd(ny, nx); U3 = MatrixXd(ny, nx);
    bctype = MatrixXi(ny, nx);

    readMatrix(folderPath + "/rho.dat", rho);
    readMatrix(folderPath + "/u.dat",   u);
    readMatrix(folderPath + "/v.dat",   v);
    readMatrix(folderPath + "/p.dat",   p);
    readMatrixInt(folderPath + "/bctype.dat", bctype);

    for (int i = 0; i < ny; ++i)
        for (int j = 0; j < nx; ++j) {
            double r=rho(i,j), ux=u(i,j), vy=v(i,j), pr=p(i,j);
            U0(i,j) = r;
            U1(i,j) = r*ux;
            U2(i,j) = r*vy;
            U3(i,j) = pr/(gamma-1.0) + 0.5*r*(ux*ux+vy*vy);
        }

    initDGModes();
}

// ============================================================================
// splitMeshVertically（与原版相同，末尾调用 initDGModes）
// ============================================================================
vector<Mesh> splitMeshVertically(const Mesh& original, int n)
{
    vector<Mesh> sub_meshes;
    const int Nx = original.nx, Ny = original.ny;

    vector<int> widths(n);
    int remain = Nx;
    for (int k = 0; k < n; ++k) { widths[k] = remain/(n-k); remain -= widths[k]; }

    int start = 0;
    for (int k = 0; k < n; ++k) {
        int real_w = widths[k];
        bool has_left  = (k > 0), has_right = (k < n-1);
        int left_ghost  = has_left  ? 3 : 0;
        int right_ghost = has_right ? 3 : 0;
        int sub_nx = real_w + left_ghost + right_ghost;

        Mesh sub(Ny, sub_nx, original.da, original.gamma);
        int real_start = left_ghost;

        for (int i = 0; i < Ny; ++i) {
            for (int j = 0; j < real_w; ++j) {
                int oj = start+j, sj = real_start+j;
                sub.rho(i,sj)=original.rho(i,oj); sub.u(i,sj)=original.u(i,oj);
                sub.v(i,sj)=original.v(i,oj);     sub.p(i,sj)=original.p(i,oj);
                sub.U0(i,sj)=original.U0(i,oj);   sub.U1(i,sj)=original.U1(i,oj);
                sub.U2(i,sj)=original.U2(i,oj);   sub.U3(i,sj)=original.U3(i,oj);
                sub.bctype(i,sj)=original.bctype(i,oj);
                // 复制 DG 模式（若原始网格已有高阶模式）
                for (int c = 0; c < 4; ++c)
                    for (int m = 0; m < DG_NM; ++m)
                        sub.dof[c][m](i,sj) = original.dof[c][m](i,oj);
            }
            if (has_left)
                for (int gj = 0; gj < left_ghost; ++gj) {
                    int oj = std::max(start-left_ghost+gj, 0), sj = gj;
                    sub.rho(i,sj)=original.rho(i,oj); sub.u(i,sj)=original.u(i,oj);
                    sub.v(i,sj)=original.v(i,oj);     sub.p(i,sj)=original.p(i,oj);
                    sub.U0(i,sj)=original.U0(i,oj);   sub.U1(i,sj)=original.U1(i,oj);
                    sub.U2(i,sj)=original.U2(i,oj);   sub.U3(i,sj)=original.U3(i,oj);
                    sub.bctype(i,sj) = -3;
                    for (int c = 0; c < 4; ++c)
                        for (int m = 0; m < DG_NM; ++m)
                            sub.dof[c][m](i,sj) = original.dof[c][m](i,oj);
                }
            if (has_right)
                for (int gj = 0; gj < right_ghost; ++gj) {
                    int oj = std::min(start+real_w+gj, Nx-1), sj = real_start+real_w+gj;
                    sub.rho(i,sj)=original.rho(i,oj); sub.u(i,sj)=original.u(i,oj);
                    sub.v(i,sj)=original.v(i,oj);     sub.p(i,sj)=original.p(i,oj);
                    sub.U0(i,sj)=original.U0(i,oj);   sub.U1(i,sj)=original.U1(i,oj);
                    sub.U2(i,sj)=original.U2(i,oj);   sub.U3(i,sj)=original.U3(i,oj);
                    sub.bctype(i,sj) = -3;
                    for (int c = 0; c < 4; ++c)
                        for (int m = 0; m < DG_NM; ++m)
                            sub.dof[c][m](i,sj) = original.dof[c][m](i,oj);
                }
        }
        sub_meshes.push_back(sub);
        start += real_w;
    }
    return sub_meshes;
}

// ============================================================================
// exchangeColumns（单矩阵版本，与原版相同）
// ============================================================================
void exchangeColumns(MatrixXd& matrix, int rank, int num_procs)
{
    const int rows  = matrix.rows();
    const int cols  = matrix.cols();
    const int count = rows * 3;

    int left_rank  = (rank == 0)             ? MPI_PROC_NULL : rank - 1;
    int right_rank = (rank == num_procs - 1) ? MPI_PROC_NULL : rank + 1;

    const double* send_left_ptr  = matrix.data() + 3 * rows;
    const double* send_right_ptr = matrix.data() + (cols-6) * rows;

    vector<double> recv_left(count), recv_right(count);

    MPI_Sendrecv(send_left_ptr,  count, MPI_DOUBLE, left_rank,  0,
                 recv_left.data(),  count, MPI_DOUBLE, left_rank,  1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(send_right_ptr, count, MPI_DOUBLE, right_rank, 1,
                 recv_right.data(), count, MPI_DOUBLE, right_rank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (left_rank  != MPI_PROC_NULL)
        std::memcpy(matrix.data(), recv_left.data(), count*sizeof(double));
    if (right_rank != MPI_PROC_NULL)
        std::memcpy(matrix.data()+(cols-3)*rows, recv_right.data(), count*sizeof(double));
}
void exchangeConservativeColumns(Mesh& mesh, int rank, int num_procs)
{
    const int rows = mesh.ny;
    const int cols = mesh.nx;
    
    // 基础单元：每个 (variable, mode) 需要交换的列数（通常是 3 列 ghost cells）
    const int ghost_width = 3;
    const size_t count_per_field = static_cast<size_t>(rows) * ghost_width;
    const size_t total_elements = 4 * DG_NM * count_per_field;

    int left_rank  = (rank == 0)             ? MPI_PROC_NULL : rank - 1;
    int right_rank = (rank == num_procs - 1) ? MPI_PROC_NULL : rank + 1;

    // 预留足够的空间
    std::vector<double> send_to_left(total_elements, 0.0);
    std::vector<double> send_to_right(total_elements, 0.0);
    std::vector<double> recv_from_left(total_elements, 0.0);
    std::vector<double> recv_from_right(total_elements, 0.0);

    // --- 1. 谨慎打包 (Packing) ---
    for (int c = 0; c < 4; ++c) {
        for (int m = 0; m < DG_NM; ++m) {
            size_t offset = (static_cast<size_t>(c) * DG_NM + m) * count_per_field;
            const double* field_data = mesh.dof[c][m].data();

            // 发送给左邻居的数据：取当前块的“内部左边界”（跳过 ghost 列）
            // 假设 matrix.data() 布局为 [col0][col1]... 
            // 如果 exchangeColumns 取的是 matrix.data() + 3*rows，这里也应保持一致
            std::memcpy(&send_to_left[offset], 
                        field_data + ghost_width * rows, 
                        count_per_field * sizeof(double));

            // 发送给右邻居的数据：取当前块的“内部右边界”
            std::memcpy(&send_to_right[offset], 
                        field_data + (cols - 2 * ghost_width) * rows, 
                        count_per_field * sizeof(double));
        }
    }

    // --- 2. 稳健通信 (Communication) ---
    // 使用明确的 Tag 区分左右方向，防止在某些 MPI 实现下发生混淆
    MPI_Request requests[4];
    
    // 向左交换
    MPI_Sendrecv(send_to_left.data(),  (int)total_elements, MPI_DOUBLE, left_rank,  10,
                 recv_from_left.data(), (int)total_elements, MPI_DOUBLE, left_rank,  11,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // 向右交换
    MPI_Sendrecv(send_to_right.data(), (int)total_elements, MPI_DOUBLE, right_rank, 11,
                 recv_from_right.data(), (int)total_elements, MPI_DOUBLE, right_rank, 10,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // --- 3. 谨慎解包 (Unpacking) ---
    for (int c = 0; c < 4; ++c) {
        for (int m = 0; m < DG_NM; ++m) {
            size_t offset = (static_cast<size_t>(c) * DG_NM + m) * count_per_field;
            double* field_data = mesh.dof[c][m].data();

            // 如果左边有邻居，将收到的数据放入最左侧 ghost 区域 (col 0~2)
            if (left_rank != MPI_PROC_NULL) {
                std::memcpy(field_data, 
                            &recv_from_left[offset], 
                            count_per_field * sizeof(double));
            }

            // 如果右边有邻居，将收到的数据放入最右侧 ghost 区域 (col cols-3~cols-1)
            if (right_rank != MPI_PROC_NULL) {
                std::memcpy(field_data + (cols - ghost_width) * rows, 
                            &recv_from_right[offset], 
                            count_per_field * sizeof(double));
            }
        }
    }

    // --- 4. 基础变量同步 --

    mesh.syncCellAverages();
}

// ============================================================================
// MUSCL 重构（保留，DG 主循环不使用）
// ============================================================================
void muscl_reconstruct(double UL2, double UL1, double UR1, double UR2,
                       int limiter, double& UL, double& UR)
{
    double dL=UL1-UL2, dC=UR1-UL1, dR=UR2-UR1;
    double sigma_L=0.0, sigma_R=0.0;
    auto minmod2=[](double a,double b)->double{
        return a*b<=0.0?0.0:(std::abs(a)<std::abs(b)?a:b);};
    auto van_leer=[](double a,double b)->double{
        double ab=a*b; return ab<=0.0?0.0:2.0*ab/(a+b);};
    auto superbee=[&](double a,double b)->double{
        double s1=minmod2(a,2.0*b),s2=minmod2(2.0*a,b);
        return s1*s2<=0.0?0.0:(std::abs(s1)>std::abs(s2)?s1:s2);};
    switch(limiter){
        case 0: sigma_L=minmod2(dL,dC); sigma_R=minmod2(dC,dR); break;
        case 1: sigma_L=van_leer(dL,dC); sigma_R=van_leer(dC,dR); break;
        case 2: sigma_L=superbee(dL,dC); sigma_R=superbee(dC,dR); break;
        default: break;
    }
    UL = UL1 + 0.5*sigma_L;
    UR = UR1 - 0.5*sigma_R;
}

// ============================================================================ 
// HLLC 数值通量（DG 面积分调用）
// ============================================================================
void hllcFlux(
    double UL0,double UL1,double UL2,double UL3,
    double UR0,double UR1,double UR2,double UR3,
    double gamma,
    double& F0,double& F1,double& F2,double& F3)
{
    const double eps = 1e-6;

    // --- 左状态 ---
    double rhoL = std::max(UL0, eps);
    double uL   = UL1 / rhoL;
    double vL   = UL2 / rhoL;
    double EL   = UL3 / rhoL; // 单位质量能量
    double pL   = std::max((gamma-1.0)*(UL3 - 0.5*rhoL*(uL*uL + vL*vL)), eps);
    double aL   = sqrt(std::max(gamma*pL/rhoL, eps));

    // --- 右状态 ---
    double rhoR = std::max(UR0, eps);
    double uR   = UR1 / rhoR;
    double vR   = UR2 / rhoR;
    double ER   = UR3 / rhoR;
    double pR   = std::max((gamma-1.0)*(UR3 - 0.5*rhoR*(uR*uR + vR*vR)), eps);
    double aR   = sqrt(std::max(gamma*pR/rhoR, eps));

    // --- 物理通量 ---
    double FL0 = rhoL*uL, FL1 = rhoL*uL*uL + pL, FL2 = rhoL*uL*vL, FL3 = (UL3+pL)*uL;
    double FR0 = rhoR*uR, FR1 = rhoR*uR*uR + pR, FR2 = rhoR*uR*vR, FR3 = (UR3+pR)*uR;

    // --- Roe 平均 ---
    double sqL = sqrt(rhoL), sqR = sqrt(rhoR);
    double den = sqL + sqR;
    if (den < eps) den = eps; // 防止除零
    double u_roe = (sqL*uL + sqR*uR)/den;
    double h_roe = (sqL*(EL+pL/rhoL) + sqR*(ER+pR/rhoR))/den;
    double c_roe = sqrt(std::max((gamma-1.0)*(h_roe - 0.5*u_roe*u_roe), eps));

    // --- 波速 ---
    double SL = std::min(uL-aL, u_roe-c_roe);
    double SR = std::max(uR+aR, u_roe+c_roe);

    if (SL >= 0.0) { F0=FL0; F1=FL1; F2=FL2; F3=FL3; return; }
    if (SR <= 0.0) { F0=FR0; F1=FR1; F2=FR2; F3=FR3; return; }

    double num = pR - pL + rhoL*uL*(SL-uL) - rhoR*uR*(SR-uR);
    double den2 = rhoL*(SL-uL) - rhoR*(SR-uR);
    if (fabs(den2) < eps) den2 = (den2>=0? eps : -eps); // 防止除零
    double Sstar = num / den2;

    auto star_state = [&](double rho,double u,double v,double E,double p,double S) -> std::array<double,4> {
        double denom = S - Sstar;
        if (fabs(denom) < eps) denom = (denom >=0 ? eps : -eps);
        double c = rho*(S-u)/denom;
        double term = (Sstar-u)*(Sstar + p/(rho*(S-u)));
        if (std::isnan(term) || std::isinf(term)) term = 0.0; // 防止 NaN
        return {c, c*Sstar, c*v, c*(E + term)};
    };

    if (Sstar >= 0.0) {
        auto UL_s = star_state(rhoL,uL,vL,EL,pL,SL);
        F0 = FL0 + SL*(UL_s[0]-UL0);
        F1 = FL1 + SL*(UL_s[1]-UL1);
        F2 = FL2 + SL*(UL_s[2]-UL2);
        F3 = FL3 + SL*(UL_s[3]-UL3);
    } else {
        auto UR_s = star_state(rhoR,uR,vR,ER,pR,SR);
        F0 = FR0 + SR*(UR_s[0]-UR0);
        F1 = FR1 + SR*(UR_s[1]-UR1);
        F2 = FR2 + SR*(UR_s[2]-UR2);
        F3 = FR3 + SR*(UR_s[3]-UR3);
    }
}

// ============================================================================
// HLLE 数值通量（高马赫数鲁棒版本）
// ============================================================================
void hlleFlux(
    double UL0, double UL1, double UL2, double UL3,
    double UR0, double UR1, double UR2, double UR3,
    double gamma,
    double& F0, double& F1, double& F2, double& F3)
{
    const double eps = 1e-6;

    // --- 左状态 ---
    double rhoL = std::max(UL0, eps);
    double uL   = UL1 / rhoL;
    double vL   = UL2 / rhoL;
    double EL   = UL3;
    double pL   = std::max((gamma-1.0)*(EL - 0.5*rhoL*(uL*uL + vL*vL)), eps);
    double aL   = sqrt(std::max(gamma*pL/rhoL, eps));
    double HL   = (EL+pL)/rhoL;

    // --- 右状态 ---
    double rhoR = std::max(UR0, eps);
    double uR   = UR1 / rhoR;
    double vR   = UR2 / rhoR;
    double ER   = UR3;
    double pR   = std::max((gamma-1.0)*(ER - 0.5*rhoR*(uR*uR + vR*vR)), eps);
    double aR   = sqrt(std::max(gamma*pR/rhoR, eps));
    double HR   = (ER+pR)/rhoR;

    // --- Roe 平均 ---
    double sqL = sqrt(rhoL), sqR = sqrt(rhoR);
    double invDen = 1.0 / std::max(sqL + sqR, eps);

    double u_roe = (sqL*uL + sqR*uR)*invDen;
    double v_roe = (sqL*vL + sqR*vR)*invDen;
    double h_roe = (sqL*HL + sqR*HR)*invDen;
    double c_roe = sqrt(std::max((gamma-1.0)*(h_roe - 0.5*(u_roe*u_roe + v_roe*v_roe)), eps));

    // --- 波速估计 ---
    double SL = std::min(uL - aL, u_roe - c_roe);
    double SR = std::max(uR + aR, u_roe + c_roe);

    // --- 物理通量 ---
    double FL0 = rhoL*uL, FL1 = rhoL*uL*uL + pL, FL2 = rhoL*uL*vL, FL3 = (EL+pL)*uL;
    double FR0 = rhoR*uR, FR1 = rhoR*uR*uR + pR, FR2 = rhoR*uR*vR, FR3 = (ER+pR)*uR;

    // --- HLLE 数值通量 ---
    if (SL >= 0.0) {
        F0=FL0; F1=FL1; F2=FL2; F3=FL3;
    } else if (SR <= 0.0) {
        F0=FR0; F1=FR1; F2=FR2; F3=FR3;
    } else {
        double invDiff = 1.0 / std::max(SR - SL, eps);
        F0 = (SR*FL0 - SL*FR0 + SL*SR*(UR0-UL0)) * invDiff;
        F1 = (SR*FL1 - SL*FR1 + SL*SR*(UR1-UL1)) * invDiff;
        F2 = (SR*FL2 - SL*FR2 + SL*SR*(UR2-UL2)) * invDiff;
        F3 = (SR*FL3 - SL*FR3 + SL*SR*(UR3-UL3)) * invDiff;
    }
}
static void fillDGGhostCells(Mesh& mesh)
{
    const int ny=mesh.ny, nx=mesh.nx;
    constexpr int OFF[6] = {-3,-2,-1,1,2,3};

    // --- Y 方向 Ghost 填充 (处理顶/底边界) ---
    for (int i = 3; i < ny-3; ++i)
        for (int j = 3; j < nx-3; ++j) {
            if (mesh.bctype(i,j) != 0) continue;
            for (int dk : OFF) {
                int ni = i+dk;
                if ((unsigned)ni >= (unsigned)ny) continue;
                
                int bt = mesh.bctype(ni,j);
                if (bt == -1) { // 零梯度 (Zero-Gradient)
                    for (int c = 0; c < 4; ++c)
                        for (int m = 0; m < DG_NM; ++m)
                            mesh.dof[c][m](ni,j) = mesh.dof[c][m](i,j);
                } 
                else if (bt == 1) { // 滑移边界 (Slip: y是法向)
                    for (int m = 0; m < DG_NM; ++m) {
                        mesh.dof[0][m](ni,j) =  mesh.dof[0][m](i,j); // rho
                        mesh.dof[1][m](ni,j) =  mesh.dof[1][m](i,j); // rho*u (切向)
                        mesh.dof[2][m](ni,j) = -mesh.dof[2][m](i,j); // rho*v (法向取反)
                        mesh.dof[3][m](ni,j) =  mesh.dof[3][m](i,j); // E
                    }
                }
            }
        }

    // --- X 方向 Ghost 填充 (处理左/右边界) ---
    for (int i = 3; i < ny-3; ++i)
        for (int j = 3; j < nx-3; ++j) {
            if (mesh.bctype(i,j) != 0) continue;
            for (int dk : OFF) {
                int nj = j+dk; // 修正：这里应使用 nj
                if ((unsigned)nj >= (unsigned)nx) continue;

                int bt = mesh.bctype(i,nj);
                if (bt == -1) { // 零梯度
                    for (int c = 0; c < 4; ++c)
                        for (int m = 0; m < DG_NM; ++m)
                            mesh.dof[c][m](i,nj) = mesh.dof[c][m](i,j);
                }
                else if (bt == 1) { // 滑移边界 (Slip: x是法向)
                    for (int m = 0; m < DG_NM; ++m) {
                        mesh.dof[0][m](i,nj) =  mesh.dof[0][m](i,j); // rho
                        mesh.dof[1][m](i,nj) = -mesh.dof[1][m](i,j); // rho*u (法向取反)
                        mesh.dof[2][m](i,nj) =  mesh.dof[2][m](i,j); // rho*v (切向)
                        mesh.dof[3][m](i,nj) =  mesh.dof[3][m](i,j); // E
                    }
                }
            }
        }
}

// ============================================================================
// computeRHS — DG 空间离散算子
//
// 对每个内部单元 (i,j)，计算 dû_m^c/dt（各模式系数的时间导数）
//
// ── 体积分（通量 ∇ 分部积分后）────────────────────────────────────────────
//
//   V_m^c = ∑_{a,b} w_a·w_b · [F_c(U(ξ_a,η_b)) ∂φ_m/∂ξ + G_c(U(ξ_a,η_b)) ∂φ_m/∂η]
//
//   其中 U(ξ,η) = ∑_k û_k φ_k(ξ,η)（由 DG_NM 个模式系数展开）
//
// ── 面积分（4 个界面，各用 NQ 个 Gauss 点）───────────────────────────────
//
//   右界面 (ξ=+1)：S_m += ∑_q w_q · F*(U_L(+1,η_q), U_R_right(-1,η_q)) · φ_m(+1,η_q)
//   左界面 (ξ=-1)：S_m -= ∑_q w_q · F*(U_L_left(+1,η_q), U_R(-1,η_q)) · φ_m(-1,η_q)
//   下界面 (η=+1)：S_m += ∑_q w_q · G*(U_L(ξ_q,+1), U_R_down(ξ_q,-1)) · φ_m(ξ_q,+1)
//   上界面 (η=-1)：S_m -= ∑_q w_q · G*(U_L_up(ξ_q,+1), U_R(ξ_q,-1)) · φ_m(ξ_q,-1)
//
// ============================================================================
void computeRHS(Mesh& mesh, double dt,
                vector<MatrixXd> dU[4])
{
    const int ny=mesh.ny, nx=mesh.nx;
    const double h=mesh.da, gam=mesh.gamma;

    // 初始化输出
    for (int c = 0; c < 4; ++c)
        for (int m = 0; m < DG_NM; ++m)
            dU[c][m] = MatrixXd::Zero(ny, nx);

    // 填充固壁 ghost（MPI ghost 已由 exchangeConservativeColumns 填充）
    fillDGGhostCells(mesh);

    // ── 在单元 (i,j) 的参考坐标 (xi,eta) 处重建 U ─────────────────────────
    // U(xi,eta) = ∑_m û_m φ_m(xi,eta)
    auto evalU = [&](int i, int j, double xi, double eta, double U[4]) __attribute__((always_inline)) {
        U[0]=U[1]=U[2]=U[3]=0.0;
        for (int m = 0; m < DG_NM; ++m) {
            double b = phi2(m, xi, eta);
            U[0] += mesh.dof[0][m](i,j)*b;
            U[1] += mesh.dof[1][m](i,j)*b;
            U[2] += mesh.dof[2][m](i,j)*b;
            U[3] += mesh.dof[3][m](i,j)*b;
        }
    };

    const double inv2h = 0.5 / h;   // invMass 归一化因子的一部分

    for (int i = 3; i < ny-3; ++i)
    {
        for (int j = 3; j < nx-3; ++j)
        {
            if (mesh.bctype(i,j) != 0) continue;

            // rhs[c][m] 为"参考空间"积分值（待乘以 invMass*2/h）
            double rhs[4][DG_NM] = {};

            // ──────────────────────────────────────────────────────────────
            // STEP 1：体积分
            // V_m^c = ∑_{a,b} w_a w_b [F_c ∂φ_m/∂ξ + G_c ∂φ_m/∂η]
            // ──────────────────────────────────────────────────────────────
            for (int qa = 0; qa < NQ; ++qa) {
                double xi=GL_PT[qa], wa=GL_WT[qa];
                for (int qb = 0; qb < NQ; ++qb) {
                    double eta=GL_PT[qb], wb=GL_WT[qb], ww=wa*wb;

                    double U[4];
                    evalU(i, j, xi, eta, U);

                    double F[4], G[4];
                    eulerFluxX(U, F, gam);
                    eulerFluxY(U, G, gam);

                    for (int m = 0; m < DG_NM; ++m) {
                        double dxi, det;
                        dphi2(m, xi, eta, dxi, det);
                        double coeff = ww;
                        for (int c = 0; c < 4; ++c)
                            rhs[c][m] += coeff * (F[c]*dxi + G[c]*det);
                    }
                }
            }

            // ──────────────────────────────────────────────────────────────
            // STEP 2：面积分
            //
            // 符号约定（外法向 n）：
            //   右界面 (ξ=+1, n=+x)：rhs -= F* φ_m(+1,η)
            //   左界面 (ξ=-1, n=-x)：rhs += F* φ_m(-1,η)   [n=-x → ∮ 号翻转]
            //   下界面 (η=+1, n=+y)：rhs -= G* φ_m(ξ,+1)
            //   上界面 (η=-1, n=-y)：rhs += G* φ_m(ξ,-1)   [n=-y → ∮ 号翻转]
            // ──────────────────────────────────────────────────────────────

            // ── 右界面 (ξ=+1)，与单元 (i, j+1) 共享 ───────────────────────
            for (int q = 0; q < NQ; ++q) {
                double eta=GL_PT[q], w=GL_WT[q];
                double UL[4], UR[4];
                evalU(i, j,   +1.0, eta, UL);   // 当前格右边界
                evalU(i, j+1, -1.0, eta, UR);   // 右邻格左边界
                double F0,F1,F2,F3;
                hllcFlux(UL[0],UL[1],UL[2],UL[3], UR[0],UR[1],UR[2],UR[3], gam, F0,F1,F2,F3);
                double Fstar[4] = {F0,F1,F2,F3};
                for (int m = 0; m < DG_NM; ++m) {
                    double bv = phi2(m, +1.0, eta) * w;
                    for (int c = 0; c < 4; ++c)
                        rhs[c][m] -= bv * Fstar[c];   // 外法向 +x：减号
                }
            }

            // ── 左界面 (ξ=-1)，与单元 (i, j-1) 共享 ───────────────────────
            for (int q = 0; q < NQ; ++q) {
                double eta=GL_PT[q], w=GL_WT[q];
                double UL[4], UR[4];
                evalU(i, j-1, +1.0, eta, UL);   // 左邻格右边界
                evalU(i, j,   -1.0, eta, UR);   // 当前格左边界
                double F0,F1,F2,F3;
                hllcFlux(UL[0],UL[1],UL[2],UL[3], UR[0],UR[1],UR[2],UR[3], gam, F0,F1,F2,F3);
                double Fstar[4] = {F0,F1,F2,F3};
                for (int m = 0; m < DG_NM; ++m) {
                    double bv = phi2(m, -1.0, eta) * w;
                    for (int c = 0; c < 4; ++c)
                        rhs[c][m] += bv * Fstar[c];   // 外法向 -x：加号
                }
            }

            // ── 下界面 (η=+1)，与单元 (i+1, j) 共享（y 方向）─────────────
            // y 方向通量 G：调用 hllcFlux 时交换 U1↔U2（使 v 成为法向速度）
            // hllcFlux 返回 (G0, G_n, G_t, G3)，其中 G_n=y动量通量，G_t=x动量通量
            for (int q = 0; q < NQ; ++q) {
                double xi=GL_PT[q], w=GL_WT[q];
                double UL[4], UR[4];
                evalU(i,   j, xi, +1.0, UL);   // 当前格下边界
                evalU(i+1, j, xi, -1.0, UR);   // 下邻格上边界
                double G0,Gn,Gt,G3;
                // 交换 U[1]↔U[2]：令 v 为法向速度（HLLC 要求法向速度在 [1] 位置）
                hllcFlux(UL[0],UL[2],UL[1],UL[3], UR[0],UR[2],UR[1],UR[3], gam, G0,Gn,Gt,G3);
                // 解映射：Gn → rho*v 通量（[2] 分量），Gt → rho*u 通量（[1] 分量）
                double Gstar[4] = {G0, Gt, Gn, G3};
                for (int m = 0; m < DG_NM; ++m) {
                    double bv = phi2(m, xi, +1.0) * w;
                    for (int c = 0; c < 4; ++c)
                        rhs[c][m] -= bv * Gstar[c];   // 外法向 +y：减号
                }
            }

            // ── 上界面 (η=-1)，与单元 (i-1, j) 共享（y 方向）─────────────
            for (int q = 0; q < NQ; ++q) {
                double xi=GL_PT[q], w=GL_WT[q];
                double UL[4], UR[4];
                evalU(i-1, j, xi, +1.0, UL);   // 上邻格下边界
                evalU(i,   j, xi, -1.0, UR);   // 当前格上边界
                double G0,Gn,Gt,G3;
                hllcFlux(UL[0],UL[2],UL[1],UL[3], UR[0],UR[2],UR[1],UR[3], gam, G0,Gn,Gt,G3);
                double Gstar[4] = {G0, Gt, Gn, G3};
                for (int m = 0; m < DG_NM; ++m) {
                    double bv = phi2(m, xi, -1.0) * w;
                    for (int c = 0; c < 4; ++c)
                        rhs[c][m] += bv * Gstar[c];   // 外法向 -y：加号
                }
            }

            // ── STEP 3：质量矩阵归一化 → 写入输出 ─────────────────────────
            // dû_m/dt = rhs_m · invMass(m) · 2/h
            for (int m = 0; m < DG_NM; ++m) {
                double scale = invMass(m) * inv2h;
                for (int c = 0; c < 4; ++c)
                    dU[c][m](i,j) = rhs[c][m] * scale;
            }
        }
    }
}
inline double basis(int m, double xi, double eta)
{
    int px = mpx(m);
    int py = mpy(m);

    // Legendre P0, P1, P2
    auto P = [](int p, double x)
    {
        if (p == 0) return 1.0;
        if (p == 1) return x;
        if (p == 2) return 0.5*(3.0*x*x-1.0);
        return 0.0;
    };

    return P(px,xi) * P(py,eta);
}
void applyLimiter(Mesh& mesh) 
{
    const int ny = mesh.ny;
    const int nx = mesh.nx;

    const double kappa = 1.0;  // smooth width
    const double s0 = -3.0*log10((double)DG_P);

    const double RHO_MIN = 1e-8;
    const double P_MIN   = 1e-8;
    const double EPS     = 1e-8;

    // ============================
    // Persson shock detector
    // ============================
    for(int i=3;i<ny-3;i++)
    for(int j=3;j<nx-3;j++)
    {
        if(mesh.bctype(i,j)!=0) continue;

        for(int c=0;c<4;c++)
        {
            double Etot = 0.0;
            double Ep   = 0.0;

            for(int px=0;px<=DG_P;px++)
            for(int py=0;py<=DG_P;py++)
            {
                int m = px*(DG_P+1)+py;
                double a = mesh.dof[c][m](i,j);

                Etot += a*a;

                if(px==DG_P || py==DG_P)
                    Ep += a*a;
            }

            double s = log10(Ep/(Etot+EPS));

            if(s > s0)
            {
                double eps;
                if(s > s0 + kappa)
                    eps = 1.0;
                else
                {
                    double r = (s - s0)/kappa;
                    eps = 0.5*(1.0 + sin(M_PI*(r-0.5)));
                }

                for(int px=0;px<=DG_P;px++)
                for(int py=0;py<=DG_P;py++)
                {
                    if(px==DG_P || py==DG_P)
                    {
                        int m = px*(DG_P+1)+py;
                        mesh.dof[c][m](i,j) *= (1.0-eps);
                    }
                }
            }
        }
    }

    // =====================================
    // Zhang-Shu positivity limiter
    // =====================================
    const double xi[4]  = {-0.8611363115940526, -0.3399810435848563,
                             0.3399810435848563, 0.8611363115940526};
    const double eta[4] = {-0.8611363115940526, -0.3399810435848563,
                             0.3399810435848563, 0.8611363115940526};

    for(int i=3;i<ny-3;i++)
    for(int j=3;j<nx-3;j++)
    {
        // ------------------------------
        // Step 0: 强制零阶模态非负
        // ------------------------------
        mesh.dof[0][0](i,j) = std::max(mesh.dof[0][0](i,j), RHO_MIN);

        double rho0 = mesh.dof[0][0](i,j);
        double mx0  = mesh.dof[1][0](i,j);
        double my0  = mesh.dof[2][0](i,j);
        double E0   = mesh.dof[3][0](i,j);

        double p0 = (mesh.gamma-1.0)*(E0 - 0.5*(mx0*mx0 + my0*my0)/rho0);
        if(p0 < P_MIN)
        {
            E0 = P_MIN/(mesh.gamma-1.0) + 0.5*(mx0*mx0 + my0*my0)/rho0;
            mesh.dof[3][0](i,j) = E0;
        }

        // ------------------------------
        // Step 1: 检查 rho
        // ------------------------------
        double theta_rho = 1.0;
        for(int a=0;a<4;a++)
        for(int b=0;b<4;b++)
        {
            double rho = 0.0;
            for(int m=0;m<DG_NM;m++)
            {
                rho += mesh.dof[0][m](i,j) * basis(m, xi[a], eta[b]);
            }
            if(rho < RHO_MIN)
            {
                double t = (rho0 - RHO_MIN)/(rho0 - rho + EPS);
                theta_rho = std::min(theta_rho, t);
            }
        }

        if(theta_rho < 1.0)
        {
            for(int m=1;m<DG_NM;m++)
            {
                mesh.dof[0][m](i,j) *= theta_rho;
                mesh.dof[1][m](i,j) *= theta_rho;
                mesh.dof[2][m](i,j) *= theta_rho;
                mesh.dof[3][m](i,j) *= theta_rho;
            }
        }

        // ------------------------------
        // Step 2: 检查 p
        // ------------------------------
        double theta_p = 1.0;
        for(int a=0;a<4;a++)
        for(int b=0;b<4;b++)
        {
            double rho=0.0, mx=0.0, my=0.0, E=0.0;
            for(int m=0;m<DG_NM;m++)
            {
                double phi = basis(m, xi[a], eta[b]);
                rho += mesh.dof[0][m](i,j)*phi;
                mx  += mesh.dof[1][m](i,j)*phi;
                my  += mesh.dof[2][m](i,j)*phi;
                E   += mesh.dof[3][m](i,j)*phi;
            }

            if(rho <= RHO_MIN) continue;

            double u = mx/rho;
            double v = my/rho;
            double p = (mesh.gamma-1.0)*(E - 0.5*rho*(u*u + v*v));

            if(p < P_MIN)
            {
                double p_avg = (mesh.gamma-1.0)*(mesh.dof[3][0](i,j)
                                -0.5*(mesh.dof[1][0](i,j)*mesh.dof[1][0](i,j)
                                      +mesh.dof[2][0](i,j)*mesh.dof[2][0](i,j))/mesh.dof[0][0](i,j));

                double t = (p_avg - P_MIN)/(p_avg - p + EPS);
                theta_p = std::min(theta_p, t);
            }
        }

        if(theta_p < 1.0)
        {
            for(int m=1;m<DG_NM;m++)
            {
                mesh.dof[0][m](i,j) *= theta_p;
                mesh.dof[1][m](i,j) *= theta_p;
                mesh.dof[2][m](i,j) *= theta_p;
                mesh.dof[3][m](i,j) *= theta_p;
            }
        }
    }
}

// ============================================================================
// updateMesh — SSP-RK3 时间推进
//
// 使用 Shu-Osher 三阶强稳定格式：
//   Stage 1: U*      = U^n + dt·L(U^n)
//   Stage 2: U**     = 3/4·U^n + 1/4·(U* + dt·L(U*))
//   Stage 3: U^{n+1} = 1/3·U^n + 2/3·(U** + dt·L(U**))
//
// 每个阶段后的操作：
//   applyLimiter → syncCellAverages → exchangeConservativeColumns（MPI ghost）
// ============================================================================
void updateMesh(Mesh& mesh, double dt, int rank, int num_procs)
{
    const int ny=mesh.ny, nx=mesh.nx;

    // 保存 U^n（所有 DG 模式）
    vector<MatrixXd> dof_n[4];
    for (int c = 0; c < 4; ++c) {
        dof_n[c].resize(DG_NM);
        for (int m = 0; m < DG_NM; ++m)
            dof_n[c][m] = mesh.dof[c][m];
    }

    // RHS 缓冲区
    vector<MatrixXd> dU[4];
    for (int c = 0; c < 4; ++c)
        dU[c].resize(DG_NM, MatrixXd::Zero(ny, nx));
    // ── Stage 1：U* = U^n + dt·L(U^n) ──────────────────────────────────
    computeRHS(mesh, dt, dU);
    exchangeConservativeColumns(mesh, rank, num_procs);
    for (int c = 0; c < 4; ++c)
        for (int m = 0; m < DG_NM; ++m)
            mesh.dof[c][m] = dof_n[c][m] + dt * dU[c][m];
      
    applyLimiter(mesh);
    exchangeConservativeColumns(mesh, rank, num_procs);
    mesh.syncCellAverages();
    exchangeConservativeColumns(mesh, rank, num_procs);

    // ── Stage 2：U** = 3/4·U^n + 1/4·(U* + dt·L(U*)) ───────────────────
    computeRHS(mesh, dt, dU);
    exchangeConservativeColumns(mesh, rank, num_procs);
    for (int c = 0; c < 4; ++c)
        for (int m = 0; m < DG_NM; ++m)
            mesh.dof[c][m] = 0.75*dof_n[c][m] + 0.25*(mesh.dof[c][m] + dt*dU[c][m]);    
    applyLimiter(mesh);
    exchangeConservativeColumns(mesh, rank, num_procs);
    mesh.syncCellAverages();
    exchangeConservativeColumns(mesh, rank, num_procs);

    // ── Stage 3：U^{n+1} = 1/3·U^n + 2/3·(U** + dt·L(U**)) ─────────────
    computeRHS(mesh, dt, dU);
    for (int c = 0; c < 4; ++c)
        for (int m = 0; m < DG_NM; ++m)
            mesh.dof[c][m] = (1.0/3.0)*dof_n[c][m] + (2.0/3.0)*(mesh.dof[c][m] + dt*dU[c][m]);
    exchangeConservativeColumns(mesh, rank, num_procs);        
    applyLimiter(mesh);
    exchangeConservativeColumns(mesh, rank, num_procs);
    mesh.syncCellAverages();
    exchangeConservativeColumns(mesh, rank, num_procs);
}

// ============================================================================
// recoverPrimitives（从单元平均守恒量恢复原始变量，与原版相同）
// ============================================================================
void recoverPrimitives(Mesh& mesh)
{
    const int ny=mesh.ny, nx=mesh.nx;
    const double g=mesh.gamma;
    for (int j = 0; j < nx; ++j)
        for (int i = 0; i < ny; ++i) {
            double rho=mesh.U0(i,j), ux=mesh.U1(i,j)/rho, vy=mesh.U2(i,j)/rho;
            mesh.rho(i,j)=rho; mesh.u(i,j)=ux; mesh.v(i,j)=vy;
            mesh.p(i,j)=(g-1.0)*(mesh.U3(i,j)-0.5*rho*(ux*ux+vy*vy));
        }
}

double recoverAndMaxSpeed(Mesh& mesh)
{
    const int ny=mesh.ny, nx=mesh.nx;
    const double g=mesh.gamma;
    double umax=0.0;
    for (int j = 0; j < nx; ++j)
        for (int i = 0; i < ny; ++i) {
            double rho=mesh.U0(i,j), ux=mesh.U1(i,j)/rho, vy=mesh.U2(i,j)/rho;
            double pr=(g-1.0)*(mesh.U3(i,j)-0.5*rho*(ux*ux+vy*vy));
            mesh.rho(i,j)=rho; mesh.u(i,j)=ux; mesh.v(i,j)=vy; mesh.p(i,j)=pr;
            double a=sqrt(g*pr/rho);
            double spd=sqrt(ux*ux+vy*vy)+a;
            if (spd > umax) umax = spd;
        }
    return umax;
}

double computeMaxSpeed(const Mesh& mesh, double gamma)
{
    double umax=0.0;
    const int ny=mesh.ny, nx=mesh.nx;
    for (int j = 0; j < nx; ++j)
        for (int i = 0; i < ny; ++i) {
            double rho=mesh.U0(i,j), u=mesh.U1(i,j)/rho, v=mesh.U2(i,j)/rho;
            double p=(gamma-1.0)*(mesh.U3(i,j)-0.5*rho*(u*u+v*v));
            double a=sqrt(gamma*p/rho);
            double spd=sqrt(u*u+v*v)+a;
            if (spd > umax) umax = spd;
        }
    return umax;
}

// ============================================================================
// saveMeshData（与原版相同）
// ============================================================================
void saveMeshData(const Mesh& mesh, int rank, const std::string& timestep_folder)
{
    try {
        fs::path dir;
        if (!timestep_folder.empty()) { dir=fs::path(timestep_folder); fs::create_directories(dir); }
        auto write=[&](const std::string& name, const auto& field){
            fs::path p=dir.empty()
                ? fs::path(name+"_"+std::to_string(rank)+".dat")
                : dir/(name+"_"+std::to_string(rank)+".dat");
            std::ofstream f(p);
            if (!f) throw std::runtime_error("无法创建文件: "+p.string());
            f<<field;
        };
        write("U0",mesh.U0); write("U1",mesh.U1);
        write("U2",mesh.U2); write("U3",mesh.U3);
    } catch(const std::exception& e){
        std::cerr<<"[Rank "<<rank<<"] 保存 Mesh 数据失败: "<<e.what()<<std::endl;
        throw;
    }
}