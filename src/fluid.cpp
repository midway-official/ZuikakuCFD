#include "fluid.h"

namespace fs = std::filesystem;
// 全局变量定义


// ============================================================================
// 工具函数：读取 double 矩阵
// ============================================================================
static void readMatrix(const std::string& filename, MatrixXd& mat)
{
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    for (int i = 0; i < mat.rows(); ++i)
        for (int j = 0; j < mat.cols(); ++j)
            file >> mat(i, j);
}

// ============================================================================
// 工具函数：读取 int 矩阵
// ============================================================================
static void readMatrixInt(const std::string& filename, MatrixXi& mat)
{
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    for (int i = 0; i < mat.rows(); ++i)
        for (int j = 0; j < mat.cols(); ++j)
            file >> mat(i, j);
}


// ============================================================================
// 构造函数1：尺寸初始化（全0场）
// ============================================================================
Mesh::Mesh(int n_y, int n_x, double da_, double gamma_)
{
    ny = n_y;
    nx = n_x;
    da = da_;
    gamma = gamma_;

    // 物理量
    rho = MatrixXd::Zero(ny, nx);
    u   = MatrixXd::Zero(ny, nx);
    v   = MatrixXd::Zero(ny, nx);
    p   = MatrixXd::Zero(ny, nx);

    // 守恒变量
    U0  = MatrixXd::Zero(ny, nx);
    U1  = MatrixXd::Zero(ny, nx);
    U2  = MatrixXd::Zero(ny, nx);
    U3  = MatrixXd::Zero(ny, nx);

    // 边界类型
    bctype = MatrixXi::Zero(ny, nx);
}


// ============================================================================
// 构造函数2：从文件夹读取网格
// ============================================================================
Mesh::Mesh(const std::string& folderPath)
{
    // -------------------------
    // 读取 params.txt
    // -------------------------
    std::ifstream paramFile(folderPath + "/params.txt");

    if (!paramFile) {
        throw std::runtime_error("Cannot open params.txt");
    }

    paramFile >> nx >> ny >> da >> gamma;

    // -------------------------
    // 初始化矩阵
    // -------------------------
    rho = MatrixXd(ny, nx);
    u   = MatrixXd(ny, nx);
    v   = MatrixXd(ny, nx);
    p   = MatrixXd(ny, nx);

    U0  = MatrixXd(ny, nx);
    U1  = MatrixXd(ny, nx);
    U2  = MatrixXd(ny, nx);
    U3  = MatrixXd(ny, nx);

    bctype = MatrixXi(ny, nx);

    // -------------------------
    // 读取数据文件
    // -------------------------
    readMatrix(folderPath + "/rho.dat", rho);
    readMatrix(folderPath + "/u.dat",   u);
    readMatrix(folderPath + "/v.dat",   v);
    readMatrix(folderPath + "/p.dat",   p);

    readMatrixInt(folderPath + "/bctype.dat", bctype);

    // -------------------------
    // 计算守恒变量
    // -------------------------
    for (int i = 0; i < ny; ++i)
    {
        for (int j = 0; j < nx; ++j)
        {
            double r = rho(i,j);
            double ux = u(i,j);
            double vy = v(i,j);
            double pressure = p(i,j);

            U0(i,j) = r;
            U1(i,j) = r * ux;
            U2(i,j) = r * vy;

            double kinetic = 0.5 * r * (ux*ux + vy*vy);
            double internal = pressure / (gamma - 1.0);

            U3(i,j) = internal + kinetic;
        }
    }
}

vector<Mesh> splitMeshVertically(const Mesh& original, int n)
{
    vector<Mesh> sub_meshes;

    const int Nx = original.nx;
    const int Ny = original.ny;

    // ===== 计算每段宽度 =====
    vector<int> widths(n);
    int remain = Nx;
    for(int k = 0; k < n; ++k)
    {
        widths[k] = remain / (n - k);
        remain -= widths[k];
    }

    int start = 0;

    for(int k = 0; k < n; ++k)
    {
        int real_w = widths[k];

        bool has_left  = (k > 0);
        bool has_right = (k < n-1);

        int left_ghost  = has_left  ? 2 : 0;
        int right_ghost = has_right ? 2 : 0;

        int sub_nx = real_w + left_ghost + right_ghost;

        Mesh sub(Ny, sub_nx, original.da, original.gamma);

        int real_start = left_ghost;

        // =====================================================
        // 复制真实数据
        // =====================================================
        for(int i = 0; i < Ny; ++i)
        {
            for(int j = 0; j < real_w; ++j)
            {
                int orig_j = start + j;
                int sub_j  = real_start + j;

                // ---- primitive variables ----
                sub.rho(i,sub_j) = original.rho(i,orig_j);
                sub.u  (i,sub_j) = original.u  (i,orig_j);
                sub.v  (i,sub_j) = original.v  (i,orig_j);
                sub.p  (i,sub_j) = original.p  (i,orig_j);

                // ---- conservative variables ----
                sub.U0(i,sub_j) = original.U0(i,orig_j);
                sub.U1(i,sub_j) = original.U1(i,orig_j);
                sub.U2(i,sub_j) = original.U2(i,orig_j);
                sub.U3(i,sub_j) = original.U3(i,orig_j);

                // ---- boundary type ----
                sub.bctype(i,sub_j) = original.bctype(i,orig_j);
            }
        }

        // =====================================================
        // 填充 ghost 区域
        // =====================================================
        for(int i = 0; i < Ny; ++i)
        {
            // 左 ghost
            if(has_left)
            {
                for(int gj = 0; gj < left_ghost; ++gj)
                {
                    int orig_j = start - left_ghost + gj; // 左邻居的真实列
                    int sub_j = gj;                        // ghost列
                    // 防止越界
                    if(orig_j < 0) {
                        orig_j = 0;
                    }

                    sub.rho(i,sub_j) = original.rho(i,orig_j);
                    sub.u  (i,sub_j) = original.u  (i,orig_j);
                    sub.v  (i,sub_j) = original.v  (i,orig_j);
                    sub.p  (i,sub_j) = original.p  (i,orig_j);

                    sub.U0(i,sub_j) = original.U0(i,orig_j);
                    sub.U1(i,sub_j) = original.U1(i,orig_j);
                    sub.U2(i,sub_j) = original.U2(i,orig_j);
                    sub.U3(i,sub_j) = original.U3(i,orig_j);

                    sub.bctype(i,sub_j) = -3;
                }
            }

            // 右 ghost
            if(has_right)
            {
                for(int gj = 0; gj < right_ghost; ++gj)
                {
                    int orig_j = start + real_w + gj;   // 右邻居的真实列
                    int sub_j = real_start + real_w + gj; // ghost列
                    if(orig_j >= Nx) {
                        orig_j = Nx-1;
                    }

                    sub.rho(i,sub_j) = original.rho(i,orig_j);
                    sub.u  (i,sub_j) = original.u  (i,orig_j);
                    sub.v  (i,sub_j) = original.v  (i,orig_j);
                    sub.p  (i,sub_j) = original.p  (i,orig_j);

                    sub.U0(i,sub_j) = original.U0(i,orig_j);
                    sub.U1(i,sub_j) = original.U1(i,orig_j);
                    sub.U2(i,sub_j) = original.U2(i,orig_j);
                    sub.U3(i,sub_j) = original.U3(i,orig_j);

                    sub.bctype(i,sub_j) = -3;
                }
            }
        }

        sub_meshes.push_back(sub);

        start += real_w;
    }

    return sub_meshes;
}


//超高性能的列交换函数
void exchangeColumns(MatrixXd& matrix, int rank, int num_procs) {
    const int rows = matrix.rows();
    const int cols = matrix.cols();
    const int count = rows * 2; // 每次交换2列

    // 确定邻居
    int left_rank  = (rank == 0) ? MPI_PROC_NULL : rank - 1;
    int right_rank = (rank == num_procs - 1) ? MPI_PROC_NULL : rank + 1;

    // 分配缓冲区：Eigen 默认是 ColMajor，列内数据连续
    VectorXd send_left = Map<VectorXd>(matrix.block(0, 2, rows, 2).data(), count);
    VectorXd send_right = Map<VectorXd>(matrix.block(0, cols - 4, rows, 2).data(), count);
    VectorXd recv_left(count), recv_right(count);

    // 使用 Sendrecv 
    // 向左发，从左收
    MPI_Sendrecv(send_left.data(),  count, MPI_DOUBLE, left_rank,  0,
                 recv_left.data(),  count, MPI_DOUBLE, left_rank,  1, 
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // 向右发，从右收
    MPI_Sendrecv(send_right.data(), count, MPI_DOUBLE, right_rank, 1,
                 recv_right.data(), count, MPI_DOUBLE, right_rank, 0, 
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // 写回数据
    if (left_rank != MPI_PROC_NULL)
        matrix.block(0, 0, rows, 2) = Map<MatrixXd>(recv_left.data(), rows, 2);
    
    if (right_rank != MPI_PROC_NULL)
        matrix.block(0, cols - 2, rows, 2) = Map<MatrixXd>(recv_right.data(), rows, 2);
}


void muscl_reconstruct(
    double UL2, double UL1, double UR1, double UR2,
    int limiter,
    double& UL, double& UR)
{
    // ====================================================
    // 斜率计算
    // ====================================================
    double dL = UL1 - UL2;   // 左侧向后差分
    double dC = UR1 - UL1;   // 中心差分（左右共用）
    double dR = UR2 - UR1;   // 右侧向前差分

    double sigma_L = 0.0;    // 左状态斜率
    double sigma_R = 0.0;    // 右状态斜率

    // ====================================================
    // 限制器
    // ====================================================
    auto minmod2 = [](double a, double b) -> double {
        if (a * b <= 0.0) return 0.0;
        return std::abs(a) < std::abs(b) ? a : b;
    };

    auto van_leer = [](double a, double b) -> double {
        double ab = a * b;
        if (ab <= 0.0) return 0.0;
        return 2.0 * ab / (a + b);
    };

    auto superbee = [&](double a, double b) -> double {
        // max(minmod(a, 2b), minmod(2a, b))
        double s1 = minmod2(a, 2.0 * b);
        double s2 = minmod2(2.0 * a, b);
        // maxmod: 同号取绝对值更大的
        if (s1 * s2 <= 0.0) return 0.0;
        return std::abs(s1) > std::abs(s2) ? s1 : s2;
    };

    switch (limiter)
    {
        case 0: // minmod
            sigma_L = minmod2(dL, dC);
            sigma_R = minmod2(dC, dR);
            break;

        case 1: // van leer
            sigma_L = van_leer(dL, dC);
            sigma_R = van_leer(dC, dR);
            break;

        case 2: // superbee
            sigma_L = superbee(dL, dC);
            sigma_R = superbee(dC, dR);
            break;

        default:
            sigma_L = 0.0;   // 退化为一阶
            sigma_R = 0.0;
            break;
    }

    // ====================================================
    // 界面重构
    //   左状态：UL1 向右半个格子外推
    //   右状态：UR1 向左半个格子外推
    // ====================================================
    UL = UL1 + 0.5 * sigma_L;
    UR = UR1 - 0.5 * sigma_R;
}


void exchangeConservativeColumns(Mesh& mesh, int rank, int num_procs)
{
    exchangeColumns(mesh.U0, rank, num_procs);
    exchangeColumns(mesh.U1, rank, num_procs);
    exchangeColumns(mesh.U2, rank, num_procs);
    exchangeColumns(mesh.U3, rank, num_procs);
}
void hllcFlux(
    double UL0, double UL1, double UL2, double UL3,
    double UR0, double UR1, double UR2, double UR3,
    double gamma,
    double& F0, double& F1, double& F2, double& F3)
{
    const double eps = 1e-12;

    double rhoL = std::max(UL0, eps);
    double uL   = UL1 / rhoL;
    double vL   = UL2 / rhoL;
    double EL   = UL3 / rhoL;
    double pL   = std::max((gamma-1.0)*(UL3 - 0.5*rhoL*(uL*uL+vL*vL)), eps);
    double aL   = sqrt(gamma * pL / rhoL);

    double rhoR = std::max(UR0, eps);
    double uR   = UR1 / rhoR;
    double vR   = UR2 / rhoR;
    double ER   = UR3 / rhoR;
    double pR   = std::max((gamma-1.0)*(UR3 - 0.5*rhoR*(uR*uR+vR*vR)), eps);
    double aR   = sqrt(gamma * pR / rhoR);

    // 物理通量
    double FL0 = rhoL*uL,  FL1 = rhoL*uL*uL+pL, FL2 = rhoL*uL*vL, FL3 = (UL3+pL)*uL;
    double FR0 = rhoR*uR,  FR1 = rhoR*uR*uR+pR, FR2 = rhoR*uR*vR, FR3 = (UR3+pR)*uR;

    // Roe 平均波速
    double sqL = sqrt(rhoL), sqR = sqrt(rhoR), den = sqL + sqR;
    double u_roe = (sqL*uL + sqR*uR) / den;
    double h_roe = (sqL*(EL+pL/rhoL) + sqR*(ER+pR/rhoR)) / den;
    double c_roe = sqrt(std::max((gamma-1.0)*(h_roe - 0.5*u_roe*u_roe), eps));

    double SL = std::min(uL - aL, u_roe - c_roe);
    double SR = std::max(uR + aR, u_roe + c_roe);

    if (SL >= 0.0) {
        F0=FL0; F1=FL1; F2=FL2; F3=FL3; return;
    }
    if (SR <= 0.0) {
        F0=FR0; F1=FR1; F2=FR2; F3=FR3; return;
    }

    // ★ 接触波速 S*
    double num = pR - pL + rhoL*uL*(SL-uL) - rhoR*uR*(SR-uR);
    double den2 = rhoL*(SL-uL) - rhoR*(SR-uR);
    double Sstar = num / den2;

    // ★ 构造中间状态 U*_K = rho_K * (S_K - u_K)/(S_K - S*) * [...]
    auto star_state = [&](double rho, double u, double v, double E,
                          double p, double S) -> std::array<double,4>
    {
        double coeff = rho * (S - u) / (S - Sstar);
        return {
            coeff,
            coeff * Sstar,
            coeff * v,
            coeff * (E + (Sstar - u)*(Sstar + p/(rho*(S-u))))
        };
    };

    if (Sstar >= 0.0) {
        // F* = FL + SL*(U*L - UL)
        auto UstarL = star_state(rhoL, uL, vL, EL, pL, SL);
        F0 = FL0 + SL*(UstarL[0] - UL0);
        F1 = FL1 + SL*(UstarL[1] - UL1);
        F2 = FL2 + SL*(UstarL[2] - UL2);
        F3 = FL3 + SL*(UstarL[3] - UL3);
    } else {
        // F* = FR + SR*(U*R - UR)
        auto UstarR = star_state(rhoR, uR, vR, ER, pR, SR);
        F0 = FR0 + SR*(UstarR[0] - UR0);
        F1 = FR1 + SR*(UstarR[1] - UR1);
        F2 = FR2 + SR*(UstarR[2] - UR2);
        F3 = FR3 + SR*(UstarR[3] - UR3);
    }
}
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
    int limiter )
{
    // --------------------------------------------------------
    // MUSCL 重构：对每个分量独立重构，得到界面左右状态
    //
    //  x方向右界面 j+1/2：模板 [ Ul1 | Up | Ur1 | Ur2 ]
    //  x方向左界面 j-1/2：模板 [ Ul2 | Ul1 | Up | Ur1 ]
    //  y方向下界面 i+1/2：模板 [ Uu1 | Up | Ud1 | Ud2 ]
    //  y方向上界面 i-1/2：模板 [ Uu2 | Uu1 | Up | Ud1 ]
    // --------------------------------------------------------
    double QL_right[4], QR_right[4];   // 右界面 j+1/2
    double QL_left [4], QR_left [4];   // 左界面 j-1/2
    double QL_down [4], QR_down [4];   // 下界面 i+1/2
    double QL_up   [4], QR_up   [4];   // 上界面 i-1/2

    for (int c = 0; c < 4; ++c)
    {
        // x 右界面：UL2=Ul1, UL1=Up, UR1=Ur1, UR2=Ur2
        muscl_reconstruct(Ul1[c], Up[c], Ur1[c], Ur2[c],
                          limiter, QL_right[c], QR_right[c]);

        // x 左界面：UL2=Ul2, UL1=Ul1, UR1=Up, UR2=Ur1
        muscl_reconstruct(Ul2[c], Ul1[c], Up[c], Ur1[c],
                          limiter, QL_left[c], QR_left[c]);

        // y 下界面：UL2=Uu1, UL1=Up, UR1=Ud1, UR2=Ud2
        muscl_reconstruct(Uu1[c], Up[c], Ud1[c], Ud2[c],
                          limiter, QL_down[c], QR_down[c]);

        // y 上界面：UL2=Uu2, UL1=Uu1, UR1=Up, UR2=Ud1
        muscl_reconstruct(Uu2[c], Uu1[c], Up[c], Ud1[c],
                          limiter, QL_up[c], QR_up[c]);
    }

    // --------------------------------------------------------
    // 计算各界面数值通量
    // --------------------------------------------------------
    double F_right[4], F_left[4], F_down[4], F_up[4];

    // x 方向：直接用重构后的 QL/QR
    hllcFlux(QL_right[0], QL_right[1], QL_right[2], QL_right[3],
             QR_right[0], QR_right[1], QR_right[2], QR_right[3],
             gamma, F_right[0], F_right[1], F_right[2], F_right[3]);

    hllcFlux(QL_left[0],  QL_left[1],  QL_left[2],  QL_left[3],
             QR_left[0],  QR_left[1],  QR_left[2],  QR_left[3],
             gamma, F_left[0],  F_left[1],  F_left[2],  F_left[3]);

    // y 方向：交换 U1/U2（法向变切向）再传入
    hllcFlux(QL_down[0], QL_down[2], QL_down[1], QL_down[3],
             QR_down[0], QR_down[2], QR_down[1], QR_down[3],
             gamma, F_down[0], F_down[1], F_down[2], F_down[3]);

    hllcFlux(QL_up[0],   QL_up[2],   QL_up[1],   QL_up[3],
             QR_up[0],   QR_up[2],   QR_up[1],   QR_up[3],
             gamma, F_up[0],   F_up[1],   F_up[2],   F_up[3]);

    // --------------------------------------------------------
    // 有限体积更新
    // --------------------------------------------------------
    std::vector<double> U_new(4);

    U_new[0] = Up[0] - (dt/dx)*(F_right[0] - F_left[0])
                     - (dt/dy)*(F_down[0]   - F_up[0]);

    U_new[1] = Up[1] - (dt/dx)*(F_right[1] - F_left[1])
                     - (dt/dy)*(F_down[2]   - F_up[2]);   // 切向

    U_new[2] = Up[2] - (dt/dx)*(F_right[2] - F_left[2])
                     - (dt/dy)*(F_down[1]   - F_up[1]);   // 法向

    U_new[3] = Up[3] - (dt/dx)*(F_right[3] - F_left[3])
                     - (dt/dy)*(F_down[3]   - F_up[3]);

    return U_new;
}
// ============================================================================
// 辅助函数：计算单步残差 dU（即 dt * L(U)）
// ============================================================================
static void computeRHS(
    const Mesh& mesh,
    double dt,
    MatrixXd& dU0, MatrixXd& dU1, MatrixXd& dU2, MatrixXd& dU3)
{
    const int ny   = mesh.ny;
    const int nx   = mesh.nx;
    const double dx = mesh.da;
    const double dy = mesh.da;
    const double gamma = mesh.gamma;

    dU0 = MatrixXd::Zero(ny, nx);
    dU1 = MatrixXd::Zero(ny, nx);
    dU2 = MatrixXd::Zero(ny, nx);
    dU3 = MatrixXd::Zero(ny, nx);

    for (int i = 2; i < ny - 2; ++i)
    {
        for (int j = 2; j < nx - 2; ++j)
        {
            if (mesh.bctype(i, j) != 0) continue;

            // ---- 取当前格 Up ----
            std::vector<double> Up = {
                mesh.U0(i,j), mesh.U1(i,j), mesh.U2(i,j), mesh.U3(i,j)
            };

            // ---- 邻格取值（与原 updateMesh 完全相同的逻辑） ----
            auto getNeighbor = [&](int ii, int jj) -> std::vector<double> {
                int bc = mesh.bctype(ii, jj);
                if (bc == 0 || bc == -3)
                    return {mesh.U0(ii,jj), mesh.U1(ii,jj),
                            mesh.U2(ii,jj), mesh.U3(ii,jj)};
                else if (bc == -1)
                    return Up;   // 零梯度
                else
                    return {mesh.U0(ii,jj), mesh.U1(ii,jj),
                            mesh.U2(ii,jj), mesh.U3(ii,jj)};
            };

            auto Ur1 = getNeighbor(i, j+1);
            auto Ur2 = getNeighbor(i, j+2);
            auto Ul1 = getNeighbor(i, j-1);
            auto Ul2 = getNeighbor(i, j-2);
            auto Uu1 = getNeighbor(i-1, j);
            auto Uu2 = getNeighbor(i-2, j);
            auto Ud1 = getNeighbor(i+1, j);
            auto Ud2 = getNeighbor(i+2, j);

            // ---- 计算更新量 ----
            std::vector<double> U_new = updateCenterCell(
                Up,
                Ur1, Ur2, Ul1, Ul2,
                Uu1, Uu2, Ud1, Ud2,
                gamma, dt, dx, dy, 2 /* limiter */);

            // dU = U_new - U_old（即 dt * L(U)）
            dU0(i,j) = U_new[0] - Up[0];
            dU1(i,j) = U_new[1] - Up[1];
            dU2(i,j) = U_new[2] - Up[2];
            dU3(i,j) = U_new[3] - Up[3];
        }
    }
}

// ============================================================================
// 2阶 SSP-RK2 时间推进
//   Stage 1: U*       = U^n + dt·L(U^n)
//   Stage 2: U^{n+1}  = 1/2·U^n + 1/2·(U* + dt·L(U*))
// ============================================================================
void updateMesh(Mesh& mesh, double dt, int rank, int num_procs)
{
    // ---------- Stage 1 ----------
    MatrixXd dU0, dU1, dU2, dU3;
    computeRHS(mesh, dt, dU0, dU1, dU2, dU3);

    // 保存 U^n
    MatrixXd U0n = mesh.U0, U1n = mesh.U1,
             U2n = mesh.U2, U3n = mesh.U3;

    // U* = U^n + dU
    mesh.U0 = U0n + dU0;
    mesh.U1 = U1n + dU1;
    mesh.U2 = U2n + dU2;
    mesh.U3 = U3n + dU3;
    exchangeConservativeColumns(mesh, rank, num_procs); // 交换 U* 的边界列
    // ---------- Stage 2 ----------
    MatrixXd dU0s, dU1s, dU2s, dU3s;
    computeRHS(mesh, dt, dU0s, dU1s, dU2s, dU3s);

    // U^{n+1} = 0.5·U^n + 0.5·(U* + dU*)
    mesh.U0 = 0.5 * U0n + 0.5 * (mesh.U0 + dU0s);
    mesh.U1 = 0.5 * U1n + 0.5 * (mesh.U1 + dU1s);
    mesh.U2 = 0.5 * U2n + 0.5 * (mesh.U2 + dU2s);
    mesh.U3 = 0.5 * U3n + 0.5 * (mesh.U3 + dU3s);
}

// ============================================================================
// 从守恒变量恢复原始变量（密度、速度、压力）
// ============================================================================
void recoverPrimitives(Mesh& mesh)
{
    const int ny = mesh.ny;
    const int nx = mesh.nx;
    const double gamma = mesh.gamma;

    for(int i = 0; i < ny; ++i)
    {
        for(int j = 0; j < nx; ++j)
        {
            double rho = mesh.U0(i,j);
            mesh.rho(i,j) = rho;
            mesh.u(i,j) = mesh.U1(i,j) / rho;
            mesh.v(i,j) = mesh.U2(i,j) / rho;

            double kinetic = 0.5 * rho * (mesh.u(i,j)*mesh.u(i,j) + mesh.v(i,j)*mesh.v(i,j));
            double internal = mesh.U3(i,j) - kinetic;
            mesh.p(i,j) = (gamma - 1.0) * internal;
        }
    }
}

void saveMeshData(
    const Mesh& mesh,
    int rank,
    const std::string& timestep_folder)
{
    try {
        fs::path dir;

        if (!timestep_folder.empty()) {
            dir = fs::path(timestep_folder);
            fs::create_directories(dir);   // MPI-safe 幂等
        }

        auto write = [&](const std::string& name, const auto& field) {
            fs::path p = dir.empty()
                ? fs::path(name + "_" + std::to_string(rank) + ".dat")
                : dir / (name + "_" + std::to_string(rank) + ".dat");

            std::ofstream f(p);
            if (!f) {
                throw std::runtime_error("无法创建文件: " + p.string());
            }
            f << field;
        };

        write("U0", mesh.U0);
        write("U1", mesh.U1);
        write("U2", mesh.U2);
        write("U3", mesh.U3);

    } catch (const std::exception& e) {
        std::cerr << "[Rank " << rank << "] 保存 Mesh 数据失败: "
                  << e.what() << std::endl;
        throw;   // 关键 
    }
}


// 在主函数中添加CFL检查
double computeMaxSpeed(const Mesh& mesh, double gamma)
{
    double umax = 0.0;
    for(int i = 0; i < mesh.ny; ++i)
    {
        for(int j = 0; j < mesh.nx; ++j)
        {
            double rho = mesh.U0(i,j);
            double u = mesh.U1(i,j) / rho;
            double v = mesh.U2(i,j) / rho;
            double p = (gamma - 1.0) * (mesh.U3(i,j) - 0.5*rho*(u*u + v*v));
            double a = sqrt(gamma * p / rho);
            double speed = sqrt(u*u + v*v) + a;
            umax = std::max(umax, speed);
        }
    }
    return umax;
}
