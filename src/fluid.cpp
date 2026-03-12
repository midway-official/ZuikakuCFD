#include "fluid.h"

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

        int left_ghost  = has_left  ? 3 : 0;
        int right_ghost = has_right ? 3 : 0;

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
    const int count = rows * 3; // 每次交换3列

    // 确定邻居
    int left_rank  = (rank == 0) ? MPI_PROC_NULL : rank - 1;
    int right_rank = (rank == num_procs - 1) ? MPI_PROC_NULL : rank + 1;

    // 分配缓冲区：Eigen 默认是 ColMajor，列内数据连续
    VectorXd send_left = Map<VectorXd>(matrix.block(0, 3, rows, 3).data(), count);
    VectorXd send_right = Map<VectorXd>(matrix.block(0, cols - 6, rows, 3).data(), count);
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
        matrix.block(0, 0, rows, 3) = Map<MatrixXd>(recv_left.data(), rows, 3);
    
    if (right_rank != MPI_PROC_NULL)
        matrix.block(0, cols - 3, rows, 3) = Map<MatrixXd>(recv_right.data(), rows, 3);
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

/**
 * WENO5 (5th-order Weighted Essentially Non-Oscillatory) 重构
 * 重构 UP 单元右界面处的左右状态
 *
 * 输入点索引示意:
 *   UL2   UL1   UP   UR1   UR2   UR3
 *   i-2   i-1   i   i+1   i+2   i+3
 *                    ^
 *                重构此界面 (i+1/2)
 *
 * @param UL2~UR3  模板点值
 * @param UL       输出：界面左状态 (来自左侧的重构)
 * @param UR       输出：界面右状态 (来自右侧的重构)
 */
inline void weno5_reconstruct(
    double UL2, double UL1, double UP,
    double UR1, double UR2, double UR3,
    double& UL, double& UR)
{
    constexpr double eps = 1e-6;

    constexpr double C1 = 13.0/12.0;
    constexpr double C2 = 0.25;

    // =====================================================
    // 左状态 UL
    // =====================================================

    // candidate polynomials
    double qL0 =  (1.0/3.0)*UL2 - (7.0/6.0)*UL1 + (11.0/6.0)*UP;
    double qL1 = -(1.0/6.0)*UL1 + (5.0/6.0)*UP  + (1.0/3.0)*UR1;
    double qL2 =  (1.0/3.0)*UP  + (5.0/6.0)*UR1 - (1.0/6.0)*UR2;

    // smoothness indicators
    double d0 = UL2 - 2.0*UL1 + UP;
    double d1 = UL2 - 4.0*UL1 + 3.0*UP;

    double d2 = UL1 - 2.0*UP + UR1;
    double d3 = UL1 - UR1;

    double d4 = UP - 2.0*UR1 + UR2;
    double d5 = 3.0*UP - 4.0*UR1 + UR2;

    double bL0 = C1*(d0*d0) + C2*(d1*d1);
    double bL1 = C1*(d2*d2) + C2*(d3*d3);
    double bL2 = C1*(d4*d4) + C2*(d5*d5);

    // nonlinear weights
    double t0 = eps + bL0;
    double t1 = eps + bL1;
    double t2 = eps + bL2;

    double aL0 = 0.1 / (t0*t0);
    double aL1 = 0.6 / (t1*t1);
    double aL2 = 0.3 / (t2*t2);

    double aLsum = aL0 + aL1 + aL2;

    double wL0 = aL0 / aLsum;
    double wL1 = aL1 / aLsum;
    double wL2 = aL2 / aLsum;

    UL = wL0*qL0 + wL1*qL1 + wL2*qL2;

    // =====================================================
    // 右状态 UR
    // =====================================================

    double qR0 = (11.0/6.0)*UR1 - (7.0/6.0)*UR2 + (1.0/3.0)*UR3;
    double qR1 = (1.0/3.0)*UP  + (5.0/6.0)*UR1 - (1.0/6.0)*UR2;
    double qR2 =-(1.0/6.0)*UL1 + (5.0/6.0)*UP  + (1.0/3.0)*UR1;

    double r0 = UR1 - 2.0*UR2 + UR3;
    double r1 = 3.0*UR1 - 4.0*UR2 + UR3;

    double r2 = UP - 2.0*UR1 + UR2;
    double r3 = UP - UR2;

    double r4 = UL1 - 2.0*UP + UR1;
    double r5 = UL1 - 4.0*UP + 3.0*UR1;

    double bR0 = C1*(r0*r0) + C2*(r1*r1);
    double bR1 = C1*(r2*r2) + C2*(r3*r3);
    double bR2 = C1*(r4*r4) + C2*(r5*r5);

    double s0 = eps + bR0;
    double s1 = eps + bR1;
    double s2 = eps + bR2;

    double aR0 = 0.1 / (s0*s0);
    double aR1 = 0.6 / (s1*s1);
    double aR2 = 0.3 / (s2*s2);

    double aRsum = aR0 + aR1 + aR2;

    double wR0 = aR0 / aRsum;
    double wR1 = aR1 / aRsum;
    double wR2 = aR2 / aRsum;

    UR = wR0*qR0 + wR1*qR1 + wR2*qR2;
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
    if (fabs(den2) < 1e-12)
    den2 = (den2>=0 ? 1e-12 : -1e-12);
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
    const std::vector<double>& Ur3,   // 新增：x方向 i+3
    const std::vector<double>& Ul1,
    const std::vector<double>& Ul2,
    const std::vector<double>& Ul3,   // 新增：x方向 i-3
    const std::vector<double>& Uu1,
    const std::vector<double>& Uu2,
    const std::vector<double>& Uu3,   // 新增：y方向 j-3
    const std::vector<double>& Ud1,
    const std::vector<double>& Ud2,
    const std::vector<double>& Ud3,   // 新增：y方向 j+3
    double gamma,
    double dt,
    double dx,
    double dy )
{
    // --------------------------------------------------------
    // WENO5 重构：对每个分量独立重构，得到界面左右状态
    //
    //  x方向右界面 j+1/2：模板 [ Ul2 | Ul1 | Up | Ur1 | Ur2 | Ur3 ]
    //  x方向左界面 j-1/2：模板 [ Ul3 | Ul2 | Ul1 | Up | Ur1 | Ur2 ]
    //  y方向下界面 i+1/2：模板 [ Uu2 | Uu1 | Up | Ud1 | Ud2 | Ud3 ]
    //  y方向上界面 i-1/2：模板 [ Uu3 | Uu2 | Uu1 | Up | Ud1 | Ud2 ]
    // --------------------------------------------------------
    double QL_right[4], QR_right[4];   // 右界面 j+1/2
    double QL_left [4], QR_left [4];   // 左界面 j-1/2
    double QL_down [4], QR_down [4];   // 下界面 i+1/2
    double QL_up   [4], QR_up   [4];   // 上界面 i-1/2

    for (int c = 0; c < 4; ++c)
    {
        // x 右界面 j+1/2：中心 = Up
        //   UL2=Ul2, UL1=Ul1, UP=Up, UR1=Ur1, UR2=Ur2, UR3=Ur3
        weno5_reconstruct(Ul2[c], Ul1[c], Up[c], Ur1[c], Ur2[c], Ur3[c],
                          QL_right[c], QR_right[c]);

        // x 左界面 j-1/2：中心 = Ul1
        //   UL2=Ul3, UL1=Ul2, UP=Ul1, UR1=Up, UR2=Ur1, UR3=Ur2
        weno5_reconstruct(Ul3[c], Ul2[c], Ul1[c], Up[c], Ur1[c], Ur2[c],
                          QL_left[c], QR_left[c]);

        // y 下界面 i+1/2：中心 = Up（沿 y 轴正方向）
        //   UL2=Uu2, UL1=Uu1, UP=Up, UR1=Ud1, UR2=Ud2, UR3=Ud3
        weno5_reconstruct(Uu2[c], Uu1[c], Up[c], Ud1[c], Ud2[c], Ud3[c],
                          QL_down[c], QR_down[c]);

        // y 上界面 i-1/2：中心 = Uu1（沿 y 轴负方向）
        //   UL2=Uu3, UL1=Uu2, UP=Uu1, UR1=Up, UR2=Ud1, UR3=Ud2
        weno5_reconstruct(Uu3[c], Uu2[c], Uu1[c], Up[c], Ud1[c], Ud2[c],
                          QL_up[c], QR_up[c]);
    }

    // --------------------------------------------------------
    // 计算各界面数值通量
    // --------------------------------------------------------
    double F_right[4], F_left[4], F_down[4], F_up[4];

    // x 方向通量（法向为 x）
    hllcFlux(QL_right[0], QL_right[1], QL_right[2], QL_right[3],
             QR_right[0], QR_right[1], QR_right[2], QR_right[3],
             gamma, F_right[0], F_right[1], F_right[2], F_right[3]);

    hllcFlux(QL_left[0],  QL_left[1],  QL_left[2],  QL_left[3],
             QR_left[0],  QR_left[1],  QR_left[2],  QR_left[3],
             gamma, F_left[0],  F_left[1],  F_left[2],  F_left[3]);

    // y 方向通量：交换 [1]/[2]（法向 y ↔ 切向 x）后传入 hllcFlux
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
                     - (dt/dy)*(F_down[2]   - F_up[2]);   // y通量中切向分量

    U_new[2] = Up[2] - (dt/dx)*(F_right[2] - F_left[2])
                     - (dt/dy)*(F_down[1]   - F_up[1]);   // y通量中法向分量

    U_new[3] = Up[3] - (dt/dx)*(F_right[3] - F_left[3])
                     - (dt/dy)*(F_down[3]   - F_up[3]);

    return U_new;
}

void computeRHS(
    Mesh& mesh,
    double dt,
    MatrixXd& dU0, MatrixXd& dU1, MatrixXd& dU2, MatrixXd& dU3)
{
    const int ny    = mesh.ny;
    const int nx    = mesh.nx;
    const double dx = mesh.da;
    const double dy = mesh.da;
    const double gamma = mesh.gamma;

    // ============================================================
    // 复制 Ghost Cell：遍历内部 bctype==0 的中心点 (i,j)
    // 检查其上下左右 ±1、±2 的邻格，若 bctype==-1 则将中心点守恒量写入
    // ============================================================

    constexpr int offsets[6] = {-3,-2, -1, 1, 2, 3};

    for (int i = 3; i < ny - 3; ++i)
    {
        for (int j = 3; j < nx - 3; ++j)
        {
            

            // 行方向（上下）邻格
            for (int dk : offsets)
            {
                int ni = i + dk;
                if (ni < 0 || ni >= ny) continue;
                if (mesh.bctype(ni, j) == -1)
                {
                    mesh.U0(ni, j) = mesh.U0(i, j);
                    mesh.U1(ni, j) = mesh.U1(i, j);
                    mesh.U2(ni, j) = mesh.U2(i, j);
                    mesh.U3(ni, j) = mesh.U3(i, j);
                }
            }

            // 列方向（左右）邻格
            for (int dk : offsets)
            {
                int nj = j + dk;
                if (nj < 0 || nj >= nx) continue;
                if (mesh.bctype(i, nj) == -1)
                {
                    mesh.U0(i, nj) = mesh.U0(i, j);
                    mesh.U1(i, nj) = mesh.U1(i, j);
                    mesh.U2(i, nj) = mesh.U2(i, j);
                    mesh.U3(i, nj) = mesh.U3(i, j);
                }
            }
        }
    }

    // ============================================================

    dU0 = MatrixXd::Zero(ny, nx);
    dU1 = MatrixXd::Zero(ny, nx);
    dU2 = MatrixXd::Zero(ny, nx);
    dU3 = MatrixXd::Zero(ny, nx);

    for (int i = 3; i < ny - 3; ++i)
    {
        for (int j = 3; j < nx - 3; ++j)
        {
            if (mesh.bctype(i, j) != 0) continue;

            std::vector<double> Up = {
                mesh.U0(i,j), mesh.U1(i,j), mesh.U2(i,j), mesh.U3(i,j)
            };

            auto getNeighbor = [&](int ii, int jj) -> std::vector<double> {
                return {mesh.U0(ii,jj), mesh.U1(ii,jj),
                        mesh.U2(ii,jj), mesh.U3(ii,jj)};
            };

            auto Ur1 = getNeighbor(i, j+1);
            auto Ur2 = getNeighbor(i, j+2);
            auto Ur3 = getNeighbor(i, j+3);
            auto Ul1 = getNeighbor(i, j-1);
            auto Ul2 = getNeighbor(i, j-2);
            auto Ul3 = getNeighbor(i, j-3);
            auto Uu1 = getNeighbor(i-1, j);
            auto Uu2 = getNeighbor(i-2, j);
            auto Uu3 = getNeighbor(i-3, j);
            auto Ud1 = getNeighbor(i+1, j);
            auto Ud2 = getNeighbor(i+2, j);
            auto Ud3 = getNeighbor(i+3, j);

            std::vector<double> U_new = updateCenterCell(
                Up,
                Ur1, Ur2, Ur3, Ul1, Ul2, Ul3,
                Uu1, Uu2, Uu3, Ud1, Ud2, Ud3,
                gamma, dt, dx, dy);

            dU0(i,j) = U_new[0] - Up[0];
            dU1(i,j) = U_new[1] - Up[1];
            dU2(i,j) = U_new[2] - Up[2];
            dU3(i,j) = U_new[3] - Up[3];
        }
    }
}
// ============================================================================
// 3阶 SSP-RK3 时间推进（适配 WENO5 空间离散）
//   Stage 1: U(1)     = U^n + dt·L(U^n)
//   Stage 2: U(2)     = 3/4·U^n + 1/4·(U(1) + dt·L(U(1)))
//   Stage 3: U^{n+1}  = 1/3·U^n + 2/3·(U(2) + dt·L(U(2)))
// ============================================================================
void updateMesh(Mesh& mesh, double dt, int rank, int num_procs)
{
    // 保存 U^n
    const MatrixXd U0n = mesh.U0, U1n = mesh.U1,
                   U2n = mesh.U2, U3n = mesh.U3;

    MatrixXd dU0, dU1, dU2, dU3;

    // -------- Stage 1: U(1) = U^n + dt·L(U^n) --------
    computeRHS(mesh, dt, dU0, dU1, dU2, dU3);

    mesh.U0 = U0n + dU0;
    mesh.U1 = U1n + dU1;
    mesh.U2 = U2n + dU2;
    mesh.U3 = U3n + dU3;
    exchangeConservativeColumns(mesh, rank, num_procs);

    // -------- Stage 2: U(2) = 3/4·U^n + 1/4·(U(1) + dt·L(U(1))) --------
    computeRHS(mesh, dt, dU0, dU1, dU2, dU3);

    mesh.U0 = 0.75*U0n + 0.25*(mesh.U0 + dU0);
    mesh.U1 = 0.75*U1n + 0.25*(mesh.U1 + dU1);
    mesh.U2 = 0.75*U2n + 0.25*(mesh.U2 + dU2);
    mesh.U3 = 0.75*U3n + 0.25*(mesh.U3 + dU3);
    exchangeConservativeColumns(mesh, rank, num_procs);

    // -------- Stage 3: U^{n+1} = 1/3·U^n + 2/3·(U(2) + dt·L(U(2))) --------
    computeRHS(mesh, dt, dU0, dU1, dU2, dU3);

    mesh.U0 = (1.0/3.0)*U0n + (2.0/3.0)*(mesh.U0 + dU0);
    mesh.U1 = (1.0/3.0)*U1n + (2.0/3.0)*(mesh.U1 + dU1);
    mesh.U2 = (1.0/3.0)*U2n + (2.0/3.0)*(mesh.U2 + dU2);
    mesh.U3 = (1.0/3.0)*U3n + (2.0/3.0)*(mesh.U3 + dU3);
    exchangeConservativeColumns(mesh, rank, num_procs);
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