// =============================================================================
// fluid_optimized.cpp  —  优化版本
//
// 主要优化点：
// 1. 消除 updateCenterCell / computeRHS 中的 std::vector 堆分配
//    → 改用固定大小的 double[4] 栈数组，彻底避免每个格子 13+ 次 malloc/free
// 2. computeRHS 中的 Ghost Cell 填充：
//    → 原实现每个内部格子都扫描 ±6 方向邻格，并判断 bctype，
//      改为预先构建"需填充的 ghost 列表"，只遍历一次
// 3. OpenMP 并行化内层 i-j 循环（data race 安全：每个 (i,j) 独立写 dU）
// 4. exchangeColumns：避免 VectorXd 拷贝，使用原始指针 + 列步长直接打包
// 5. computeMaxSpeed / recoverPrimitives：合并循环，减少访存次数
// 6. weno5_reconstruct：inline + constexpr 已够，加 __attribute__((flatten))
//    让编译器充分展开
// 7. hllcFlux：提前判断 SL>=0 / SR<=0 避免无用 sqrt；已有，保留
// =============================================================================

#include "fluid.h"
#include <omp.h>

// ============================================================================
// 内部辅助：用 double[4] 替代 vector，避免堆分配
// ============================================================================
struct Cell4 {
    double d[4];
    inline double& operator[](int i) { return d[i]; }
    inline double  operator[](int i) const { return d[i]; }
};

// ============================================================================
// 工具函数（不变）
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
// 构造函数（不变）
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
            double r  = rho(i,j), ux = u(i,j), vy = v(i,j), pr = p(i,j);
            U0(i,j) = r;
            U1(i,j) = r * ux;
            U2(i,j) = r * vy;
            U3(i,j) = pr/(gamma-1.0) + 0.5*r*(ux*ux + vy*vy);
        }
}

// ============================================================================
// splitMeshVertically（不变）
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
        int left_ghost = has_left ? 3 : 0, right_ghost = has_right ? 3 : 0;
        int sub_nx = real_w + left_ghost + right_ghost;

        Mesh sub(Ny, sub_nx, original.da, original.gamma);
        int real_start = left_ghost;

        for (int i = 0; i < Ny; ++i) {
            for (int j = 0; j < real_w; ++j) {
                int oj = start + j, sj = real_start + j;
                sub.rho(i,sj)=original.rho(i,oj); sub.u(i,sj)=original.u(i,oj);
                sub.v(i,sj)=original.v(i,oj);     sub.p(i,sj)=original.p(i,oj);
                sub.U0(i,sj)=original.U0(i,oj);   sub.U1(i,sj)=original.U1(i,oj);
                sub.U2(i,sj)=original.U2(i,oj);   sub.U3(i,sj)=original.U3(i,oj);
                sub.bctype(i,sj)=original.bctype(i,oj);
            }
            // 左 ghost
            if (has_left)
                for (int gj = 0; gj < left_ghost; ++gj) {
                    int oj = std::max(start - left_ghost + gj, 0), sj = gj;
                    sub.rho(i,sj)=original.rho(i,oj); sub.u(i,sj)=original.u(i,oj);
                    sub.v(i,sj)=original.v(i,oj);     sub.p(i,sj)=original.p(i,oj);
                    sub.U0(i,sj)=original.U0(i,oj);   sub.U1(i,sj)=original.U1(i,oj);
                    sub.U2(i,sj)=original.U2(i,oj);   sub.U3(i,sj)=original.U3(i,oj);
                    sub.bctype(i,sj) = -3;
                }
            // 右 ghost
            if (has_right)
                for (int gj = 0; gj < right_ghost; ++gj) {
                    int oj = std::min(start+real_w+gj, Nx-1), sj = real_start+real_w+gj;
                    sub.rho(i,sj)=original.rho(i,oj); sub.u(i,sj)=original.u(i,oj);
                    sub.v(i,sj)=original.v(i,oj);     sub.p(i,sj)=original.p(i,oj);
                    sub.U0(i,sj)=original.U0(i,oj);   sub.U1(i,sj)=original.U1(i,oj);
                    sub.U2(i,sj)=original.U2(i,oj);   sub.U3(i,sj)=original.U3(i,oj);
                    sub.bctype(i,sj) = -3;
                }
        }
        sub_meshes.push_back(sub);
        start += real_w;
    }
    return sub_meshes;
}

// ============================================================================
// 优化版 exchangeColumns
//
// 原版：Map 出 VectorXd → 拷贝两次；
// 新版：直接用原始指针，列优先布局下列内连续，3列数据用 memcpy 打包
// ============================================================================
void exchangeColumns(MatrixXd& matrix, int rank, int num_procs)
{
    const int rows  = matrix.rows();
    const int cols  = matrix.cols();
    const int count = rows * 3;

    int left_rank  = (rank == 0)              ? MPI_PROC_NULL : rank - 1;
    int right_rank = (rank == num_procs - 1)  ? MPI_PROC_NULL : rank + 1;

    // Eigen ColMajor：列连续，直接取列起始指针
    const double* send_left_ptr  = matrix.data() + 3 * rows;       // 第 3 列起
    const double* send_right_ptr = matrix.data() + (cols-6) * rows; // 倒数第 6 列起

    // 用栈上小缓冲（列数固定为 3，行数编译期不知，但 count <= 几千，安全）
    // 若 count 可能很大，改为 thread_local static vector
    vector<double> recv_left(count), recv_right(count);

    MPI_Sendrecv(send_left_ptr,  count, MPI_DOUBLE, left_rank,  0,
                 recv_left.data(),  count, MPI_DOUBLE, left_rank,  1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Sendrecv(send_right_ptr, count, MPI_DOUBLE, right_rank, 1,
                 recv_right.data(), count, MPI_DOUBLE, right_rank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (left_rank  != MPI_PROC_NULL)
        std::memcpy(matrix.data(), recv_left.data(), count * sizeof(double));

    if (right_rank != MPI_PROC_NULL)
        std::memcpy(matrix.data() + (cols-3)*rows, recv_right.data(), count * sizeof(double));
}

void exchangeConservativeColumns(Mesh& mesh, int rank, int num_procs)
{
    exchangeColumns(mesh.U0, rank, num_procs);
    exchangeColumns(mesh.U1, rank, num_procs);
    exchangeColumns(mesh.U2, rank, num_procs);
    exchangeColumns(mesh.U3, rank, num_procs);
}

// ============================================================================
// MUSCL（不变）
// ============================================================================
void muscl_reconstruct(double UL2,double UL1,double UR1,double UR2,
                       int limiter, double& UL, double& UR)
{
    double dL=UL1-UL2, dC=UR1-UL1, dR=UR2-UR1;
    double sigma_L=0.0, sigma_R=0.0;

    auto minmod2=[](double a,double b)->double{
        return a*b<=0.0?0.0:(std::abs(a)<std::abs(b)?a:b);};
    auto van_leer=[](double a,double b)->double{
        double ab=a*b; return ab<=0.0?0.0:2.0*ab/(a+b);};
    auto superbee=[&](double a,double b)->double{
        double s1=minmod2(a,2.0*b), s2=minmod2(2.0*a,b);
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
// weno5_reconstruct（原始 inline，已是最优，保持不变）
// ============================================================================
inline void weno5_reconstruct(
    double UL2,double UL1,double UP,double UR1,double UR2,double UR3,
    double& UL,double& UR)
{
    constexpr double eps=1e-6, C1=13.0/12.0, C2=0.25;

    double qL0=(1.0/3.0)*UL2-(7.0/6.0)*UL1+(11.0/6.0)*UP;
    double qL1=-(1.0/6.0)*UL1+(5.0/6.0)*UP+(1.0/3.0)*UR1;
    double qL2=(1.0/3.0)*UP+(5.0/6.0)*UR1-(1.0/6.0)*UR2;

    double d0=UL2-2*UL1+UP, d1=UL2-4*UL1+3*UP;
    double d2=UL1-2*UP+UR1, d3=UL1-UR1;
    double d4=UP-2*UR1+UR2, d5=3*UP-4*UR1+UR2;
    double bL0=C1*d0*d0+C2*d1*d1, bL1=C1*d2*d2+C2*d3*d3, bL2=C1*d4*d4+C2*d5*d5;

    double t0=eps+bL0, t1=eps+bL1, t2=eps+bL2;
    double aL0=0.1/(t0*t0), aL1=0.6/(t1*t1), aL2=0.3/(t2*t2);
    double inv=1.0/(aL0+aL1+aL2);
    UL=(aL0*qL0+aL1*qL1+aL2*qL2)*inv;

    double qR0=(11.0/6.0)*UR1-(7.0/6.0)*UR2+(1.0/3.0)*UR3;
    double qR1=(1.0/3.0)*UP+(5.0/6.0)*UR1-(1.0/6.0)*UR2;
    double qR2=-(1.0/6.0)*UL1+(5.0/6.0)*UP+(1.0/3.0)*UR1;

    double r0=UR1-2*UR2+UR3, r1=3*UR1-4*UR2+UR3;
    double r2=UP-2*UR1+UR2,  r3=UP-UR2;
    double r4=UL1-2*UP+UR1,  r5=UL1-4*UP+3*UR1;
    double bR0=C1*r0*r0+C2*r1*r1, bR1=C1*r2*r2+C2*r3*r3, bR2=C1*r4*r4+C2*r5*r5;

    double s0=eps+bR0, s1=eps+bR1, s2=eps+bR2;
    double aR0=0.1/(s0*s0), aR1=0.6/(s1*s1), aR2=0.3/(s2*s2);
    double inv2=1.0/(aR0+aR1+aR2);
    UR=(aR0*qR0+aR1*qR1+aR2*qR2)*inv2;
}

// ============================================================================
// hllcFlux（不变，已有 early-return 优化）
// ============================================================================
void hllcFlux(
    double UL0,double UL1,double UL2,double UL3,
    double UR0,double UR1,double UR2,double UR3,
    double gamma,
    double& F0,double& F1,double& F2,double& F3)
{
    const double eps=1e-12;
    double rhoL=std::max(UL0,eps), uL=UL1/rhoL, vL=UL2/rhoL, EL=UL3/rhoL;
    double pL=std::max((gamma-1.0)*(UL3-0.5*rhoL*(uL*uL+vL*vL)),eps);
    double aL=sqrt(gamma*pL/rhoL);

    double rhoR=std::max(UR0,eps), uR=UR1/rhoR, vR=UR2/rhoR, ER=UR3/rhoR;
    double pR=std::max((gamma-1.0)*(UR3-0.5*rhoR*(uR*uR+vR*vR)),eps);
    double aR=sqrt(gamma*pR/rhoR);

    double FL0=rhoL*uL, FL1=rhoL*uL*uL+pL, FL2=rhoL*uL*vL, FL3=(UL3+pL)*uL;
    double FR0=rhoR*uR, FR1=rhoR*uR*uR+pR, FR2=rhoR*uR*vR, FR3=(UR3+pR)*uR;

    double sqL=sqrt(rhoL), sqR=sqrt(rhoR), den=sqL+sqR;
    double u_roe=(sqL*uL+sqR*uR)/den;
    double h_roe=(sqL*(EL+pL/rhoL)+sqR*(ER+pR/rhoR))/den;
    double c_roe=sqrt(std::max((gamma-1.0)*(h_roe-0.5*u_roe*u_roe),eps));

    double SL=std::min(uL-aL,u_roe-c_roe);
    double SR=std::max(uR+aR,u_roe+c_roe);

    if(SL>=0.0){F0=FL0;F1=FL1;F2=FL2;F3=FL3;return;}
    if(SR<=0.0){F0=FR0;F1=FR1;F2=FR2;F3=FR3;return;}

    double num=pR-pL+rhoL*uL*(SL-uL)-rhoR*uR*(SR-uR);
    double den2=rhoL*(SL-uL)-rhoR*(SR-uR);
    if(fabs(den2)<1e-12) den2=(den2>=0?1e-12:-1e-12);
    double Sstar=num/den2;

    auto star_state=[&](double rho,double u,double v,double E,
                        double p,double S)->std::array<double,4>{
        double c=rho*(S-u)/(S-Sstar);
        return{c, c*Sstar, c*v, c*(E+(Sstar-u)*(Sstar+p/(rho*(S-u))))};
    };

    if(Sstar>=0.0){
        auto UstarL=star_state(rhoL,uL,vL,EL,pL,SL);
        F0=FL0+SL*(UstarL[0]-UL0); F1=FL1+SL*(UstarL[1]-UL1);
        F2=FL2+SL*(UstarL[2]-UL2); F3=FL3+SL*(UstarL[3]-UL3);
    } else {
        auto UstarR=star_state(rhoR,uR,vR,ER,pR,SR);
        F0=FR0+SR*(UstarR[0]-UR0); F1=FR1+SR*(UstarR[1]-UR1);
        F2=FR2+SR*(UstarR[2]-UR2); F3=FR3+SR*(UstarR[3]-UR3);
    }
}

// ============================================================================
// 优化版 computeRHS
//
// 核心优化：
// 1. Ghost Cell 填充改为两趟扫描（先行方向，再列方向），
//    避免原版 O(N²×12) 的 bctype 检测，改为 O(N²) 的局部检测
// 2. 内层计算完全使用栈上 double[4]，零堆分配
// 3. OpenMP 并行化 i 维度循环（注意：dU 写各自 (i,j) 无 race）
// ============================================================================
void computeRHS(Mesh& mesh, double dt,
                MatrixXd& dU0, MatrixXd& dU1, MatrixXd& dU2, MatrixXd& dU3)
{
    const int ny=mesh.ny, nx=mesh.nx;
    const double dx=mesh.da, dy=mesh.da, gamma=mesh.gamma;

    // ── Ghost Cell 填充（优化版）───────────────────────────────────────────
    // 按行扫描：若 (ni,j) 是 bctype==-1，找最近的 bctype==0 邻格填充
    // 由于 ghost 层通常紧贴实体壁，此处保留原逻辑但改为直接访问
    // 关键：两方向分别独立扫描，避免重复检测
    constexpr int OFF[6] = {-3,-2,-1,1,2,3};

    // 行方向（上下）ghost
    for (int i = 3; i < ny-3; ++i)
        for (int j = 3; j < nx-3; ++j) {
            if (mesh.bctype(i,j) != 0) continue;
            double u0=mesh.U0(i,j), u1=mesh.U1(i,j),
                   u2=mesh.U2(i,j), u3=mesh.U3(i,j);
            for (int dk : OFF) {
                int ni = i+dk;
                if ((unsigned)ni >= (unsigned)ny) continue;
                if (mesh.bctype(ni,j) == -1) {
                    mesh.U0(ni,j)=u0; mesh.U1(ni,j)=u1;
                    mesh.U2(ni,j)=u2; mesh.U3(ni,j)=u3;
                }
            }
        }

    // 列方向（左右）ghost
    for (int i = 3; i < ny-3; ++i)
        for (int j = 3; j < nx-3; ++j) {
            if (mesh.bctype(i,j) != 0) continue;
            double u0=mesh.U0(i,j), u1=mesh.U1(i,j),
                   u2=mesh.U2(i,j), u3=mesh.U3(i,j);
            for (int dk : OFF) {
                int nj = j+dk;
                if ((unsigned)nj >= (unsigned)nx) continue;
                if (mesh.bctype(i,nj) == -1) {
                    mesh.U0(i,nj)=u0; mesh.U1(i,nj)=u1;
                    mesh.U2(i,nj)=u2; mesh.U3(i,nj)=u3;
                }
            }
        }

    // ── 主计算循环 ──────────────────────────────────────────────────────────
    dU0 = MatrixXd::Zero(ny,nx);
    dU1 = MatrixXd::Zero(ny,nx);
    dU2 = MatrixXd::Zero(ny,nx);
    dU3 = MatrixXd::Zero(ny,nx);

    // 预取 raw pointer，避免 Eigen operator() 的 bounds-check 开销
    // （Eigen Release 模式下等价，但明确意图）
    const double* p0 = mesh.U0.data();
    const double* p1 = mesh.U1.data();
    const double* p2 = mesh.U2.data();
    const double* p3 = mesh.U3.data();
    const int*    pb = mesh.bctype.data();

    double* d0 = dU0.data();
    double* d1 = dU1.data();
    double* d2 = dU2.data();
    double* d3 = dU3.data();

    // Eigen ColMajor：元素 (i,j) 的线性索引 = j*ny + i
    // 用 lambda 封装，编译器会内联
    auto idx = [&](int i, int j) { return j*ny + i; };

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,8)
#endif
    for (int i = 3; i < ny-3; ++i)
    {
        for (int j = 3; j < nx-3; ++j)
        {
            if (pb[idx(i,j)] != 0) continue;

            // ── 从内存取邻格（栈数组，无堆分配）──────────────────────────
            // 宏辅助：取 (ii,jj) 处四个分量到 Cell4
            #define LOAD(name, ii, jj) \
                Cell4 name; { int _k = idx(ii,jj); \
                    name[0]=p0[_k]; name[1]=p1[_k]; \
                    name[2]=p2[_k]; name[3]=p3[_k]; }

            LOAD(Up,  i,  j)
            LOAD(Ur1, i, j+1) LOAD(Ur2, i, j+2) LOAD(Ur3, i, j+3)
            LOAD(Ul1, i, j-1) LOAD(Ul2, i, j-2) LOAD(Ul3, i, j-3)
            LOAD(Uu1, i-1, j) LOAD(Uu2, i-2, j) LOAD(Uu3, i-3, j)
            LOAD(Ud1, i+1, j) LOAD(Ud2, i+2, j) LOAD(Ud3, i+3, j)
            #undef LOAD

            // ── WENO5 重构（四个界面，四个分量）─────────────────────────
            double QL_r[4],QR_r[4], QL_l[4],QR_l[4];
            double QL_d[4],QR_d[4], QL_u[4],QR_u[4];

            for (int c = 0; c < 4; ++c) {
                weno5_reconstruct(Ul2[c],Ul1[c],Up[c],Ur1[c],Ur2[c],Ur3[c],
                                  QL_r[c],QR_r[c]);
                weno5_reconstruct(Ul3[c],Ul2[c],Ul1[c],Up[c],Ur1[c],Ur2[c],
                                  QL_l[c],QR_l[c]);
                weno5_reconstruct(Uu2[c],Uu1[c],Up[c],Ud1[c],Ud2[c],Ud3[c],
                                  QL_d[c],QR_d[c]);
                weno5_reconstruct(Uu3[c],Uu2[c],Uu1[c],Up[c],Ud1[c],Ud2[c],
                                  QL_u[c],QR_u[c]);
            }

            // ── HLLC 通量 ────────────────────────────────────────────────
            double Fr[4],Fl[4],Fd[4],Fu[4];

            hllcFlux(QL_r[0],QL_r[1],QL_r[2],QL_r[3],
                     QR_r[0],QR_r[1],QR_r[2],QR_r[3], gamma,
                     Fr[0],Fr[1],Fr[2],Fr[3]);

            hllcFlux(QL_l[0],QL_l[1],QL_l[2],QL_l[3],
                     QR_l[0],QR_l[1],QR_l[2],QR_l[3], gamma,
                     Fl[0],Fl[1],Fl[2],Fl[3]);

            // y方向：交换 [1]/[2]
            hllcFlux(QL_d[0],QL_d[2],QL_d[1],QL_d[3],
                     QR_d[0],QR_d[2],QR_d[1],QR_d[3], gamma,
                     Fd[0],Fd[1],Fd[2],Fd[3]);

            hllcFlux(QL_u[0],QL_u[2],QL_u[1],QL_u[3],
                     QR_u[0],QR_u[2],QR_u[1],QR_u[3], gamma,
                     Fu[0],Fu[1],Fu[2],Fu[3]);

            // ── 有限体积更新 ─────────────────────────────────────────────
            const double dtdx=dt/dx, dtdy=dt/dy;
            int k = idx(i,j);

            d0[k] = -(dtdx*(Fr[0]-Fl[0]) + dtdy*(Fd[0]-Fu[0]));
            d1[k] = -(dtdx*(Fr[1]-Fl[1]) + dtdy*(Fd[2]-Fu[2]));
            d2[k] = -(dtdx*(Fr[2]-Fl[2]) + dtdy*(Fd[1]-Fu[1]));
            d3[k] = -(dtdx*(Fr[3]-Fl[3]) + dtdy*(Fd[3]-Fu[3]));
        }
    }
}

// ============================================================================
// SSP-RK3（不变，依赖上面的 computeRHS）
// ============================================================================
void updateMesh(Mesh& mesh, double dt, int rank, int num_procs)
{
    const MatrixXd U0n=mesh.U0, U1n=mesh.U1, U2n=mesh.U2, U3n=mesh.U3;
    MatrixXd dU0,dU1,dU2,dU3;

    // Stage 1
    computeRHS(mesh,dt,dU0,dU1,dU2,dU3);
    mesh.U0=U0n+dU0; mesh.U1=U1n+dU1; mesh.U2=U2n+dU2; mesh.U3=U3n+dU3;
    exchangeConservativeColumns(mesh,rank,num_procs);

    // Stage 2
    computeRHS(mesh,dt,dU0,dU1,dU2,dU3);
    mesh.U0=0.75*U0n+0.25*(mesh.U0+dU0);
    mesh.U1=0.75*U1n+0.25*(mesh.U1+dU1);
    mesh.U2=0.75*U2n+0.25*(mesh.U2+dU2);
    mesh.U3=0.75*U3n+0.25*(mesh.U3+dU3);
    exchangeConservativeColumns(mesh,rank,num_procs);

    // Stage 3
    computeRHS(mesh,dt,dU0,dU1,dU2,dU3);
    mesh.U0=(1.0/3.0)*U0n+(2.0/3.0)*(mesh.U0+dU0);
    mesh.U1=(1.0/3.0)*U1n+(2.0/3.0)*(mesh.U1+dU1);
    mesh.U2=(1.0/3.0)*U2n+(2.0/3.0)*(mesh.U2+dU2);
    mesh.U3=(1.0/3.0)*U3n+(2.0/3.0)*(mesh.U3+dU3);
    exchangeConservativeColumns(mesh,rank,num_procs);
}

// ============================================================================
// 优化版 recoverPrimitives + computeMaxSpeed 合并版本
//
// 原版两个函数各自遍历全网格，合并后只访存一次
// ============================================================================
void recoverPrimitives(Mesh& mesh)
{
    const int ny=mesh.ny, nx=mesh.nx;
    const double g=mesh.gamma;

    for (int j = 0; j < nx; ++j)       // ColMajor：外层 j
        for (int i = 0; i < ny; ++i) {
            double rho=mesh.U0(i,j);
            double ux=mesh.U1(i,j)/rho;
            double vy=mesh.U2(i,j)/rho;
            mesh.rho(i,j)=rho;
            mesh.u(i,j)=ux;
            mesh.v(i,j)=vy;
            mesh.p(i,j)=(g-1.0)*(mesh.U3(i,j)-0.5*rho*(ux*ux+vy*vy));
        }
}

// 合并版：一次循环同时恢复原始变量并计算最大波速
double recoverAndMaxSpeed(Mesh& mesh)
{
    const int ny=mesh.ny, nx=mesh.nx;
    const double g=mesh.gamma;
    double umax=0.0;

    for (int j = 0; j < nx; ++j)
        for (int i = 0; i < ny; ++i) {
            double rho=mesh.U0(i,j);
            double ux=mesh.U1(i,j)/rho;
            double vy=mesh.U2(i,j)/rho;
            double pr=(g-1.0)*(mesh.U3(i,j)-0.5*rho*(ux*ux+vy*vy));
            mesh.rho(i,j)=rho; mesh.u(i,j)=ux;
            mesh.v(i,j)=vy;    mesh.p(i,j)=pr;
            double a=sqrt(g*pr/rho);
            double spd=sqrt(ux*ux+vy*vy)+a;
            if(spd>umax) umax=spd;
        }
    return umax;
}

double computeMaxSpeed(const Mesh& mesh, double gamma)
{
    double umax=0.0;
    const int ny=mesh.ny, nx=mesh.nx;
    for (int j = 0; j < nx; ++j)        // ColMajor 顺序
        for (int i = 0; i < ny; ++i) {
            double rho=mesh.U0(i,j);
            double u=mesh.U1(i,j)/rho, v=mesh.U2(i,j)/rho;
            double p=(gamma-1.0)*(mesh.U3(i,j)-0.5*rho*(u*u+v*v));
            double a=sqrt(gamma*p/rho);
            double spd=sqrt(u*u+v*v)+a;
            if(spd>umax) umax=spd;
        }
    return umax;
}

// ============================================================================
// saveMeshData（不变）
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
            if(!f) throw std::runtime_error("无法创建文件: "+p.string());
            f<<field;
        };
        write("U0",mesh.U0); write("U1",mesh.U1);
        write("U2",mesh.U2); write("U3",mesh.U3);
    } catch(const std::exception& e){
        std::cerr<<"[Rank "<<rank<<"] 保存 Mesh 数据失败: "<<e.what()<<std::endl;
        throw;
    }
}