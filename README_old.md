# ZuikakuCFD — 二维欧拉方程高阶并行求解器

![License](https://img.shields.io/badge/license-MIT-blue)
![Language](https://img.shields.io/badge/language-C%2B%2B17-brightgreen)
![MPI](https://img.shields.io/badge/parallel-MPI-orange)

## 项目概览

ZuikakuCFD 是一个专为**二维可压缩 Euler 方程**设计的高性能并行求解框架，提供两种互补的数值方法实现：

| 求解方法 | 空间离散 | 时间推进 | 精度 | 特色 | 位置 |
|---------|--------|--------|------|------|------|
| **WENO5 FV** | WENO5 重构 | SSP-RK3 | 5 阶 | 经典稳定，激波锐利 | [src_WENO/](src_WENO/) |
| **DG** | Legendre 多项式 | SSP-RK3 | 2/3/4 阶 | 无模板依赖，通量局部 | [src_DG/](src_DG/) |

两个求解器均采用 **HLLC 近似 Riemann 求解器** 处理数值通量，使用 **列方向 MPI 域共分解** 实现并行。

---

## 核心特性

### 物理模型
- **守恒律系统**：$$\frac{\partial \mathbf{U}}{\partial t} + \nabla \cdot (\mathbf{F}, \mathbf{G}) = 0$$
- **守恒变量**：$$\mathbf{U} = [\rho, \rho u, \rho v, E]^T$$
- **二维笛卡尔网格**：均匀间距（$$\Delta x = \Delta y = h$$）
- **周期或固壁边界**：支持边界类型标记

### 计算特性
| 特性 | 描述 |
|------|------|
| **MPI 并行** | x 方向列分割，ghost 列通信，批量打包优化 |
| **时间推进** | SSP-RK3（三阶强稳定性保持）|
| **数值通量** | HLLC Riemann 求解器（处理接触间断与激波）|
| **自适应 CFL** | 自动计算安全时间步（基于 lambda-max）|
| **I/O** | 每 100 步输出 `.dat` 格式结果 |

### 精度对比

| 方法 | 最高阶数 | ghost 层数 | 模板宽度 | 适用场景 |
|------|---------|----------|---------|---------|
| MUSCL | 2 阶 | 3 | 4 点 | 低成本、光滑解 |
| **WENO5** | **5 阶** | **3** | **6 点** | **激波、接触间断** |
| DG(P=1) | 2 阶 | 1 | 双线性 | 光滑解或弱间断 |
| DG(P=2) | 3 阶 | 1 | 双二次 | 均衡精度与成本 |
| DG(P=3) | 4 阶 | 1 | 双三次 | 高精度、需小时间步 |

---

## 快速开始

### 编译

```bash
# 编译全部（WENO5 + DG P=1,2,3）
make all

# 仅编译 WENO5
make weno

# 仅编译 DG（默认 P=2）
make dg
```

### 运行

#### WENO5 有限体积格式
```bash
mpirun -np 4 ./solver_WENO ./2D-Riemann 1e-4 1000
# 参数：网格目录 dt 时间步数
```

#### DG 间断伽辽金格式（P=2）
```bash
mpirun -np 4 ./solver_DG2 ./2D-Riemann auto 1000
# 参数：网格目录 {auto|dt值} 时间步数
```

### 可视化

```bash
python plot.ipynb
# 加载结果/result/*/, 绘制密度、压力、速度、Mach数
# 参考详见下文"典型测试用例"的结果展示
```

---

## 方法细节

### WENO5 FV（有限体积法）

**离散形式（有限体积）**：
$$\frac{d\bar{\mathbf{U}}_{i,j}}{dt} = -\frac{1}{\Delta x}\left(F^*_{i+1/2,j} - F^*_{i-1/2,j}\right) - \frac{1}{\Delta y}\left(G^*_{i,j+1/2} - G^*_{i,j-1/2}\right)$$

**重构策略**：
- **MUSCL（二阶）**：4 点模板 + minmod/van Leer/superbee 斜率限制器
- **WENO5（五阶）**：6 点模板，Jiang-Shu 光滑指示子
  - 光滑区：完全 5 阶精度
  - 间断区：自动退化为低阶保证 TVD 性质

**关键函数**：
- [src_WENO/fluid.cpp](src_WENO/fluid.cpp#L234)：`weno5_reconstruct()`
- [src_WENO/fluid.cpp](src_WENO/fluid.cpp#L274)：`hllcFlux()`
- [src_WENO/solver.cpp](src_WENO/solver.cpp)：主求解循环

**优势**：
- 经典成熟、生产级代码
- 激波捕捉能力强
- 需要 3 层 ghost（通信开销大）

---

### DG 间断伽辽金法

#### 基本理论

**弱形式（每单元 K）**：
$$\int_K \phi_k \frac{\partial \mathbf{U}}{\partial t} d\mathbf{x} = -\int_K \left(F \frac{\partial \phi_k}{\partial x} + G \frac{\partial \phi_k}{\partial y}\right) d\mathbf{x} + \oint_{\partial K} \mathbf{F}^*(\mathbf{U}^-, \mathbf{U}^+) \cdot \mathbf{n} \phi_k d\gamma$$

**张量积 Legendre 基函数**：
- **定义**：$$\varphi_{m}(\xi, \eta) = L_{p_x}(\xi) \cdot L_{p_y}(\eta)$$，其中 $$\xi, \eta \in [-1, 1]$$，$$L_p$$ 为 $p$ 阶 Legendre 多项式
- **模式索引映射**：$$m = p_x(P+1) + p_y$$，其中 $$p_x, p_y \in [0, P]$$
- **总自由度**（每变量/单元）：$$N_M = (P+1)^2$$
- **正交性质**：$$\int_{-1}^{1} L_i(x) L_j(x) dx = \frac{2}{2i+1} \delta_{ij}$$

**高斯积分**：
体积分和面积分均采用 Gauss-Legendre 求积，积分精度满足：
- 体积分（使用 3×3 点）：精确到 5 次多项式
- 面积分（沿界面用 3 点）：与多项式度数匹配

| P | 阶数 | 模式/单元 | CFL 安全限 | 推荐应用 |
|---|------|---------|----------|---------|
| 1 | 2 阶 | 4 | 0.167 | 激波问题 |
| 2 | 3 阶 | 9 | 0.100 | **均衡型**（默认）|
| 3 | 4 阶 | 16 | 0.067 | 高精度需求 |

#### DG 限制器（Cockburn-Shu）

**问题背景**：
DG 方法在间断附近容易产生 Gibbs 振荡，特别是在激波或陡峭梯度处。多项式高阶项的过度振荡会导致数值解不物理（例如密度或压力负值）。Cockburn-Shu 限制器的目的是在保持整体高阶精度的前提下，在间断处自动降阶以消除振荡。

**核心思想**：
- 只对**线性模式**（$$p_x + p_y \leq 1$$）应用 minmod 限制器
- 若线性部分被限制，则**清零所有高阶模式**（$$p_x + p_y \geq 2$$）
- 这样激波附近自动退化为分段线性，光滑区保持高阶

**算法步骤**（伪代码）**：

对每个单元 $K_{i,j}$ 的每个守恒变量 $U_l$：

```
1. 提取单元平均值与线性模式系数：
   ū = coeff[0]     # 常数项 φ₀
   ux = coeff[1]    # x方向一阶系数 φ₁(ξ)
   uy = coeff[2]    # y方向一阶系数 φ₂(η)

2. 计算单元与相邻单元的平均值：
   u_L = ū - coeff[1]     # 左相邻单元估计
   u_R = ū + coeff[1]     # 右相邻单元估计
   u_B = ū - coeff[2]     # 下相邻单元估计
   u_T = ū + coeff[2]     # 上相邻单元估计

3. 应用 minmod 函数限制线性系数：
   ux_lim = minmod(ux, u_L - ū, u_R - ū)
   uy_lim = minmod(uy, u_B - ū, u_T - ū)

4. 检测是否激活限制器：
   if (ux_lim != ux) OR (uy_lim != uy):
      # 线性模式被修改，清零所有高阶模式
      for m in [3, ..., N_M-1]:
         coeff[m] = 0
      coeff[1] = ux_lim
      coeff[2] = uy_lim
```

**Minmod 函数定义**：
$$\text{minmod}(a, b, c) = \begin{cases}
\text{sign}(a) \min(|a|, |b|, |c|) & \text{if } a,b,c \text{ 同号} \\
0 & \text{otherwise}
\end{cases}$$

**物理意义**：
- **光滑区**：线性限制器不激活，高阶模式保留 → 仍是高阶精度
- **激波附近**：线性限制器激活，高阶模式置零 → 退化为 MUSCL (2 阶)，保证单调性
- **自动识别**：无需显式的激波检测器，梯度本身就指示何时需要限制

**典型行为**：
| 区域 | 限制状态 | 多项式度数 | 精度 |
|------|---------|---------|------|
| 光滑内部 | 不激活 | P | P+1 阶 |
| 温和梯度 | 部分激活 | < P | 2-3 阶（过渡）|
| 激波/间断 | 激活 | 1 | 2 阶（保证单调性）|

**关键函数**：
- [src_DG/fluid.cpp](src_DG/fluid.cpp#L79)：`applyCockburnShuLimiter()`
- [src_DG/fluid.cpp](src_DG/fluid.cpp#L120)：`minmod3()` 
- [src_DG/solver.cpp](src_DG/solver.cpp#L180)：每个 RK 阶段后调用限制器

**优势与权衡**：
✓ **优势**：
- 自动、无参数的激波捕捉
- 全面保证物理约束（密度、压力非负）
- 激波附近无振荡，保持单调

✗ **权衡**：
- 激波附近精度降至 2 阶
- 高度弯曲区域可能过度限制
- 计算开销（每个 RK 阶段都需要调用）

#### DG 求解流程

**关键函数**：
- [src_DG/fluid.cpp](src_DG/fluid.cpp#L79)：`invMass()`、`dphi2()`
- [src_DG/solver.cpp](src_DG/solver.cpp#L13)：`dgCFLLimit()`、`computeMaxSpeedDG_omp()`

**优势**：
- 仅单层 ghost，通信量（1/3 of WENO5）
- 多项式表示便于后处理（微分、采样、数据压缩）
- 与限制器结合实现自适应精度

---

## 网格格式

### 参数文件 (`params.txt`)
```
nx ny da gamma
```
| 字段 | 含义 | 示例 |
|------|------|------|
| `nx` | 全局 x 方向单元数 | 512 |
| `ny` | 全局 y 方向单元数 | 512 |
| `da` | 均匀网格间距 | 1.0 |
| `gamma` | 比热比 | 1.4 |

### 数据文件（`.dat` ）
纯文本矩阵格式，行优先，尺寸 `ny × nx`：

| 文件 | 变量 | 备注 |
|------|------|------|
| `rho.dat` | 密度 $$\rho$$ | 标量场 |
| `u.dat` | x 速度 $$u$$ | m/s |
| `v.dat` | y 速度 $$v$$ | m/s |
| `p.dat` | 压力 $$p$$ | Pa |
| `bctype.dat` | 边界标记 | 整数矩阵 |

### 边界类型编码
```
 0 : 内部单元（参与计算）
-1 : 固壁 ghost（镜像复制）
-3 : MPI 进程间 ghost（通信填充）
```

---

## 编译与配置

### 系统要求
- **C++ 编译器**：支持 C++17（GCC 7+、Clang 5+、MSVC 2017+）
- **MPI**：OpenMPI、MPICH 等（带 mpic++ 前端编译器）
- **Eigen3**：仅需头文件（已包含或系统库）

### 编译选项

```bash
# 标准编译
make all

# 调试模式（需修改 Makefile，去掉 -O3）
make clean && make all

# 清理
make clean              # 删除 build/ 和可执行文件
make distclean          # 同上 + 删除 report/
```

### DG 多版本生成

Makefile 自动为 $$P = 1, 2, 3$$ 各生成一个可执行文件：

```bash
./solver_DG2   # P=1，二阶
./solver_DG3   # P=2，三阶
./solver_DG4   # P=3，四阶
```

编译时通过 `-DDG_P_VAL=<P>` 宏传递阶数参数。

---

## 运行示例

### 例 1：WENO5 求解 Riemann 问题

```bash
# 编译
make weno

# 在 4 进程上运行，固定时间步
mpirun -np 4 ./solver_WENO ./2D-Riemann 1e-4 5000
# 输出：result/step_<k>/ 中各进程的 U0_*.dat 等文件
```

### 例 2：DG(P=2) 自适应 CFL

```bash
# 编译（默认 P=2）
make dg

# 自动计算 CFL，运行 1000 步
mpirun -np 4 ./solver_DG3 ./2D-Riemann auto 1000
# CFL 由 dgCFLLimit() 确定，dt 自动调整
```

### 例 3：性能对标

```bash
# WENO5
mpirun -np 8 ./solver_WENO ./2D-Riemann 1e-4 10000

# DG(P=1)：同等精度，但模式数少
mpirun -np 8 ./solver_DG2 ./2D-Riemann auto 10000

# 输出：性能统计报告（wall time, load balance）
```

---

## 输出与可视化

### 结果文件结构
```
result/
├── step_0000/
│   ├── U0_0.dat      # rank 0 的密度
│   ├── U1_0.dat      # rank 0 的 x 动量
│   ├── U2_0.dat      # rank 0 的 y 动量
│   └── U3_0.dat      # rank 0 的能量
├── step_0100/
└── ...
```

### 后处理脚本
[plot.ipynb](plot.ipynb)：
1. 加载并拼接多进程数据
2. 恢复原始变量 $$(\rho, u, v, p)$$
3. 计算 Mach 数、速度幅值
4. 绘制 contour 图

**使用**：
```bash
jupyter notebook plot.ipynb
# 修改 data_dir='result', n_splits=<进程数>-1
```

---

## 典型测试用例与算例

### 高马赫数喷流问题（High Mach Jet）

**物理背景**：
超高速喷流在天体物理极端条件下的传播特性研究。喷流从左边界以极音速（Mach ≈ 2700）射出，与周围低速背景气体相互作用，产生复杂的激波、膨胀波和接触间断。

**计算域与边界条件**：

| 参数 | 值 | 备注 |
|------|-----|------|
| 计算域 | $[0, 1] \times [-0.25, 0.25]$ | 长宽比 2:1 |
| 喷流出口 | x=0, $y \in [-0.05, 0.05]$ | 左边界条件 |
| 背景气体 | $(\rho, u, v, p) = (0.5, 0, 0, 0.4127)$ | 整体初值 |
| 喷流射流 | $(\rho, u, v, p) = (5, 800, 0, 0.4127)$ | 出口条件 |
| 比热比 | $\gamma = 5/3$ | 单原子气体 |
| 网格 | $nx=400, ny=200$ | 推荐 CFL $< 0.01$ |

**求解器运行**：
```bash

# DG(P=2) 求解器
mpirun -np 32 ./solver_DG2 HighMachJet 1e-8 100000
```

**结果展示**：

| WENO5 | DG(P=2) |
|-------|---------|
| ![HMJP.png](img/HMJP.png) | ![HMJV.png](img/HMJV.png) |

**特征分析**：
- **激波结构**：喷流头部形成强激波（压力突跃）
- **膨胀扇**：喷流与背景的接触面两侧出现膨胀波
- **接触间断**：密度与压力在接触面处不连续
- **涡量爆发**：速度剪切导致剧烈涡量生成
- **数值鲁棒性**：极速流动考验激波捕捉与 TVD 稳定性

---

## 代码组织

### 文件结构
```
ZuikakuCFD/
├── src_WENO/
│   ├── fluid.h          # 数据结构、API 声明
│   ├── fluid.cpp        # WENO5 FV 实现
│   └── solver.cpp       # MPI 主循环、时间推进
├── src_DG/
│   ├── fluid.h          # DG 网格类、模态系数
│   ├── fluid.cpp        # Legendre 基、弱形式离散
│   └── solver.cpp       # DG 求解器、自适应 CFL
├── 2D-Riemann/          # 示例网格数据
├── HighMachJet/         # 高马赫数喷流算例
├── img/                 # 测试结果展示图
│   ├── HMJP.png         # 喷流压力分布（WENO5）
│   ├── HMJV.png         # 喷流速度幅值（DG）
│   ├── WENO5.png        # Riemann 问题 WENO5 结果
│   ├── DG2.png          # Riemann 问题 DG(P=2) 结果
│   └── DG3.png          # Riemann 问题 DG(P=3) 结果
├── Makefile             # 构建脚本（支持多版本 DG）
├── plot.ipynb           # 可视化后处理脚本
└── README.md            # 本文件
```

### 核心类与接口

#### Mesh 类（共有）
```cpp
class Mesh {
  int nx, ny;              // 网格尺寸（含 ghost）
  double da, gamma;        // 间距、比热比
  
  MatrixXd rho, u, v, p;  // 原始变量
  MatrixXd U0, U1, U2, U3; // 守恒变量平均值
  
  // DG 专用
  vector<MatrixXd> dof[4]; // dof[变量][模式]
  MatrixXi bctype;        // 边界标记
};
```

#### WENO5 重构
```cpp
void weno5_reconstruct(
  double UL2, UL1, UP, UR1, UR2, UR3,
  double& UL, UR);  // 5 点输出左右界面值
```

#### DG 基函数
```cpp
double phi2(int m, double xi, double eta);  // φ_m(ξ,η)
void dphi2(int m, double xi, double eta, 
           double& dxi, double& deta);      // ∇φ_m
```

---

## 性能优化

### WENO5 版本
- **数据结构**：栈分配 `double[4]` 替代 `vector`，消除堆分配
- **优化的 ghost 填充**：两遍扫描而非 O(N²) bctype 检测
- **OpenMP 并行**：内层 i-j 循环并行化，无 data race
- **内联优化**：关键函数 WENO5、HLLC 采用 `inline + __attribute__`

### DG 版本
- **矩阵预计算**：phi(m, ξ_k, η_l) 预缓存
- **批量通信**：所有 4 变量×DG_NM 模式一次通信打包
- **OpenMP 并行**：主循环并行，合并 recover+maxspeed 减访存

### MPI 通信
- **列分割**：沿 x 方向均分，最小化通信边界
- **ghost 层优化**：
  - WENO5：3 层 ghost，但支持批量拷贝
  - DG：1 层 ghost，通信量减 67%
- **中阶段通信**：SSP-RK3 各阶段后调用 ghost 交换

---

## 许可与参考

- **License**：[MIT License](LICENSE)
- **作者**：midway酱
- **参考文献**：
  - Cockburn, B., & Shu, C. W. (2001). Runge–Kutta Discontinuous Galerkin Methods for Convection-Dominated Problems. *Journal of Scientific Computing*, 16(3).
  - Jiang, G. S., & Shu, C. W. (1996). Efficient Implementation of Weighted ENO Schemes. *Journal of Computational Physics*, 126(1).
  - Toro, E. F. (2009). *Riemann Solvers and Numerical Methods for Fluid Dynamics*. Springer.

---

## 常见问题

**Q: 如何选择 WENO5 还是 DG？**  
A: WENO5 适合复杂激波问题；DG 适合光滑或多项式适配场景。精度相当时，DG 通信少。

**Q: DG 为什么需要限制器？**  
A: 间断处多项式通常产生 Gibbs 振荡。Cockburn-Shu 限制器在检测到陡峭梯度时自动激活。

**Q: 如何调整时间步？**  
A: WENO5 使用固定 dt；DG 支持 `auto` 自动 CFL。手动调需修改 `dgCFLLimit()` 或 `dgSafeDt()`。

**Q: ghost 列交换的性能瓶颈？**  
A: 3 层 ghost（WENO5）vs 1 层（DG）。若通信占比 >20%，考虑更多进程并减少边界比。

**Q: 结果数值爆炸或出现 NaN?**  
A: WENO5：检查 CFL 数是否过大（应 ≤ 0.4），尝试减小 `dt`；DG：若手动指定 `dt`，保证满足限制条件，或使用 `auto` 模式。

**Q: DG 版本如何改变求解精度?**  
编辑 [src_DG/fluid.h](src_DG/fluid.h#L19)：
```cpp
static constexpr int DG_P = 2;  // 改为 1, 2 或 3
```
然后重新编译：
```bash
make distclean && make dg
```

---

## 扩展方向

可考虑加入以下功能：
- **2D 域分解**：同时沿 x、y 方向分割，进一步扩展扩展性
- **自适应时间步**：根据局部 CFL 动态调整（特别是 DG）
- **OpenMP 线程化**：内层循环并行化加速
- **GPU 加速**：CUDA/HIP 实现通量计算与限制器
- **复杂边界条件**：非等熵边界、流入流出处理、反射壁面
- **源项与化学反应**：扩展到反应流、燃烧模型
- **DG P4+**：更高阶精度
- **适应性网格细化**（AMR）：动态网格调整

---

## 典型测试用例与算例

### 高马赫数喷流问题（High Mach Jet）

**物理背景**：
超高速喷流在天体物理极端条件下的传播特性研究。喷流从左边界以极音速（Mach ≈ 2700）射出，与周围低速背景气体相互作用，产生复杂的激波、膨胀波和接触间断。

**计算域**：
- 求解域：$[0, 1] \times [-0.25, 0.25]$（长宽比 4:1）
- 喷流出口：x=0, $y \in [-0.05, 0.05]$
- 背景条件：$(\rho, u, v, p) = (0.5, 0, 0, 0.4127)$
- 喷流条件：$(\rho, u, v, p) = (5, 800, 0, 0.4127)$
- 比热比：$\gamma = 5/3$（单原子气体）

**运行示例**：
```bash
# WENO5 求解器（固定时间步）
mpirun -np 4 ./solver_WENO ./HighMachJet 1e-5 5000

# DG(P=2) 求解器（自动 CFL）
mpirun -np 4 ./solver_DG3 ./HighMachJet auto 5000
```

**结果展示**：

| 压力分布（WENO5） | 速度幅值（DG） |
|----------|----------|
| ![HMJP.png](img/HMJP.png) | ![HMJV.png](img/HMJV.png) |

**物理特征**：
- **激波结构**：喷流头部形成强激波（极高压力比）
- **膨胀扇**：喷流与背景接触面两侧产生膨胀波
- **接触间断**：密度跳变但压力連续的接触间断
- **涡量生成**：速度剪切层产生集中涡量
- **数值考验**：极速流动对激波捕捉精度极致要求

---

### 二维黎曼问题（2D Riemann Problem）

**计算设置**：
标准二维黎曼问题是多维激波相互作用的经典test case，网格规模 **512×512**，四个象限独立初始化：

| 象限 | 区域 | $\rho$ | $u$ | $v$ | $p$ |
|------|------|--------|--------|--------|---------|
| 右上 | $x \geq 0.5, y \geq 0.5$ | 1.5 | 0 | 0 | 1.5 |
| 左上 | $x < 0.5, y \geq 0.5$ | 0.5323 | 1.206 | 0 | 0.3 |
| 左下 | $x < 0.5, y < 0.5$ | 0.138 | 1.206 | -1.206 | 0.029 |
| 右下 | $x \geq 0.5, y < 0.5$ | 0.5323 | 0 | -1.206 | 0.3 |

**物理特征**：
中心的四条间断线各产生不同的激波、膨胀波或接触间断，是检验多维格式**对称性、激波捕捉、间断分辨**的理想用例。

**求解器对比**：

```bash
# 编译与运行
make all
mpirun -np 4 ./solver_WENO ./2D-Riemann 1e-4 1000
mpirun -np 4 ./solver_DG3 ./2D-Riemann auto 1000
mpirun -np 4 ./solver_DG4 ./2D-Riemann auto 1000

# 可视化
jupyter notebook plot.ipynb
```

**结果展示**：

| WENO5（5阶） | DG(P=2)（3阶） | DG(P=3)（4阶） |
|-----|-----|-----|
| ![WENO5](img/WENO5.png) | ![DG2](img/DG2.png) | ![DG3](img/DG3.png) |

| 方法 | 激波锐利度 | Ghost 开销 | DOF/单元 | 推荐场景 |
|------|---------|---------|---------|---------|
| **WENO5** | 优秀 | 3层（大） | 1 | 强激波问题 |
| **DG(P=2)** | 良好 | 1层（小） | 9 | 均衡（默认） |
| **DG(P=3)** | 优秀 | 1层（小） | 16 | 高精度需求 |

---



## 常见问题

**Q: 如何选择 WENO5 还是 DG？**  
A: WENO5 适合复杂激波问题；DG 适合光滑或多项式适配场景。精度相当时，DG 通信少。

**Q: DG 为什么需要限制器？**  
A: 间断处多项式通常产生 Gibbs 振荡。Cockburn-Shu 限制器在检测到陡峭梯度时自动激活。

**Q: 如何调整时间步？**  
A: WENO5 使用固定 dt；DG 支持 `auto` 自动 CFL。手动调需修改 `dgCFLLimit()` 或 `dgSafeDt()`。

**Q: ghost 列交换的性能瓶颈？**  
A: 3 层 ghost（WENO5）vs 1 层（DG）。若通信占比 >20%，考虑更多进程并减少边界比。

---

## 常见问题

### Q：结果数值爆炸或出现 NaN
- **WENO5**：检查 CFL 数是否过大（应 ≤ 0.4），尝试减小 `dt`
- **DG**：若手动指定 `dt`，保证满足限制条件；或使用 `auto` 模式让程序自动计算

### Q：DG 版本如何改变求解精度？
编辑 [src_DG/fluid.h](src_DG/fluid.h#L19)：
```cpp
static constexpr int DG_P = 2;  // 改为 1, 2 或 3
```
然后重新编译：
```bash
make distclean
make dg
```

### Q：DG vs WENO5 选择
| 问题 | DG | WENO5 |
|------|-----|--------|
| 激波捕捉 | 需限制器 | 自适应 |
| 光滑区精度 | 多项式优秀 | 高阶精度 |
| 编程复杂度 | 相对简洁 | 较为复杂 |
| 计算成本 | 高阶模式多 | 单值存储 |
| 教学/研究 | 清晰结构 | 工业标准 |

---

##  项目结构

```
ZuikakuCFD/
├── src_WENO/              # 有限体积 WENO5 版本
│   ├── fluid.h            # Mesh 类定义、通用函数声明
│   ├── fluid.cpp          # WENO5 重构、HLLC、updateMesh 实现
│   └── solver.cpp         # 主程序入口、MPI 初始化、时间循环
├── src_DG/                # 间断伽辽金版本
│   ├── fluid.h            # DG Mesh 扩展：模态系数存储
│   ├── fluid.cpp          # Legendre 基、弱形式、限制器
│   └── solver.cpp         # DG 主程序（自动CFL、性能报告）
├── Makefile               # 联合构建脚本（solver_WENO + solver_DG）
├── README.md              # 本文件
├── gen.ipynb              # 网格生成笔记本
├── plot.ipynb             # 结果可视化笔记本
└── 2D-Riemann/            # 示例网格数据
    └── params.txt         # nx ny da gamma

result/                    # 输出目录（运行时创建）
├── 100/
│   ├── U0_0.dat
│   ├── U1_0.dat
│   └── ...
└── 200/
    └── ...
```

---

##  扩展方向

可考虑加入以下功能：

- **2D 域分解**：同时沿 x、y 方向分割，进一步扩展扩展性
- **自适应时间步**：根据局部 CFL 动态调整（特别是 DG）
- **OpenMP 线程化**：内层循环并行化加速
- **GPU 加速**：CUDA/HIP 实现通量计算与限制器
- **复杂边界条件**：非等熵边界、流入流出处理、反射壁面
- **源项与化学反应**：扩展到反应流、燃烧模型
- **DG P4+**：更高阶精度
- **适应性网格细化**（AMR）：动态网格调整
