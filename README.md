# ZuikakuCFD — 2D Euler 方程并行求解器

本项目提供两种高阶求解器实现，均基于 MPI + Eigen，支持在多进程（列方向域分解）上并行求解：

| 格式 | 描述 | 精度 | 特点 |
|------|------|------|------|
| **WENO5 FV** | 有限体积 + WENO5 重构 | 5阶 | 不连续解、激波锐利、经典稳定 |
| **DG** | 间断伽辽金 + Legendre 基 | 可调(P=1,2,3) | 无需外部模板、通量局部化、多项式表示 |

两个求解器都使用 **HLLC Riemann 通量** + **SSP-RK3 时间推进**。

---

##  求解器选择指南

### 🔵 使用 WENO5 FV 版本（`src_WENO/`）
**推荐情景**：
- 经典 FV 方法，生产级代码
- 不连续解（激波、接触间断）问题
- 需要高阶精度且稳定性保证
- 熟悉有限体积方法的用户

**编译与运行**：
```bash
make weno
mpirun -np 4 ./solver_WENO ./mesh_data 1e-4 1000
```

### 🔴 使用 DG 版本（`src_DG/`）
**推荐情景**：
- 研究或对比多项式方法
- 光滑或弱间断解
- 无需宽模板依赖（仅 1 层 ghost）
- 研究高阶方法的优缺点

**编译与运行**：
```bash
make dg
# DG: P=1 (4模式/格) 或 P=2 (9模式/格)，见 src_DG/fluid.h 中 DG_P
mpirun -np 4 ./solver_DG ./mesh_data auto 1000
```

---

##  核心特性

### 共有特性
- **MPI 并行**：按列（x 方向）分域，各进程维护 ghost 列并通信
- **稳定时间推进**：SSP-RK3（Strong Stability Preserving）三阶格式
- **通量求解**：HLLC Riemann 求解器，处理 2D 守恒律
- **网格输入**：从磁盘读取参数与初始条件（`params.txt` + `.dat` 文件）
- **周期输出**：每 100 步输出结果到 `result/<step>/` 目录

### WENO5 FV 特有
- **五阶精度**：WENO5（Weighted Essentially Non-Oscillatory）重构
- **灵活重构**：支持 MUSCL（二阶）与 WENO5（五阶）切换
- **3 层 ghost**：WENO5 需要 6 点模板

### DG 特有
- **可调阶数**：P=1（二阶）/ P=2（三阶）/ P=3（四阶）
- **多项式表示**：Legendre 正交基，无需外部重构模板
- **1 层 ghost**：DG 界面处直接多项式求值
- **Cockburn-Shu 限制器**：消除 Gibbs 振荡

---

##  网格数据格式

### 参数文件 (`mesh/params.txt`)
```
nx ny da gamma
```

- `nx`, `ny`：全局格子数（不含 ghost）
- `da`：均匀网格间距（$\Delta x = \Delta y = da$）
- `gamma`：比热比（理想气体，通常 1.4）

### 数据文件（纯文本矩阵，行优先）
每个 `.dat` 文件的行数为 `ny`，列数为 `nx`：

| 文件 | 含义 |
|------|------|
| `rho.dat` | 密度 $\rho$ |
| `u.dat` | x 方向速度 $u$ |
| `v.dat` | y 方向速度 $v$ |
| `p.dat` | 压力 $p$ |
| `bctype.dat` | 边界类型标记（整数） |

**边界类型编码**：
- `0`：内部流体单元
- `-1`：固壁 ghost（复制边界条件）
- `-3`：MPI 进程间 ghost（通信填充）

---

##  依赖 (Requirements)

- C++17 兼容编译器
- MPI：`mpic++` / `mpirun`
- Eigen3（头文件即可）

>  请确保系统已安装 MPI（OpenMPI、MPICH 等）并可从命令行使用 `mpic++`。

---

##  构建（Build）

### 构建两个版本（推荐）
```bash
make all
```

### 仅构建 WENO5 版本
```bash
make weno
```
输出：`solver_WENO`

### 仅构建 DG 版本
```bash
make dg
```
输出：`solver_DG`（默认 P=1）

### 配置 DG 阶数
编辑 [src_DG/fluid.h](src_DG/fluid.h#L19)，修改 `DG_P` 值：
```cpp
static constexpr int DG_P = 1;  // 改为 1, 2 或 3
```
然后重新编译：`make clean && make dg`

| DG_P | 阶数 | 模式数 | CFL 限制 |
|------|------|--------|---------|
| 1 | 2 阶 | 4 | ~0.167 |
| 2 | 3 阶 | 9 | ~0.100 |
| 3 | 4 阶 | 16 | ~0.071 |

---

##  运行（Run）

### WENO5 FV 版本

```bash
mpirun -np <进程数> ./solver_WENO <网格目录> <dt> <时间步数>
```

示例：
```bash
mpirun -np 4 ./solver_WENO ./2D-Riemann 1e-4 1000
```

**参数说明**：
- `<网格目录>`：包含 `params.txt` 和 `*.dat` 的目录
- `<dt>`：固定时间步长
- `<时间步数>`：总迭代步数

### DG 版本

```bash
mpirun -np <进程数> ./solver_DG <网格目录> <dt_或_auto> <时间步数>
```

示例（自动 CFL）：
```bash
mpirun -np 4 ./solver_DG ./2D-Riemann auto 1000
```

示例（固定时间步）：
```bash
mpirun -np 4 ./solver_DG ./2D-Riemann 1e-5 1000
```

**参数说明**：
- `<dt_或_auto>`：`auto`（自动 CFL 计算）或指定数值
- DG auto-CFL 基于 Cockburn-Shu 稳定准则：$\Delta t \leq \frac{C}{2P+1} \cdot \frac{h}{\lambda_{\max}}$

---

## 关键模块说明

### 共有模块

#### 网格分割与通信
- **`splitMeshVertically()`**：沿 x 方向均匀分割网格给 n 个进程，包含 ghost 层
- **`exchangeColumns()` 与 `exchangeConservativeColumns()`**：各进程交换 ghost 列数据（支持守恒变量并发通信）

#### 数值通量
- **HLLC 求解器**：两版本通用，处理守恒律 Riemann 问题

### WENO5 FV 版本专用

#### 重构格式
| 格式 | 阶数 | 模板大小 | 函数 | 特点 |
|------|------|---------|------|------|
| MUSCL | 2 阶 | 4 点 | `muscl_reconstruct()` | 低耗、通用限制器 |
| **WENO5** | **5 阶** | **6 点** | `weno5_reconstruct()` | **高精度、自适应光滑性** |

WENO5 采用 Jiang-Shu 光滑指示子，自动在光滑区使用全 5 阶精度，在间断处退化为低阶保证稳定性。

#### 空间离散
- **有限体积法**：格心点处迹近，4 个界面通量求和
- **界面通量**：WENO5 重构左右状态 → HLLC 通量 → 体积平均

#### 时间推进
- **SSP-RK3**：三阶强稳定性保持格式
  ```
  Stage 1: U(1)     = U^n + dt·L(U^n)
  Stage 2: U(2)     = 3/4·U^n + 1/4·(U(1) + dt·L(U(1)))
  Stage 3: U^{n+1}  = 1/3·U^n + 2/3·(U(2) + dt·L(U(2)))
  ```

### DG 版本专用

#### Legendre 正交基
- **基函数**：$\varphi_m(\xi,\eta) = L_{p_x}(\xi) \cdot L_{p_y}(\eta)$，其中 $\xi,\eta \in [-1,1]$
- **索引映射**：$m = p_x(P+1) + p_y$
- **正交性**：$\int_{-1}^{1} L_i(x) L_j(x) dx = \frac{2}{2i+1} \delta_{ij}$

#### 弱形式离散
对每个单元的每个测试函数 $\varphi_k$：

$$M_k \frac{d\hat{u}_k}{dt} = \int_K \left( F \frac{\partial \varphi_k}{\partial \xi} + G \frac{\partial \varphi_k}{\partial \eta} \right) d\xi d\eta - \oint_{\partial K} F^* \cdot \hat{n} \varphi_k d\gamma$$

- **体积分**：Gauss-Legendre 求积（3点，精度到 5次）
- **面积分**：界面处多项式直接求值 + HLLC 数值通flux

#### 限制器
- **Cockburn-Shu 限制器**：用 minmod 作用在线性模式（$p_x+p_y\leq 1$）
- **自适应衰减**：若线性模式被限制，高阶模式置零（退化为分段线性）
- **作用**：消除激波/间断附近的 Gibbs 振荡

#### CFL 准则（自动 DT）
- **理论限制**：$\Delta t \leq \frac{C}{2P+1} \cdot \frac{h}{\lambda_{\max}}$（C=0.9 安全系数）
- **函数**：`computeMaxSpeedDG()` 对每格 4 个角点求多项式值 → 波速

---

##  输出格式

每 100 步时程序自动输出结果：

```
result/<step>/
  ├── U0_0.dat     # Rank 0 的守恒变量 ρ
  ├── U1_0.dat     # Rank 0 的动量 ρu
  ├── U2_0.dat     # Rank 0 的动量 ρv
  ├── U3_0.dat     # Rank 0 的总能量 E
  ├── U0_1.dat     # Rank 1 的 ρ
  └── ...          # 其他进程的输出
```

每个文件为纯文本矩阵，可用 Python/NumPy、Matlab 等工具读取并可视化。

---

##  物理模型

### 守恒律系统（2D 欧拉方程）
$$\frac{\partial \mathbf{U}}{\partial t} + \frac{\partial \mathbf{F}}{\partial x} + \frac{\partial \mathbf{G}}{\partial y} = 0$$

其中守恒变量 $\mathbf{U} = [\rho, \rho u, \rho v, E]^T$，物理通量通过原始变量 $(p, u, v)$ 计算。

### CFL 限制

#### WENO5 FV 版本
推荐 CFL ≤ 0.4：
$$\text{CFL} = \frac{(|u| + |v| + a) \cdot dt}{da}$$

#### DG 版本
自动计算（基于 Cockburn-Shu 准则）或用户指定时间步。若指定 `auto`，程序计算安全值：
$$\Delta t_{\text{safe}} = \frac{0.9}{2(2P+1)} \cdot \frac{h}{\lambda_{\max}}$$

其中 $a = \sqrt{\gamma p / \rho}$ 为声速，$\lambda_{\max}$ 为最大特征速度。

##  编译与系统要求

- **编译器**：C++17 兼容（GCC 7+、Clang 5+ 等）
- **MPI**：OpenMPI 3.0+ 或 MPICH 3.0+（需 `mpic++` 在 PATH）
- **Eigen3**：仅需头文件（已在 `#include <eigen3/Eigen/...>` 中）
- **OS**：Linux / macOS / Windows (WSL/Cygwin)

### 验证环境
```bash
which mpic++         # 验证 mpic++ 可用
mpirun --version     # 显示 MPI 版本
```

---

##  后处理与可视化

### Python 后处理示例
```python
import numpy as np
import matplotlib.pyplot as plt

# 读取输出数据
U0 = np.loadtxt('result/100/U0_0.dat')  # 密度
U1 = np.loadtxt('result/100/U1_0.dat')  # x 动量
U2 = np.loadtxt('result/100/U2_0.dat')  # y 动量

rho = U0
u = U1 / rho
v = U2 / rho

# 绘制等高线
plt.figure(figsize=(10, 8))
plt.contourf(rho, levels=20, cmap='RdYlBu_r')
plt.colorbar(label='Density')
plt.title('Density at t=100')
plt.show()
```

可使用 `plot.ipynb` 或 ParaView、Tecplot 等通用可视化工具。

---

##  重要注意

### 通用要求
1. **网格划分**：x 方向列分解，要求 `nx ≥ num_procs`，建议 `nx % num_procs == 0`
2. **输出同步**：Rank 0 创建目录，其他进程通过 `MPI_Barrier` 同步后写入

### WENO5 FV 特定
- **Ghost 层**：每个子域含 **3 层**左右 ghost（WENO5 需 6 点模板）
- **CFL 稳定性**：建议 CFL ≤ 0.4

### DG 特定
- **Ghost 层**：每个子域含 **1 层**左右 ghost（多项式直接求值）
- **阶数配置**：修改 `DG_P` 后必须重新编译
- **CFL 自动**：DG 版本支持 `dt` 参数为 `auto`，自动计算安全时间步

---

##  常见问题与调试

### Q：编译时找不到 Eigen3
```bash
sudo apt-get install libeigen3-dev    # Ubuntu/Debian
brew install eigen                     # macOS
```
确保 `#include <eigen3/Eigen/...>` 的路径正确。

### Q：运行时 MPI 初始化失败
```bash
# 确保 OpenMPI 已正确安装
mpicc --version
mpirun -np 2 ./solver_WENO ./mesh 1e-4 100
```

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
| 激波捕捉 | ○ 需限制器 | ✓ 自适应 |
| 光滑区精度 | ✓ 多项式 | ✓ 高阶 |
| 编程复杂度 | ✓ 相对简洁 | ✗ 复杂 |
| 计算成本 | ✗ 高阶模式多 | ✓ 单值 |
| 教学/研究 | ✓ 清晰结构 | ✓ 工业标准 |

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
