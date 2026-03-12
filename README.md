# ZuikakuCFD — 2D Euler 方程并行求解器

本项目是一个基于 MPI + Eigen 的二维欧拉方程有限体积求解器，支持在多进程（列方向域分解）上并行求解，采用 **WENO5 五阶重构** + **HLLC 通量** + **SSP-RK3 三阶时间推进**。

---

##  核心特性

- **MPI 并行**：按列（x 方向）分域，各进程维护 3 层 ghost 列并通信
- **五阶精度**：WENO5（Weighted Essentially Non-Oscillatory）重构，保持激波锐利性
- **稳定时间推进**：SSP-RK3（Strong Stability Preserving）三阶格式
- **通量求解**：HLLC Riemann 求解器，处理 2D 守恒律
- **灵活重构**：同步支持 MUSCL（二阶）与 WENO5（五阶）
- **网格输入**：从磁盘读取参数与初始条件（`mesh/params.txt` + `.dat` 文件）
- **周期输出**：每 100 步输出结果到 `result/<step>/` 目录

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

在项目根目录下运行：

```sh
make
```

构建成功后，会生成可执行文件 `solver`，以及 `build/` 目录下的中间文件。

---

##  运行（Run）

```
mpirun -np <进程数> ./solver <网格目录> <dt> <时间步数>
```

示例：

```sh
mpirun -np 4 ./solver ./mesh 1e-4 1000
```

- `<网格目录>`：包含 `params.txt` 和 `*.dat` 的目录
- `<dt>`：时间步长
- `<时间步数>`：总迭代步数

---

## 关键模块说明

### 网格分割与通信
- **`splitMeshVertically()`**：沿 x 方向均匀分割网格给 n 个进程，每个子域包含 3 层左右 ghost 层
- **`exchangeColumns()`**：各进程发送内侧 3 列，接收外侧 3 列 ghost 数据（支持 4 个守恒变量并发通信）
- **`exchangeConservativeColumns()`**：批量交换 U0、U1、U2、U3 四个矩阵

### 重构格式
| 格式 | 阶数 | 模板大小 | 函数 | 特点 |
|------|------|---------|------|------|
| MUSCL | 2 阶 | 4 点 | `muscl_reconstruct()` | 低耗、通用限制器 |
| **WENO5** | **5 阶** | **6 点** | `weno5_reconstruct()` | **高精度、自适应光滑性** |

WENO5 采用 Jiang-Shu 光滑指示子，自动在光滑区使用全 5 阶精度，在间断处退化为低阶保证稳定性。

### 通量与时间推进
- **HLLC 通量**：左右状态经过 WENO5 重构后，通过 HLLC Riemann 求解器计算数值通量
- **SSP-RK3**：三阶强稳定性保持格式
  ```
  Stage 1: U(1)     = U^n + dt·L(U^n)
  Stage 2: U(2)     = 3/4·U^n + 1/4·(U(1) + dt·L(U(1)))
  Stage 3: U^{n+1}  = 1/3·U^n + 2/3·(U(2) + dt·L(U(2)))
  ```

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
程序启动时自动计算 CFL 数：
$$\text{CFL} = \frac{(|u| + |v| + a) \cdot dt}{da}$$

建议 CFL ≤ 0.4 确保稳定性（其中 $a = \sqrt{\gamma p / \rho}$ 为声速）。

---

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

1. **网格划分**：x 方向列分解，要求 `nx ≥ num_procs`，建议 `nx % num_procs == 0`
2. **Ghost 层**：每个子域含 3 层左右 ghost，内部计算区域为中央部分
3. **输出同步**：Rank 0 创建目录，其他进程通过 `MPI_Barrier` 同步后写入
4. **WENO5 要求**：需要 3 层 ghost 来获取 6 点重构模板

---

##  常见问题与调试

### Q：编译时找不到 Eigen3
```bash
sudo apt-get install libeigen3-dev    # Ubuntu/Debian
brew install eigen                     # macOS
```
确保 `#include <eigen3/Eigen/...>` 的路径正确。

### 运行时 MPI 初始化失败
```bash
# 确保 OpenMPI 已正确安装
mpicc --version
mpirun -np 2 ./solver ./mesh 1e-4 100
```

### 结果数值爆炸或出现 NaN
- 检查 CFL 数是否过大（应 ≤ 0.4）
- 确认初始条件文件有效
- 尝试减小 `dt` 值

---

##  扩展方向

可考虑加入以下功能：

- **2D 域分解**：同时沿 x、y 方向分割
- **自适应时间步**：根据局部 CFL 动态调整
- **OpenMP 线程化**：内层循环并行化加速
- **GPU 加速**：CUDA/HIP 实现通量计算
- **复杂边界条件**：非等熵边界、流入流出处理
- **源项与化学反应**：扩展到反应流模型
