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

**弱形式（每单元 K）**：
$$\int_K \phi_k \frac{\partial \mathbf{U}}{\partial t} d\mathbf{x} = -\int_K \left(F \frac{\partial \phi_k}{\partial x} + G \frac{\partial \phi_k}{\partial y}\right) d\mathbf{x} + \oint_{\partial K} \mathbf{F}^*(\mathbf{U}^-, \mathbf{U}^+) \cdot \mathbf{n} \phi_k d\gamma$$

**基函数与自由度**：
- **张量积 Legendre 基**：$$\phi_{m}(\xi, \eta) = L_{p_x}(\xi) \cdot L_{p_y}(\eta)$$，其中 $$\xi, \eta \in [-1, 1]$$
- **模式编号**：$$m = p_x(P+1) + p_y$$，其中 $$p_x, p_y \in [0, P]$$
- **总自由度**（每变量/单元）：$$N_M = (P+1)^2$$

| P | 阶数 | 模式/单元 | CFL 安全限 | 推荐应用 |
|---|------|---------|----------|---------|
| 1 | 2 阶 | 4 | 0.167 | 激波问题 |
| 2 | 3 阶 | 9 | 0.100 | **均衡型**（默认）|
| 3 | 4 阶 | 16 | 0.067 | 光滑问题 |

**限制器（Cockburn-Shu）**：
在激波/陡峭梯度处自动激活，
1. 对线性模式应用 minmod
2. 若激活，清零高阶模式

**关键函数**：
- [src_DG/fluid.cpp](src_DG/fluid.cpp#L79)：`invMass()`、`dphi2()`
- [src_DG/solver.cpp](src_DG/solver.cpp#L13)：`dgCFLLimit()`、`computeMaxSpeedDG_omp()`

**优势**：
- 仅单层 ghost，通信少
- 多项式表示便于后处理（微分、采样）
- 每单元自由度多（存储、计算成本）

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

## 典型测试用例

### 高马赫数喷流问题（High Mach Jet）

**物理背景**：
超高速喷流在天体物理极端条件下的传播特性研究。喷流从左边界以极音速（Mach ≈ 2700）射出，与周围低速背景气体相互作用，产生复杂的激波、膨胀波和接触间断。

**计算域与边界条件**：

| 参数 | 值 | 备注 |
|------|-----|------|
| 计算域 | $[0, 1] \times [-0.25, 0.25]$ | 长宽比 4:1 |
| 喷流出口 | x=0, $y \in [-0.05, 0.05]$ | 左边界条件 |
| 背景气体 | $(\rho, u, v, p) = (0.5, 0, 0, 0.4127)$ | 整体初值 |
| 喷流射流 | $(\rho, u, v, p) = (5, 800, 0, 0.4127)$ | 出口条件 |
| 比热比 | $\gamma = 5/3$ | 单原子气体 |
| 网格 | $nx=400, ny=200$ | 推荐 CFL $< 0.01$ |

**网格初始化代码**：
```python
import numpy as np

def generate_jet_mesh(nx=400, ny=200, gamma=5/3, output_dir="HighMachJet"):
    """生成高马赫数喷流网格与初值"""
    x_min, x_max, y_min, y_max = 0.0, 1.0, -0.25, 0.25
    dx = (x_max - x_min) / nx
    dy = (y_max - y_min) / ny
    
    rho = np.full((ny, nx), 0.5)
    u = np.zeros((ny, nx))
    v = np.zeros((ny, nx))
    p = np.full((ny, nx), 0.4127)
    
    # 设置喷流出口区域 (x=0, y∈[-0.05, 0.05])
    for i in range(ny):
        yc = y_min + (i + 0.5) * dy
        for j in range(min(3, nx)):  # 左边界前 3 层网格
            if -0.05 <= yc <= 0.05:
                rho[i,j], u[i,j] = 5.0, 800.0
    
    # 保存格式（与求解器兼容）
    np.savetxt(f"{output_dir}/params.txt", [[nx, ny, dx, gamma]])
    np.savetxt(f"{output_dir}/rho.dat", rho)
    np.savetxt(f"{output_dir}/u.dat", u)
    np.savetxt(f"{output_dir}/v.dat", v)
    np.savetxt(f"{output_dir}/p.dat", p)

if __name__ == "__main__":
    generate_jet_mesh()
```

**求解器运行**：
```bash
# WENO5 求解器
mpirun -np 4 ./solver_WENO ./HighMachJet 1e-5 5000

# DG(P=2) 求解器
mpirun -np 4 ./solver_DG3 ./HighMachJet auto 5000
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

## 典型测试用例：二维黎曼问题

### 四象限初始条件

本项目已预装的标准测试为**二维黎曼问题**（2D Riemann Problem），网格规模 **512×512**，四个象限各具不同的初始状态：

| 象限 | 区域 | 密度 $$\rho$$ | 速度 $$u$$ | 速度 $$v$$ | 压力 $$p$$ |
|------|------|--------|--------|--------|---------|
| **右上** | $$x \geq 0.85, y \geq 0.85$$ | 1.5 | 0 | 0 | 1.5 |
| **左上** | $$x < 0.85, y \geq 0.85$$ | 0.5323 | 1.206 | 0 | 0.3 |
| **左下** | $$x < 0.85, y < 0.85$$ | 0.138 | 1.206 | -1.206 | 0.029 |
| **右下** | $$x \geq 0.85, y < 0.85$$ | 0.5323 | 0 | -1.206 | 0.3 |

**特征**：
- 中心四条界面处各形成激波、膨胀波、接触间断等复杂流动结构
- 考验求解器的**多维激波捕捉、间断识别、不对称性处理**能力
- 标准测试用例（见 LeVeque《数值方法》、Lax-Liu 等文献）

### 求解方法对比

下面展示三种求解器在相同初始条件、相同网格、不同精度下的密度云图对比：

#### WENO5 有限体积（5阶精度）
![WENO5 密度分布](img/WENO5.png)
- **激波捕捉**：锐利清晰，无振荡
- **接触间断**：分辨率高
- **计算特点**：3层 ghost，全网格 WENO5 重构

#### DG(P=2) 间断伽辽金（3阶精度）
![DG P=2 密度分布](img/DG2.png)
- **精度适中**：消耗少，结果光滑
- **激波处理**：限制器自动激活，无显著 Gibbs 振荡
- **计算特点**：单层 ghost，9 个模式/单元，自适应 CFL

#### DG(P=3) 间断伽辽金（4阶精度）  
![DG P=3 密度分布](img/DG3.png)
- **精度最高**：多项式度数高，局部光滑区分辨细节
- **细节保留**：激波附近结构清晰
- **计算特点**：单层 ghost，16 个模式/单元，时间步严格

### 复现步骤

```bash
# 1. 编译所有版本
make all

# 2. 运行 WENO5（在 4 进程上，1000 步）
mpirun -np 4 ./solver_WENO ./2D-Riemann 1e-4 1000

# 3. 运行 DG(P=2) 自动 CFL
mpirun -np 4 ./solver_DG3 ./2D-Riemann auto 1000

# 4. 运行 DG(P=3)
mpirun -np 4 ./solver_DG4 ./2D-Riemann auto 1000

# 5. 可视化（修改 plot.ipynb 中的 data_dir）
jupyter notebook plot.ipynb
```

### 结果对比

| 指标 | WENO5 | DG(P=2) | DG(P=3) |
|------|-------|---------|---------|
| 空间精度 | 5 阶 | 3 阶 | 4 阶 |
| 激波锐利度 | 优秀 | 良好 | 优秀 |
| Ghost 通信开销 | 大（3层） | 小（1层） | 小（1层） |
| 单元 DOF | 1 | 9 | 16 |
| 推荐应用 | 强激波问题 | 均衡型 | 高精度需求 |

---

## 实验与验证

### 推荐测试用例
1. **2D Riemann Problem**（已实装）：四象限多灰激波结构，见上文
2. **Sod Shock Tube**：1D Riemann 问题，激波+膨胀波+接触间断
3. **Lax Problem**：高强度激波，考验激波捕捉
4. **Double Mach Reflection**：2D 复杂结构，验证多维格式
5. **Smooth Flow**：平缓流动，验证精度阶

### 精度验证（收敛性）
```bash
# 生成不同网格尺寸 (nx, ny)
# 计算相同物理时间 t_end 的 L2 误差
# 图示化：log(h) vs log(error) 
# 斜率应接近理论阶数
```

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

---

## 联系与反馈

- Issue Tracker：[GitHub Issues](../../issues)
- Pull Request：欢迎贡献优化或新特性

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
- **基函数**：$$\varphi_m(\xi,\eta) = L_{p_x}(\xi) \cdot L_{p_y}(\eta)$$，其中 $$\xi,\eta \in [-1,1]$$
- **索引映射**：$$m = p_x(P+1) + p_y$$
- **正交性**：$$\int_{-1}^{1} L_i(x) L_j(x) dx = \frac{2}{2i+1} \delta_{ij}$$

#### 弱形式离散
对每个单元的每个测试函数 $$\varphi_k$$：

$$M_k \frac{d\hat{u}_k}{dt} = \int_K \left( F \frac{\partial \varphi_k}{\partial \xi} + G \frac{\partial \varphi_k}{\partial \eta} \right) d\xi d\eta - \oint_{\partial K} F^* \cdot \hat{n} \varphi_k d\gamma$$

- **体积分**：Gauss-Legendre 求积（3点，精度到 5次）
- **面积分**：界面处多项式直接求值 + HLLC 数值通flux

#### 限制器
- **Cockburn-Shu 限制器**：用 minmod 作用在线性模式（$$p_x+p_y\leq 1$$）
- **自适应衰减**：若线性模式被限制，高阶模式置零（退化为分段线性）
- **作用**：消除激波/间断附近的 Gibbs 振荡

#### CFL 准则（自动 DT）
- **理论限制**：$$\Delta t \leq \frac{C}{2P+1} \cdot \frac{h}{\lambda_{\max}}$$（C=0.9 安全系数）
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

其中守恒变量 $$\mathbf{U} = [\rho, \rho u, \rho v, E]^T$$，物理通量通过原始变量 $$(p, u, v)$$ 计算。

### CFL 限制

#### WENO5 FV 版本
推荐 CFL ≤ 0.4：
$$\text{CFL} = \frac{(|u| + |v| + a) \cdot dt}{da}$$

#### DG 版本
自动计算（基于 Cockburn-Shu 准则）或用户指定时间步。若指定 `auto`，程序计算安全值：
$$\Delta t_{\text{safe}} = \frac{0.9}{2(2P+1)} \cdot \frac{h}{\lambda_{\max}}$$

其中 $$a = \sqrt{\gamma p / \rho}$$ 为声速，$$\lambda_{\max}$$ 为最大特征速度。

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
