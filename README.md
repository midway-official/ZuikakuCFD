# ZuikakuCFD：二维欧拉方程高阶并行求解器

![License](https://img.shields.io/badge/license-MIT-blue)
![Language](https://img.shields.io/badge/language-C%2B%2B17-brightgreen)
![MPI](https://img.shields.io/badge/parallel-MPI-orange)

## 概览

ZuikakuCFD 是一个为**二维可压缩 Euler 方程**设计的高阶数值求解框架。项目提供两种互补的有限元方法实现：

| 方法 | 离散格式 | 精度 | 特色 | 位置 |
|------|--------|------|------|------|
| **WENO5 有限体积** | WENO5 重构 | 5阶 | 激波锐利、经典成熟 | [src_WENO/](src_WENO/) |
| **DG 间断伽辽金** | Legendre多项式 | 2/3/4阶 | 通信高效、自适应精度 | [src_DG/](src_DG/) |

两种求解器均采用 **HLLC Riemann求解器** 处理数值通量，使用 **列方向 MPI 域分解** 实现并行计算。

---

## 前置知识

本项目面向有限元方法、CFD 数值求解器开发人员。建议了解：
- 二维可压缩 Euler 方程基本形式
- 有限体积法和间断伽辽金法的基本原理
- MPI 并行编程基础
- Riemann 求解器概念

---

## 理论基础

### 守恒律系统

求解的二维可压缩 Euler 方程为：

$$\frac{\partial \mathbf{U}}{\partial t} + \frac{\partial \mathbf{F}(\mathbf{U})}{\partial x} + \frac{\partial \mathbf{G}(\mathbf{U})}{\partial y} = 0$$

其中守恒变量向量和物理通量为：

$$\mathbf{U} = \begin{bmatrix} \rho \\ \rho u \\ \rho v \\ E \end{bmatrix}, \quad \mathbf{F} = \begin{bmatrix} \rho u \\ \rho u^2 + p \\ \rho uv \\ (E+p)u \end{bmatrix}, \quad \mathbf{G} = \begin{bmatrix} \rho v \\ \rho uv \\ \rho v^2 + p \\ (E+p)v \end{bmatrix}$$

其中 $\rho$ 为密度，(u,v) 为速度分量，$E$ 为总能量，压力由理想气体状态方程：

$$p = (\gamma - 1)\left(E - \frac{1}{2}\rho(u^2 + v^2)\right)$$

确定，$\gamma$ 为比热比。

### 时间离散：SSP-RK3

采用三阶强稳定性保持 Runge-Kutta (SSP-RK3) 格式：

$$\mathbf{U}^{(1)} = \mathbf{U}^n + \Delta t \mathcal{L}(\mathbf{U}^n)$$
$$\mathbf{U}^{(2)} = \frac{3}{4}\mathbf{U}^n + \frac{1}{4}\left(\mathbf{U}^{(1)} + \Delta t \mathcal{L}(\mathbf{U}^{(1)})\right)$$
$$\mathbf{U}^{n+1} = \frac{1}{3}\mathbf{U}^n + \frac{2}{3}\left(\mathbf{U}^{(2)} + \Delta t \mathcal{L}(\mathbf{U}^{(2)})\right)$$

其中 $\mathcal{L}(\mathbf{U})$ 为空间离散算子。

---
## 数值方法

## 方法一：WENO5 有限体积格式

### 空间离散

在统一网格 $[x_i, x_{i+1}] \times [y_j, y_{j+1}]$ 上定义单元均值：

```math
\bar{\mathbf{U}}_{i,j}^{n}
=
\frac{1}{\Delta x \Delta y}
\int_{x_i}^{x_{i+1}}
\int_{y_j}^{y_{j+1}}
\mathbf{U}(x,y,t^n)\,dy\,dx
```

有限体积更新格式：

```math
\frac{d\bar{\mathbf{U}}_{i,j}}{dt}
=
-\frac{1}{\Delta x}
\left(
\mathbf{F}_{i+1/2,j}^{*}-\mathbf{F}_{i-1/2,j}^{*}
\right)
-
\frac{1}{\Delta y}
\left(
\mathbf{G}_{i,j+1/2}^{*}-\mathbf{G}_{i,j-1/2}^{*}
\right)
```

数值通量由 Riemann 求解器计算得到：
```math
\text{其中 } \mathbf{F}_{i+1/2,j}^{*} \text{ 和 } \mathbf{G}_{i,j+1/2}^{*} \text{ 为数值通量}
```

### HLLC Riemann 求解器

在界面处求解 Riemann 问题，得到数值通量：

```math
\mathbf{F}^*=
\begin{cases}
\mathbf{F}_L & S_L \ge 0 \\
\mathbf{F}_L^* & S_L < 0 \le S^* \\
\mathbf{F}_R^* & S^* \le 0 < S_R \\
\mathbf{F}_R & S_R \le 0
\end{cases}
```

---

### 方法二：DG 间断伽辽金格式

#### 基本理论

DG 方法在每个单元 $K$ 上用多项式族逼近解。对可压缩 Euler 方程 $\frac{\partial \mathbf{U}}{\partial t} + \nabla \cdot \mathbf{F} = 0$ 在物理单元 $K$（正方形，边长为 $h$）上投影到基函数 $\phi_m$，得到弱形式：

$$\int_K \phi_m \frac{\partial \mathbf{U}}{\partial t} d\mathbf{x} = \int_K \mathbf{F} \cdot \nabla \phi_m d\mathbf{x} + \int_K \mathbf{G} \cdot \nabla \phi_m d\mathbf{x} - \oint_{\partial K} (\mathbf{F}^* \cdot \mathbf{n}) \phi_m \, ds$$

#### 坐标变换与参考元

为提高计算效率，引入参考坐标 $(\xi, \eta) \in [-1,1]^2$。物理坐标与参考坐标的映射为：

$$x = x_c + \frac{h}{2}\xi, \quad y = y_c + \frac{h}{2}\eta$$

其中 $(x_c, y_c)$ 为单元中心。由此导出的微分关系为：

$$\frac{\partial}{\partial x} = \frac{2}{h}\frac{\partial}{\partial \xi}, \quad \frac{\partial}{\partial y} = \frac{2}{h}\frac{\partial}{\partial \eta}$$

积分测度变换：

$$dx\,dy = \frac{h^2}{4} d\xi\,d\eta, \quad ds = \frac{h}{2} d\Gamma_{\text{ref}}$$

其中 $d\Gamma_{\text{ref}}$ 为参考边界上的弧长元。

#### Legendre 正交基

选择二维张量积 Legendre 多项式基：

$$\varphi_m(\xi,\eta) = L_{p_x}(\xi) \cdot L_{p_y}(\eta), \quad \xi,\eta \in [-1,1]$$

其中 $m = p_{x}(P+1) + p_{y}$，$p_{x}, p_{y} \in [0, P]$，$L_{p}$ 为 $p$ 阶 Legendre 多项式，满足正交性：

$$\int_{-1}^1 L_i(\xi) L_j(\xi) d\xi = \frac{2}{2i+1}\delta_{ij}$$

每个单元的解表示为：

$$\mathbf{U}_h(\xi,\eta,t) = \sum_{m=0}^{N_M-1} \hat{\mathbf{U}}_m(t) \varphi_m(\xi,\eta)$$

其中 $N_M = (P+1)^2$ 为该单元的自由度总数。

#### 质量矩阵与逆矩阵推导

解展开为

$$
\boldsymbol{U}_h(\xi, \eta, t) = \sum_{n=0}^{N_M-1} \hat{\boldsymbol{u}}_n(t) \varphi_n(\xi, \eta)
$$

，代入弱形式左端：

$$(\text{LHS})_m = \int_K \varphi_m \left( \sum_n \frac{d\hat{\mathbf{u}}_n}{dt} \varphi_n \right) d\mathbf{x} = \sum_n \frac{d\hat{\mathbf{u}}_n}{dt} \underbrace{\int_K \varphi_m \varphi_n d\mathbf{x}}_{M_{mn}}$$

利用坐标变换，质量矩阵为：

$$M_{mn} = \frac{h^2}{4} \int_{-1}^1 \int_{-1}^1 [L_{p_{x,m}}(\xi)L_{p_{y,m}}(\eta)] \cdot [L_{p_{x,n}}(\xi)L_{p_{y,n}}(\eta)] d\xi d\eta$$

由于 Legendre 多项式的正交性 $\int_{-1}^1 L_i(\xi)L_j(\xi)d\xi = \frac{2}{2i+1}\delta_{ij}$，质量矩阵为**对角阵**：

$$M_{mm} = \frac{h^2}{4} \cdot \frac{2}{2p_{x,m}+1} \cdot \frac{2}{2p_{y,m}+1} = \frac{h^2}{(2p_{x,m}+1)(2p_{y,m}+1)}$$

其逆矩阵系数为：

$$M_{mm}^{-1} = \frac{(2p_{x,m}+1)(2p_{y,m}+1)}{h^2}$$

代码中的 `invMass(m)` 函数对应数学公式：
```cpp
static inline double invMass(int m) noexcept {
    return (double)((2*mpx(m)+1) * (2*mpy(m)+1));
}
```

#### 数值积分

Gauss-Legendre 求积用于计算弱形式右端的体积分和面积分。

**体积分（Volume Integral）**

通过分部积分处理通量项：

$$\int_K \mathbf{F} \cdot \nabla \varphi_m d\mathbf{x} + \int_K \mathbf{G} \cdot \nabla \varphi_m d\mathbf{x}$$

进行坐标变换和引入梯度关系 $\nabla = (\frac{2}{h}\frac{\partial}{\partial\xi}, \frac{2}{h}\frac{\partial}{\partial\eta})$，得到参考空间积分：

$$\int_K \mathbf{F} \cdot \nabla \varphi_m d\mathbf{x} = \frac{h}{2} \int_{-1}^1 \int_{-1}^1 \mathbf{F} \cdot \nabla_\xi \varphi_m \, d\xi d\eta$$

使用 $N_Q = P+1$ 个 Gauss-Legendre 点进行二维求积：

$$(\text{RHS}_{\text{vol}})_m \approx \frac{h}{2} \sum_{qa=0}^{N_Q-1} \sum_{qb=0}^{N_Q-1} w_a w_b \left[ F_m(qa,qb) \frac{\partial \varphi_m}{\partial \xi}\bigg|_{qa,qb} + G_m(qa,qb) \frac{\partial \varphi_m}{\partial \eta}\bigg|_{qa,qb} \right]$$

其中 $w_a, w_b$ 为 Gauss 权重，$F_m, G_m$ 为在求积点处的通量值。代码实现对应：

```cpp
for (int qa=0; qa<NQ; ++qa) {
    for (int qb=0; qb<NQ; ++qb) {
        double ww = GL_WT[qa] * GL_WT[qb];
        // 计算通量 F, G 和基函数梯度 dphi_dxi, dphi_deta
        rhs[c][m] += ww * (F[c]*dphi_dxi + G[c]*dphi_deta);
    }
}
```

**面积分（Surface Integral）**

边界项在参考元的四条边上离散。以右边界（$\xi=1$）为例：

$$-\oint_{\text{right}} (F^* \cdot \mathbf{n}) \varphi_m \, ds = -\frac{h}{2} \int_{-1}^1 F^*(\xi=1, \eta) \varphi_m(1, \eta) \, d\eta$$

使用 $N_Q$ 点 Gauss 积分：

$$(\text{RHS}_{\text{edge}})_m \approx -\frac{h}{2} \sum_{q=0}^{N_Q-1} w_q \, F^*_q \varphi_m(1, \eta_q)$$

代码实现对应四条边：

```cpp
// 右边界 (xi = +1)
for (int q=0; q<NQ; ++q) {
    double phi_val = phi2(m, +1.0, GL_PT[q]);
    rhs[c][m] -= phi_val * GL_WT[q] * Fstar_right[c];
}
// 左边界 (xi = -1)
for (int q=0; q<NQ; ++q) {
    double phi_val = phi2(m, -1.0, GL_PT[q]);
    rhs[c][m] += phi_val * GL_WT[q] * Fstar_left[c];
}
// 类似处理上下边界（y 方向）
```

**求积精度**

代码使用 $N_Q = P+1$ 个求积点：
- 体积分（二维）：精确至 $2P+1$ 次多项式
- 面积分（一维）：精确至 $P$ 次多项式

对于 DG(P=2) ，使用 $N_Q=3$ 点 Gauss 积分可精确积分 5 次多项式，足以处理线性/二次通量近似。

#### 最终系数更新公式

综合质量矩阵逆、体积分和面积分，模态系数的演化方程为：

$$\frac{d\hat{\mathbf{u}}_m}{dt} = M_{mm}^{-1} \left[ (\text{RHS}_{\text{vol}})_m + (\text{RHS}_{\text{surf}})_m \right]$$

代入前述结果：

$$\frac{d\hat{\mathbf{u}}_m}{dt} = \frac{(2p_x+1)(2p_y+1)}{h^2} \left[ \frac{h}{2} \cdot (\text{rhs}_{\text{code}}) \right] = (\text{rhs}_{\text{code}}) \cdot \underbrace{\frac{(2p_x+1)(2p_y+1)}{2h}}_{\text{Scaling factor}}$$

代码中的完整实现为：

```cpp
const double inv2h = 0.5 / h;  // 对应公式中的 1/(2h)

for (int m=0; m<N_modes; ++m) {
    double scale = invMass(m) * inv2h;  // (2px+1)*(2py+1) / (2h)
    dU[c][m] = rhs[c][m] * scale;
}
```

这样保证了从参考空间计算的 `rhs` 值正确转换回物理空间的时间导数。

#### 数值稳定性

为了处理激波捕捉和极值区域的数值稳定性，本程序实现了双重限制机制：

##### 1. Persson 激波检测器 (Artificial Viscosity / Modal Decay)
基于 Legendre 多项式的模态能量谱（Modal Energy Spectrum）来识别解的间断性。
- **原理**：检查最高阶模态能量 $E_p$ 与总能量 $E_{tot}$ 的占比。若占比超过阈值 $s_0$，说明高阶项存在 Gibbs 振荡，解不再光滑。
- **平滑策略**：引入自适应算子 $\epsilon$，对最高阶模态进行衰减：
  $$s = \log_{10} \left( \frac{E_p}{E_{tot} + \epsilon} \right)$$
- **代码实现**：参考 `applyLimiter()` 中的 `Persson shock detector` 部分。



##### 2. Zhang-Shu 保持正性限制器 (Positivity-Preserving)
在极低密度或极高马赫数流场中，数值振荡可能导致负密度或负压力，引发模拟崩溃。
- **实现步骤**：
  1. **单元平均值强制正性**：确保 $\bar{\rho} > \rho_{min}$ 且 $\bar{p} > p_{min}$。
  2. **密度限制**：在 Gauss-Legendre 求积点处检查 $\rho_h(x_{G})$，若小于阈值，则将高阶模态向平均值压缩（Scaling）：
     $$\theta_\rho = \min \left( 1, \frac{\bar{\rho} - \rho_{min}}{\bar{\rho} - \rho_{min, G}} \right)$$
  3. **压力限制**：在保证密度的前提下，对压力进行类似的压缩变换，确保全场 $p > 0$。
- **数学特性**：该限制器是 **保守的 (Conservative)** 且 **高阶不失真 (High-order accurate)**，因为它不改变单元平均值。

##### 3. 边界处理与 Ghost Cells
- 采用 3 层 Ghost Cells 处理对称、周期或流出边界。
- 在 `applyLimiter` 循环中，自动跳过边界单元（`bctype != 0`），由专门的边界条件函数更新。

---

##### 实现参考函数
| 模块 | 函数名 | 位置 |
| :--- | :--- | :--- |
| **激波检测** | `applyLimiter()` | `src_DG/fluid.cpp` |
| **正性修正** | `Zhang-Shu algorithm` | `src_DG/fluid.cpp` |
| **基函数求值** | `basis(m, xi, eta)` | `src_DG/basis.cpp` |
#### CFL 稳定条件

$$\Delta t_{\max} = \frac{C}{2(2P+1)} \cdot \frac{h}{\lambda_{\max}}, \quad C = 0.9$$

其中 $\lambda_{\max}$ 为最大特征速度，$h$ 为网格尺寸。

---

## 系统架构

### 并行策略

采用**列方向域分解**：沿 $x$ 方向将全局网格均匀切分为 $n_{\text{proc}}$ 个子域，每个 MPI 进程负责一个子域。

**Ghost 层配置**：
- **WENO5**：每侧 3 层 ghost（支持 6 点重构模板）
- **DG**：每侧 1 层 ghost（多项式界面求值）

**通信模式**：SSP-RK3 每个阶段后进行一次 MPI 通信，交换所有 4 个守恒变量的 ghost 层。

### 文件结构

```
ZuikakuCFD/
├── src_WENO/              # 有限体积 WENO5 实现
│   ├── fluid.h            # 网格数据结构与 API
│   ├── fluid.cpp          # 重构、通量、时间推进
│   └── solver.cpp         # MPI 主程序、时间循环
├── src_DG/                # DG 间断伽辽金实现
│   ├── fluid.h            # DG 网格扩展（模态系数存储）
│   ├── fluid.cpp          # Legendre 基、弱形式、限制器
│   └── solver.cpp         # DG 主程序、自适应 CFL
├── 2D-Riemann/            # 标准测试网格数据
├── HighMachJet/           # 天体物理喷流算例
├── img/                   # 结果展示图片
│   ├── HMJP.png           # 喷流压力分布
│   ├── HMJV.png           # 喷流速度幅值
│   ├── WENO5.png          # Riemann 问题 WENO5 结果
│   ├── DG2.png            # Riemann 问题 DG(P=2) 结果
│   └── DG3.png            # Riemann 问题 DG(P=3) 结果
├── Makefile               # 构建脚本
├── gen.ipynb              # 网格生成示例
├── plot.ipynb             # 结果后处理脚本
├── LICENSE                # MIT 许可证
└── README.md              # 本文件
```

---

## 快速开始

### 系统要求

- **C++ 编译器**：GCC 7+、Clang 5+ 或 MSVC 2017+（需 C++17 支持）
- **MPI**：OpenMPI 3.0+ 或 MPICH 3.0+（需 `mpic++` 可用）
- **Eigen3**：仅需头文件（通常已在系统库中）

**验证环境**：
```bash
which mpic++         # 可用性检查
mpic++ --version     # 显示 MPI 版本
```

### 编译

```bash
# 编译所有版本（WENO5 + DG P=1,2,3）
make all

# 仅编译 WENO5
make weno

# 仅编译 DG（生成 solver_DG2, solver_DG3, solver_DG4）
make dg

# 清理
make clean          # 删除中间文件和可执行文件
make distclean       # 同上 + 清理输出
```

### 基本运行

#### WENO5 有限体积格式
```bash
# 固定时间步、指定步数
mpirun -np 4 ./solver_WENO ./2D-Riemann 1e-4 1000
```

参数说明：网格文件夹、时间步 $\Delta t$、总步数

#### DG 间断伽辽金格式
```bash
# 自适应 CFL、指定步数
mpirun -np 4 ./solver_DG3 ./2D-Riemann auto 1000
```

参数说明：网格文件夹、`auto`（自动CFL）或时间步值、总步数

### 可视化

使用 Jupyter notebook 后处理：

```bash
jupyter notebook plot.ipynb
```

脚本将加载 `result/` 目录中的并行输出数据，拼接为全局场景，计算原始变量和 Mach 数，绘制等高线图。

---

## 详细文档

### 网格格式

#### 参数文件 `params.txt`
```
nx ny da gamma
```

| 字段 | 含义 | 类型 | 示例 |
|------|------|------|------|
| nx | x 方向单元数 | int | 512 |
| ny | y 方向单元数 | int | 512 |
| da | 均匀网格尺寸 | double | 0.002 |
| gamma | 比热比 | double | 1.4 |

#### 初值文件（纯文本矩阵格式，行优先，尺寸 `ny × nx`）

| 文件 | 物理量 | 单位 |
|------|--------|------|
| `rho.dat` | 密度 $\rho$ | kg/m³ |
| `u.dat` | x 速度 $u$ | m/s |
| `v.dat` | y 速度 $v$ | m/s |
| `p.dat` | 压力 $p$ | Pa |
| `bctype.dat` | 边界标记（整数） | — |

#### 边界标记编码

```
  0 : 内部单元（参与计算）
 -1 : 固壁 ghost 单元（镜像复制）
 -3 : MPI 间 ghost 单元（通信填充）
```

### 数值精度与参数选择

| 方法 | 阶数 | ghost层 | 模式/单元 | CFL上限 | 推荐场景 |
|------|-----|--------|---------|--------|---------|
| MUSCL | 2 | 3 | 1 | 0.4 | 低成本快速测试 |
| **WENO5** | **5** | **3** | **1** | **0.4** | **强激波、工业应用** |
| DG(P=1) | 2 | 1 | 4 | 0.167 | 激波主导问题 |
| **DG(P=2)** | **3** | **1** | **9** | **0.100** | **均衡精度/成本** |
| DG(P=3) | 4 | 1 | 16 | 0.067 | 高精度光滑场 |

### 性能优化

#### WENO5 版本
- 栈分配 `double[4]` 替代 `vector`，消除堆分配开销
- Ghost 填充采用两遍扫描避免 O(N²) bctype 检测
- OpenMP 并行化内层循环（行向量化）
- 关键函数（WENO5、HLLC）inline 优化

#### DG 版本
- Legendre 多项式值预计算与缓存
- 所有 4 个变量×所有模式一次 MPI 通信打包
- 恢复原始变量与最大波速计算合并循环
- 单层 ghost 通信量相比 WENO5 减少 67%

---

## 典型测试用例

### 用例一：二维黎曼问题（Riemann Problem）

#### 问题设置

标准四象限黎曼问题，全局网格 **512×512**，四个象限独立初始化：

| 象限 | 位置 | $\rho$ | $u$ | $v$ | $p$ |
|------|------|--------|--------|---------|---------|
| 右上 | $x \geq 0.5, y \geq 0.5$ | 1.5 | 0.0 | 0.0 | 1.5 |
| 左上 | $x < 0.5, y \geq 0.5$ | 0.5323 | 1.206 | 0.0 | 0.3 |
| 左下 | $x < 0.5, y < 0.5$ | 0.138 | 1.206 | -1.206 | 0.029 |
| 右下 | $x \geq 0.5, y < 0.5$ | 0.5323 | 0.0 | -1.206 | 0.3 |

#### 物理特征

中心的四条间断线各形成激波、膨胀波或接触间断的复杂相互作用。是验证**多维激波捕捉、格点对称性、间断分辨率**的经典 test case。

#### 运行与对比

```bash
make all

# WENO5（5阶精度）
mpirun -np 4 ./solver_WENO ./2D-Riemann time 0.8

# DG(P=2)（3阶精度）
mpirun -np 4 ./solver_DG3 ./2D-Riemann time 0.8

# DG(P=3)（4阶精度）
mpirun -np 4 ./solver_DG4 ./2D-Riemann time 0.8

# 后处理
jupyter notebook plot.ipynb
```

#### 结果展示

| WENO5 | DG(P=2) | DG(P=3) |
|-------|---------|---------|
| ![WENO5.png](img/WENO5.png) | ![DG2.png](img/DG2.png) | ![DG3.png](img/DG3.png) |

**特点**：
- **WENO5**：激波锐利、无振荡，但需 3 层 ghost 通信
- **DG(P=2)**：光滑且均衡，自动限制器消除 Gibbs 振荡
- **DG(P=3)**：高精度刻画细化结构，但时间步最严格

---

### 用例二：天体物理超高马赫数喷流

#### 问题描述

极音速喷流（ Mach ≈ 2700）从左边界射入低速背景，激发复杂的激波-膨胀-接触间断相互作用。

#### 计算参数

| 参数 | 值 | 说明 |
|------|-----|------|
| 计算域 | $[0, 1] \times [-0.25, 0.25]$ | 长宽比 |
| 喷流出口 | $x=0, y \in [-0.05, 0.05]$ | 左边界条件 |
| 喷流状态 | $\rho=5, u=800, v=0, p=0.4127$ | 出口条件 |
| 背景状态 | $\rho=0.5, u=0, v=0, p=0.4127$ | 初值 |
| 比热比 | $\gamma = 5/3$ | 单原子气体 |
| 网格 | $400 \times 200$ | 推荐 |

#### 求解与结果

```bash
# DG P=2 求解
mpirun -np 4 ./solver_DG2 HighMachJet 1e-8 100000

```

| 压力分布（DG P=2） | 速度幅值（DG P=2） |
|----------|----------|
| ![HMJP.png](img/HMJP.png) | ![HMJV.png](img/HMJV.png) |

#### 物理现象

1. **喷流头部激波**：极高压力比、密度比
2. **膨胀扇**：接触面两侧的膨胀波
3. **接触间断**：密度跳跃但压力连续
4. **剪切层涡量**：速度剪切诱发的集中涡量
5. **数值鲁棒性考验**：极端流动参数对激波捕捉精度的极致要求

---

## 常见问题与故障排除

**Q: 何时选用 WENO5 vs DG？** 
A: WENO5 适合强激波问题（如爆炸、撞击）；DG 适合光滑流或多项式适配（如层流、稀薄气体）。在精度相当时，DG 通信开销更小。

**Q: DG 为什么需要限制器？**
A: DG 多项式在间断处通常产生 Gibbs 振荡。Cockburn-Shu 限制器在检测到陡峭梯度时自动激活，强制降阶以消除振荡。

**Q: 如何改变 DG 精度？**
编辑 [src_DG/fluid.h](src_DG/fluid.h) 并重新编译：
```cpp
static constexpr int DG_P = 2;  // 改为 1, 2 或 3
```
```bash
make distclean && make dg
```

**Q: 结果出现 NaN 或数值爆炸？**
- **WENO5**：检查 CFL 数是否过大（应 ≤ 0.4），减小 `dt`
- **DG**：若手动指定 `dt`，需满足 CFL 条件；建议使用 `auto` 模式

**Q: 不同进程数结果不一致？**
A: 浮点数舍入誤差在并行中累积。网格切分点处应接近，但不会完全相同。若误差 > 1%，检查是否有 race condition 或通信错误。

---

## 许可与引用

**许可证**：MIT License（见 [LICENSE](LICENSE)）

**引用**：若在学术出版中使用本代码，请参考以下关键文献：

1. **DG 方法与限制器**：
   - Cockburn, B., & Shu, C. W. (2001). *Runge–Kutta Discontinuous Galerkin Methods for Convection-Dominated Problems*. Journal of Scientific Computing, 16(3), 173-261.

2. **WENO 重构**：
   - Jiang, G. S., & Shu, C. W. (1996). *Efficient Implementation of Weighted ENO Schemes*. Journal of Computational Physics, 126(1), 202-228.

3. **Riemann 求解器**：
   - Toro, E. F. (2009). *Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction* (3rd ed.). Springer-Verlag.

4. **激波捕捉与限制器**：
   - Barth, T. J., & Jespersen, D. C. (1989). *The Design and Application of Upwind Schemes on Unstructured Meshes*. AIAA Paper 89-0366.

---

## 扩展方向

可进一步开发的方向包括：

- 二维同时域分解（减少通信边界）
- 自适应网格细化（AMR）
- 反应涡燃烧模型耦合
- GPU 加速（CUDA/HIP）
- 更高阶 DG（P≥4）
- 稀有气体 BGK 碰撞模型
- 磁流体动力学（MHD）扩展

---

## 致谢与支持

- 基于经典有限体积和间断伽辽金理论
- 利用 Eigen 高性能线性代数库
- OpenMPI/MPICH 提供稳定的并行基础
- 感谢哈工大威海周军伟老师提供的指导
---

**项目信息**  
- **作者**：midway酱  
- **创建时间**：2024  
- **最后更新**：2026年3月  
- **版本**：1.0
