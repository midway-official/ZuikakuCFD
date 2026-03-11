# ZuikakuCFD — 2D Euler 方程并行求解器

本项目是一个基于 MPI + Eigen 的二维欧拉方程有限体积求解器，支持在多进程（列方向域分解）上并行求解，使用 HLLC 通量与 MUSCL 结构保持（MUSCL）格式。

---

## ✅ 特性

- MPI 并行：按列（x 方向）分域并交换守恒变量的 ghost 列
- 高精度：采用 MUSCL 重构 + HLLC 通量
- RK2 时间推进（SSP-RK2）
- 支持从磁盘读取网格与初始条件数据（`mesh/params.txt` + `.dat` 文件）
- 周期性输出结果到 `result/<step>/` 目录

---

## 📦 目录结构

```
./
├── build/            # 构建输出目录（目标文件、依赖文件）
├── mesh/             # 示例输入网格（参数 + 数据）
│   └── params.txt
├── report/           # 生成的分析/可视化文件（构建时创建）
├── src/              # 源代码
│   ├── fluid.cpp
│   ├── fluid.h
│   └── solver.cpp
├── gen.ipynb         # （可选）网格生成 / 后处理 Jupyter Notebook
├── plot.ipynb        # （可选）后处理绘图 Notebook
└── Makefile
```

---

## 🧩 依赖 (Requirements)

- C++17 兼容编译器
- MPI：`mpic++` / `mpirun`
- Eigen3（头文件即可）

> ⚠️ 请确保系统已安装 MPI（OpenMPI、MPICH 等）并可从命令行使用 `mpic++`。

---

## 🛠️ 构建（Build）

在项目根目录下运行：

```sh
make
```

构建成功后，会生成可执行文件 `solver`，以及 `build/` 目录下的中间文件。

---

## ▶️ 运行（Run）

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

## 🗂️ 网格数据格式说明

`mesh/params.txt` 结构（4 个数值）：

```
nx ny da gamma
```

其中：

- `nx`, `ny`：格点/单元数
- `da`：网格间距（单元边长）
- `gamma`：比热比

对应的 `.dat` 文件为纯文本矩阵，行优先存储，行数为 `ny`，列数为 `nx`：

- `rho.dat`：密度
- `u.dat`：x 方向速度
- `v.dat`：y 方向速度
- `p.dat`：压力
- `bctype.dat`：边界类型（整数编码）

---

## 📤 输出结果

程序会在每 100 步时将结果写入：

```
result/<step>/
```

每个目录包含当前子域进程对应的 `rank` 输出（例如 `rank_0.bin`、`rank_1.bin` 等），保存的是守恒变量或原始变量数据。

---

## 🚀 后处理与可视化

可使用 `plot.ipynb` 或其他常用工具（Python/Matplotlib、ParaView 等）读取输出结果并可视化。

---

## 📌 注意

- 目前网格分域仅针对 x 方向（列方向）分割，需确保 `nx` 大于进程数并能均匀分配。
- 输出目录会在 Rank 0 中创建，其他 Rank 通过 `MPI_Barrier` 同步。

---

## 🧪 进一步扩展

你可以考虑加入：

- 更灵活的网格划分（2D 域分解）
- 自适应时间步长（CFL 控制）
- 向量化/线程优化（OpenMP）、GPU 加速
- 更复杂的边界条件与源项
