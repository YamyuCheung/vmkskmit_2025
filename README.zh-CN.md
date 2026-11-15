# 三维波动方程求解器 — OpenMP、MPI 与混合并行

[English](README.md) | **简体中文**

本项目是莫斯科国立大学计算数学与控制论系（CMC / ВМК）**2025 年《超级计算机建模与技术》课程**的实践项目，使用 C++ 实现三维波动方程的并行数值求解。

程序针对 **Variant 7（第 7 变体）**，采用显式有限差分方法，分别实现 **OpenMP 共享内存并行**、**MPI 分布式内存并行**和 **MPI + OpenMP 混合并行**，并提供解析解误差验证及面向 **Polus / IBM POWER8 集群**的批量实验脚本。

## 项目特点

- 空间与时间均采用二阶中心差分，空间离散使用七点模板。
- **x、y 方向周期边界**，**z 方向齐次 Dirichlet 边界**。
- 连续内存存储三个时间层，通过指针交换避免整块网格复制。
- MPI 使用三维笛卡尔进程拓扑，支持非均匀网格划分；采用非阻塞边界交换，在等待通信完成前计算内部网格点。
- 混合版本请求 `MPI_THREAD_FUNNELED`：主线程负责 MPI 调用，进程内部由 OpenMP 分担计算。
- 输出 CFL 数、最终最大绝对误差和计算耗时。
- 提供 LSF 作业脚本及混合版本结果汇总工具，方便开展集群实验。

## 数学模型

求解的三维波动方程为：

$$
\frac{\partial^2 u}{\partial t^2}=a^2\left(\frac{\partial^2 u}{\partial x^2}+\frac{\partial^2 u}{\partial y^2}+\frac{\partial^2 u}{\partial z^2}\right),
\qquad a^2=\frac{1}{4/L_x^2+4/L_y^2+1/L_z^2}.
$$

用于验证的解析解为：

$$
u(t,x,y,z)=\sin\left(\frac{2\pi x}{L_x}+3\pi\right)\sin\left(\frac{2\pi y}{L_y}+2\pi\right)\sin\left(\frac{\pi z}{L_z}\right)\cos(\pi t+\pi).
$$

初始位移取解析解在 `t = 0` 时的值，初始速度为零。程序通过 Taylor 展开初始化第一个时间层，再使用显式三层差分格式推进。

当前源码固定 `Lx = Ly = Lz = 1`、终止时刻 `T = 0.01`，因此 `a² = 1/9`。对于 `N³` 个网格点：

- `dx = Lx/N`、`dy = Ly/N`：周期方向不重复存储末端点。
- `dz = Lz/(N−1)`：包含两端的 Dirichlet 边界点。
- `dt = T/steps`。

程序检查的稳定性条件为：

$$
a\,\Delta t\sqrt{\Delta x^{-2}+\Delta y^{-2}+\Delta z^{-2}}\leq 1.
$$

**不满足条件时，程序只输出警告，不会停止运行。** 必要时应增加时间步数。

## 仓库结构

| 路径 | 用途 |
| --- | --- |
| `wave_eq_variant7_omp.cpp` | OpenMP 求解器 |
| `wave_eq_variant7_mpi.cpp` | MPI 求解器，包含单进程串行路径 |
| `wave_eq_variant7_hybrid.cpp` | MPI + OpenMP 求解器；单 MPI 进程时仍使用 OpenMP |
| `makefile` | 面向 Polus 的 IBM XL / Spectrum MPI 编译规则 |
| `script_polus_omp_version/` | OpenMP LSF 实验脚本 |
| `script_polus_mpi_version/` | MPI LSF 实验脚本 |
| `script_polus_hybrid_version/` | 混合作业、批量提交及结果汇总脚本 |
| `smit_step1.pdf`、`smit_step2.pdf`、`smit_step3.pdf` | 俄文课程任务材料 |
| `итоговый_отчет(edu-cmc-sqi22-18).pdf` | 俄文原始项目报告 |
| `README.txt` | 原始俄文项目介绍 |

## 编译

### Polus 集群：IBM XL 与 Spectrum MPI

仓库内的 makefile 面向 POWER8，并加载 `SpectrumMPI/10.1.0`。请在 Bash 能使用 `module` 命令的集群环境中执行：

```bash
make all
# 也可以单独编译：
make omp
make mpi
make hybrid
```

生成的可执行文件会放到对应脚本目录：

```text
script_polus_omp_version/wave_eq_variant7
script_polus_mpi_version/wave_mpi
script_polus_hybrid_version/wave_hybrid
```

`make clean` 会删除上述三个可执行文件。

### 普通 Linux 环境：GCC 与 MPI

需要支持 C++11 的编译器、OpenMP，以及 MPI／混合版本所需的 MPI 开发环境（提供 `mpicxx` 和 `mpirun`）。混合源码使用 OpenMP SIMD 指令，建议使用支持 OpenMP 4.0 或更新版本的 GCC。MPI 相关实现还使用 POSIX 对齐内存分配。

默认 makefile 专用于集群，本地可以直接编译：

```bash
g++ -O3 -std=c++11 -fopenmp wave_eq_variant7_omp.cpp -o wave_eq_variant7
mpicxx -O3 -std=c++11 wave_eq_variant7_mpi.cpp -o wave_mpi
mpicxx -O3 -std=c++11 -fopenmp wave_eq_variant7_hybrid.cpp -o wave_hybrid
```

## 运行

以下示例使用本地编译后位于仓库根目录的可执行文件：

```bash
# OpenMP：网格大小、时间步数、线程数
./wave_eq_variant7 128 100 4

# MPI：4 个进程；网格大小、时间步数
mpirun -np 4 ./wave_mpi 128 100

# 混合并行：2 个 MPI 进程，每个进程 2 个 OpenMP 线程
OMP_NUM_THREADS=2 mpirun -np 2 ./wave_hybrid 128 100
```

| 参数 | OpenMP | MPI | 混合版本 |
| --- | --- | --- | --- |
| 第一个参数：`N` | 默认 `128` | 默认 `256` | 默认 `256` |
| 第二个参数：时间步数 | 默认 `20` | 默认 `20` | 默认 `20` |
| 线程数 | 第三个参数，默认 `4` | 不适用 | `OMP_NUM_THREADS`／运行时默认值 |
| MPI 进程数 | 不适用 | 由 MPI 启动器指定 | 由 MPI 启动器指定 |

独立 OpenMP 程序会根据第三个参数显式设置线程数；仅设置 `OMP_NUM_THREADS` 不会覆盖其默认的 4 线程。混合版本**没有用于设置线程数的第三个命令行参数**。

请输入整数，并确保 `N >= 2`、`steps >= 1`，线程数和进程数均为正数。当前参数校验较少。应避免产生空子域的进程划分，初次实验建议保留足够大的局部网格。进程拓扑维度由 `MPI_Dims_create` 自动选择。

### 小规模验证

```bash
./wave_eq_variant7 32 100 2
```

一次 GCC/OpenMP 冒烟测试得到：

```text
3D CFL number = 0.00182848 (should be <= 1 for stability)
Final max error: 1.452267e-06
```

耗时随硬件和负载变化。该结果仅为小规模正确性冒烟测试，**不是 Polus 性能基准，也不代表 MPI／混合可执行文件已经通过实测**。请在目标环境使用相同的 `N` 和 `steps` 对比三个版本的最终误差。

## Polus 批量实验

从包含可执行文件和作业脚本的目录提交：

```bash
(cd script_polus_omp_version && bsub < wave.lsf)
(cd script_polus_mpi_version && bsub < wavempi.lsf)
(cd script_polus_hybrid_version && bsub < hybrid_128_4_1.lsf)
```

提交前请检查队列、主机限制、运行时限、核心分配、程序路径及 MPI 模块环境。脚本包含 Polus 专用设置，并非通用调度模板。

批量运行仓库提供的混合配置：

1. 将 `script_polus_hybrid_version/submit_all_jobs.sh` 中的 `LSF_DIR` 改为实际的混合脚本目录。
2. 在该目录运行 `bash submit_all_jobs.sh`。
3. 等待作业完成后，在同一目录运行 `bash collect_results_clean.sh > results.txt`。

汇总脚本读取 `run_*.out`，需要支持关联数组的 Bash 和 `bc`，默认按相对中位数偏差超过 50% 的规则筛选耗时异常值。请同时保留原始日志和清洗后的汇总。现有混合作业覆盖 `N = 128、256、512` 和多种进程／线程组合。

## 结果解读与当前限制

- **误差含义：** `Final max error` 是终止时刻 `T` 的数值解与解析解之差的最大绝对值，不是相对误差或残差。分布式计算通过 `MPI_MAX` 汇总局部最大误差。
- **计时范围：** 计算耗时覆盖主要时间推进循环，不包含初始化／第一步计算和最终误差评估。MPI 多进程版本报告最慢进程的墙钟时间；纯 MPI 的单进程路径使用 `clock()` CPU 时间，OpenMP／混合路径使用墙钟计时。比较基线时应考虑这一差异。
- **内存需求：** 三层双精度网格约需 `24 × N³` 字节，不含 MPI 边界缓冲区和运行时开销；`N = 512` 时，全网格合计约 **3 GiB**。小内存机器请先运行小规模测试。
- **实验可比性：** 仓库各脚本的时间步数并不完全一致。比较性能前应统一网格、步数、编译选项、进程／线程绑定和计时口径。
- **汇总效率指标：** 混合汇总脚本选取首个可用配置作为基线，该基线可能已经是并行配置；输出的效率直接将加速比除以目标并行单元数，未按基线并行单元数归一化。相对效率应按 `(t_base/t_run) × (workers_base/workers_run)` 重新计算，不宜将现有列直接理解为相对串行的并行效率。
- **线程支持：** 当 MPI 实现提供的线程级别低于 `MPI_THREAD_FUNNELED` 时，混合代码只警告而继续运行。请使用满足所需线程支持级别的 MPI 环境。
- **项目范围：** 本项目面向课程数值求解和性能实验，输出摘要结果；暂不提供场数据文件导出、可视化或通用 PDE 配置接口。

## 作者与课程材料

- **作者：** Чжан Синьюй（Zhang Xinyu），623 组。
- **指导教师：** Жуков К. А.（K. A. Zhukov）。
- **院系：** 莫斯科国立大学计算数学与控制论系。
- **课程：** 超级计算机建模与技术，2025 年。

课程材料及提交的分析请参阅仓库内的任务 PDF 和原始项目报告。仓库目前未附带许可证文件；本 README 不额外授予使用或再分发许可。
