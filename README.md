# 3D Wave Equation Solver — OpenMP, MPI & Hybrid

**English** | [简体中文](README.zh-CN.md)

A C++ implementation of a three-dimensional wave equation solver, developed for the **2025 Supercomputer Modeling and Technologies course** at the Faculty of Computational Mathematics and Cybernetics (CMC), Lomonosov Moscow State University.

The project implements **Variant 7** using an explicit finite-difference scheme and compares three parallel programming models: shared-memory **OpenMP**, distributed-memory **MPI**, and **MPI + OpenMP**. It includes analytical-solution validation and batch experiment scripts for the **Polus / IBM POWER8** cluster.

## Highlights

- Second-order central differences in space and time with a seven-point spatial stencil.
- Periodic boundaries in **x/y** and homogeneous Dirichlet boundaries in **z**.
- Three rolling time levels in contiguous arrays; pointer swaps avoid copying whole grids.
- MPI 3D Cartesian decomposition, including uneven grid partitions, with nonblocking halo exchange and interior computation before communication completion.
- Hybrid implementation requests `MPI_THREAD_FUNNELED`, with MPI calls on the main thread and OpenMP work sharing within each process.
- CFL reporting, final maximum absolute error, and calculation timing.
- LSF job scripts and hybrid result collection for reproducible cluster experiments.

## Mathematical model

The solver approximates

$$
\frac{\partial^2 u}{\partial t^2}=a^2\left(\frac{\partial^2 u}{\partial x^2}+\frac{\partial^2 u}{\partial y^2}+\frac{\partial^2 u}{\partial z^2}\right),
\qquad a^2=\frac{1}{4/L_x^2+4/L_y^2+1/L_z^2}.
$$

The analytical reference solution is

$$
u(t,x,y,z)=\sin\left(\frac{2\pi x}{L_x}+3\pi\right)\sin\left(\frac{2\pi y}{L_y}+2\pi\right)\sin\left(\frac{\pi z}{L_z}\right)\cos(\pi t+\pi).
$$

The initial displacement is this solution at `t = 0`; initial velocity is zero. The first time level is initialized using a Taylor expansion, followed by explicit three-level time stepping.

Current source defaults fix `Lx = Ly = Lz = 1` and final time `T = 0.01`, giving `a² = 1/9`. For `N³` grid points:

- `dx = Lx/N`, `dy = Ly/N`: periodic endpoints are not duplicated.
- `dz = Lz/(N−1)`: both Dirichlet endpoints are included.
- `dt = T/steps`.

The reported stability condition is

$$
a\,\Delta t\sqrt{\Delta x^{-2}+\Delta y^{-2}+\Delta z^{-2}}\leq 1.
$$

**The program only warns when this condition fails; it does not stop.** Increase the number of time steps when needed.

## Repository layout

| Path | Purpose |
| --- | --- |
| `wave_eq_variant7_omp.cpp` | OpenMP solver |
| `wave_eq_variant7_mpi.cpp` | MPI solver with a single-process serial path |
| `wave_eq_variant7_hybrid.cpp` | MPI + OpenMP solver; one MPI process still uses OpenMP |
| `makefile` | IBM XL / Spectrum MPI build targets for Polus |
| `script_polus_omp_version/` | OpenMP LSF experiment script |
| `script_polus_mpi_version/` | MPI LSF experiment script |
| `script_polus_hybrid_version/` | Hybrid jobs, submission helper, and result collector |
| `smit_step1.pdf`, `smit_step2.pdf`, `smit_step3.pdf` | Course assignment materials (Russian) |
| `итоговый_отчет(edu-cmc-sqi22-18).pdf` | Original project report (Russian) |
| `README.txt` | Original project introduction (Russian) |

## Build

### Polus: IBM XL and Spectrum MPI

The supplied makefile targets POWER8 and loads `SpectrumMPI/10.1.0`. Run it in a cluster environment where the `module` command is available to Bash:

```bash
make all
# Or build one implementation:
make omp
make mpi
make hybrid
```

Executables are placed in their corresponding script directories:

```text
script_polus_omp_version/wave_eq_variant7
script_polus_mpi_version/wave_mpi
script_polus_hybrid_version/wave_hybrid
```

`make clean` removes these three executables.

### Linux workstation: GCC and an MPI implementation

Use a C++11-capable compiler, OpenMP support, and (for MPI/hybrid) an MPI development installation providing `mpicxx` and `mpirun`. The hybrid source uses OpenMP SIMD directives; GCC with OpenMP 4.0 or newer support is appropriate. The MPI implementations also use POSIX aligned allocation.

The default makefile is cluster-specific. Compile directly for a local build:

```bash
g++ -O3 -std=c++11 -fopenmp wave_eq_variant7_omp.cpp -o wave_eq_variant7
mpicxx -O3 -std=c++11 wave_eq_variant7_mpi.cpp -o wave_mpi
mpicxx -O3 -std=c++11 -fopenmp wave_eq_variant7_hybrid.cpp -o wave_hybrid
```

## Run

The following examples use local-build executables in the repository root:

```bash
# OpenMP: grid size, time steps, thread count
./wave_eq_variant7 128 100 4

# MPI: four processes; grid size, time steps
mpirun -np 4 ./wave_mpi 128 100

# Hybrid: two MPI processes, two OpenMP threads per process
OMP_NUM_THREADS=2 mpirun -np 2 ./wave_hybrid 128 100
```

| Setting | OpenMP | MPI | Hybrid |
| --- | --- | --- | --- |
| First argument: `N` | Default `128` | Default `256` | Default `256` |
| Second argument: time steps | Default `20` | Default `20` | Default `20` |
| Thread count | Third argument, default `4` | Not applicable | `OMP_NUM_THREADS` / runtime default |
| MPI process count | Not applicable | MPI launcher | MPI launcher |

The standalone OpenMP program explicitly sets its thread count from the third argument; setting only `OMP_NUM_THREADS` does not override its default of four. The hybrid program has **no third thread-count argument**.

Use integers with `N >= 2`, `steps >= 1`, and positive thread/process counts. Input validation is currently limited. Avoid process grids that create empty subdomains; use comfortably sized local blocks for initial experiments. Process-grid dimensions are selected automatically by `MPI_Dims_create`.

### Small validation run

```bash
./wave_eq_variant7 32 100 2
```

A GCC/OpenMP smoke test produced:

```text
3D CFL number = 0.00182848 (should be <= 1 for stability)
Final max error: 1.452267e-06
```

Timing depends on hardware and load. This is a small correctness smoke test, **not a Polus benchmark or validation of the MPI/hybrid executables**. Compare the final errors of all three implementations with identical `N` and `steps` on your target system.

## Polus batch experiments

Submit from the directory containing the executable and job script:

```bash
(cd script_polus_omp_version && bsub < wave.lsf)
(cd script_polus_mpi_version && bsub < wavempi.lsf)
(cd script_polus_hybrid_version && bsub < hybrid_128_4_1.lsf)
```

Before submitting, review queue/host constraints, wall time, core allocation, executable paths, and the MPI module environment. These scripts contain site-specific settings, not portable scheduler defaults.

For all supplied hybrid configurations:

1. Edit `LSF_DIR` in `script_polus_hybrid_version/submit_all_jobs.sh` to your actual hybrid script directory.
2. From that directory, run `bash submit_all_jobs.sh`.
3. After jobs finish, run `bash collect_results_clean.sh > results.txt` in the same directory.

The collector reads `run_*.out`, needs Bash associative arrays and `bc`, and filters timings using a median-based threshold (50% by default). Keep raw logs alongside cleaned summaries. The supplied hybrid jobs cover `N = 128, 256, 512` and multiple process/thread combinations.

## Interpreting results and current limitations

- **Error:** `Final max error` is the maximum absolute difference from the analytical solution at `T`, not a relative error or residual. Distributed runs reduce local maxima with `MPI_MAX`.
- **Timing:** calculation timings cover the main time-stepping loop, excluding initialization/first-step setup and final error evaluation. MPI multi-process runs report the slowest rank's wall time; the MPI-only single-process path uses CPU time (`clock()`), while OpenMP/hybrid paths use wall-clock timers. Account for this when comparing baselines.
- **Memory:** three double-precision time levels need about `24 × N³` bytes before MPI face buffers and runtime overhead—approximately **3 GiB at N = 512** across the grid. Do not start with a large grid on a small workstation.
- **Benchmark comparability:** supplied scripts do not all use the same time-step count. Align problem size, steps, compiler flags, process/thread placement, and timing definitions before comparing performance.
- **Collector efficiency:** the hybrid collector chooses the first available configuration as its baseline, which may already be parallel. Its printed efficiency divides speedup by the target worker count without normalizing by baseline workers. Recalculate relative efficiency as `(t_base/t_run) × (workers_base/workers_run)` when appropriate; do not interpret that column as serial-baseline efficiency.
- **Thread support:** the hybrid code warns but continues if the MPI library provides less than `MPI_THREAD_FUNNELED`. Use an MPI runtime that satisfies the requested level.
- **Scope:** this is a course solver and benchmarking project. It prints summary results; it does not provide field-file export, visualization, or a general-purpose PDE configuration interface.

## Credits and project materials

- **Author:** Чжан Синьюй (Zhang Xinyu), group 623.
- **Instructor:** Жуков К. А. (K. A. Zhukov).
- **Institution:** Faculty of CMC, Lomonosov Moscow State University.
- **Course:** Supercomputer Modeling and Technologies, 2025.

See the assignment PDFs and original report in this repository for the course materials and submitted analysis. No license file is currently included; this README does not grant additional reuse rights.
