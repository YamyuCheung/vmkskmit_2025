Реализация седьмого практического задания в рамках курса "Суперкомпьютерное моделирование и технологии" 2025 г. ВМК МГУ


## Сборка и Запуск (Build and Run)

For convenience, a `Makefile` is provided to automate the compilation process on the Polus cluster.

### Быстрая сборка (Quick Build)

To compile **all versions** (OpenMP, MPI, Hybrid) at once:
```bash
make all



### OpenMP-реализация
xlc++_r -std=c++11 -O3 -qarch=pwr8 -qtune=pwr8 -qsmp=omp \
        -qhot -qipa -o wave_eq_variant7 wave_eq_variant7.cpp

Then run the script in folder script_polus_omp_version

...

### MPI-реализация
mpixlC -O3 -qarch=pwr8 -qtune=pwr8 -qhot -qipa -std=c++11 wave_eq_variant7_mpi.cpp -o wave_mpi

Then run the script in folder script_polus_mpi_version


### MPI+OpenMP-реализация

mpixlC -O3 -qarch=pwr8 -qtune=pwr8 -qsmp=omp -qhot -qipa -std=c++11 wave_eq_variant7_hybrid.cpp -o wave_hybrid

Then run the scripts in folder script_polus_hybrid_version
...
