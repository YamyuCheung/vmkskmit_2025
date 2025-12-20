#include <iostream>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <mpi.h>
#include <cstring>
#include <cstdlib>

using namespace std;

const double PI = 3.14159265358979323846;


// 解析解函数
// Функция аналитического решения
double u_analytical(double t, double x, double y, double z,
                    double Lx, double Ly, double Lz) {
    double at = PI;  // 变体7中 at = π
                     // В варианте 7 at = π
    return sin(2*PI/Lx * x + 3*PI) *
           sin(2*PI/Ly * y + 2*PI) *
           sin(PI/Lz * z) *
           cos(at * t + PI);
}

// 初始条件 u(0, x, y, z) = φ(x, y, z)
// Начальное условие u(0, x, y, z) = φ(x, y, z)
double phi(double x, double y, double z, double Lx, double Ly, double Lz) {
    return u_analytical(0, x, y, z, Lx, Ly, Lz);
}

// 初始速度 ∂u/∂t|_{t=0} = 0 (根据题目要求)
// Начальная скорость ∂u/∂t|_{t=0} = 0 (согласно условию задачи)
double phi_t(double x, double y, double z, double Lx, double Ly, double Lz) {
    return 0.0;
}

// 计算a²值
// Вычисление значения a²
double calculate_a_squared(double Lx, double Ly, double Lz) {
    return 1.0 / (4.0/(Lx*Lx) + 4.0/(Ly*Ly) + 1.0/(Lz*Lz));
}

// 1D 数组索引宏
// Макросы для индексации 1D массива
// (i, j, k) -> i * Ny*Nz + j * Nz + k
#define IDX_GLOBAL(i, j, k, count_y, count_z) ((size_t)(i) * (count_y) * (count_z) + (size_t)(j) * (count_z) + (k))
#define IDX_LOCAL(i, j, k, local_count_y, local_count_z) ((size_t)(i) * (local_count_y) * (local_count_z) + (size_t)(j) * (local_count_z) + (k))


// 计算误差
// Вычисление ошибки
double calculate_error(double *u, double t,
                       double delta_x, double delta_y, double delta_z,
                       int count_x, int count_y, int count_z,
                       double Lx, double Ly, double Lz, double a) {
    double max_error = 0.0;

    for(int i = 0; i < count_x; i++) {
        for(int j = 0; j < count_y; j++) {
            for(int k = 0; k < count_z; k++) {
                double x = i * delta_x;
                double y = j * delta_y;
                double z = k * delta_z;
                double exact = u_analytical(t, x, y, z, Lx, Ly, Lz);
                double error = fabs(u[IDX_GLOBAL(i, j, k, count_y, count_z)] - exact);
                if(error > max_error) max_error = error;
            }
        }
    }
    return max_error;
}

// 计算本地子区域的误差 (并行专用)
// Вычисление локальной ошибки (для параллельной версии)
double calculate_local_error(const double *local_u, double t,
                             double delta_x, double delta_y, double delta_z,
                             int local_count_x, int local_count_y, int local_count_z,
                             int start_i_global, int start_j_global, int start_k_global,
                             double Lx, double Ly, double Lz) {
    double max_error = 0.0;

    for(int i = 0; i < local_count_x; i++) {
        for(int j = 0; j < local_count_y; j++) {
            for(int k = 0; k < local_count_z; k++) {
                // 关键点：物理坐标 = (全局偏移量 + 本地索引) * 步长
                double x = (start_i_global + i) * delta_x;
                double y = (start_j_global + j) * delta_y;
                double z = (start_k_global + k) * delta_z;

                double exact = u_analytical(t, x, y, z, Lx, Ly, Lz);

                // 注意：这里使用 IDX_LOCAL 宏
                double val = local_u[IDX_LOCAL(i, j, k, local_count_y, local_count_z)];

                double error = fabs(val - exact);
                if(error > max_error) max_error = error;
            }
        }
    }
    return max_error;
}

// 检查稳定性条件（CFL条件）
// Проверка условия устойчивости (условие CFL)
bool check_stability(double a, double delta_x, double delta_y, double delta_z, double delta_t, int rank) {

    double cfl_3D = a * delta_t * sqrt( 1.0/(delta_x*delta_x) +
                                       1.0/(delta_y*delta_y) +
                                       1.0/(delta_z*delta_z) );
    if (rank == 0) {
        cout << "3D CFL number = " << cfl_3D << " (should be <= 1 for stability)" << endl;
    }
    return cfl_3D <= 1.0;
}


struct BorderFlags {
    bool isOnXMinusBorder; // X-Dirichlet (i=0)
    bool isOnXPlusBorder;  // X-Dirichlet (i=Lx)
    bool isOnYMinusBorder; // Y-Dirichlet (j=0)
    bool isOnYPlusBorder;  // Y-Dirichlet (j=Ly)
    bool isOnZMinusBorder; // Z-Dirichlet (k=0)
    bool isOnZPlusBorder;  // Z-Dirichlet (k=Lz)
};

void setBorderFlags(BorderFlags& flags, int xCoord, int yCoord, int zCoord, const int dims[3]) {
    flags.isOnXMinusBorder = (xCoord == 0);
    flags.isOnXPlusBorder  = (xCoord == dims[0] - 1);
    flags.isOnYMinusBorder = (yCoord == 0);
    flags.isOnYPlusBorder  = (yCoord == dims[1] - 1);
    flags.isOnZMinusBorder = (zCoord == 0);
    flags.isOnZPlusBorder  = (zCoord == dims[2] - 1);
}

// 提取 Y 軸負方向邊界面 (j=0) 的數據
// Извлечение данных граничной плоскости в отрицательном направлении оси Y (j=0)
void Get_Y_Face_Minus(double* mainData, double* faceData, int local_count_x, int local_count_y, int local_count_z) {
    for (int i = 0; i < local_count_x; ++i) {
        for (int k = 0; k < local_count_z; ++k) {
            faceData[i * local_count_z + k] = mainData[IDX_LOCAL(i, 0, k, local_count_y, local_count_z)];
        }
    }
}

// 提取 Y 軸正方向邊界面 (j=local_count_y-1) 的數據
// Извлечение данных граничной плоскости в положительном направлении оси Y (j=local_count_y-1)
void Get_Y_Face_Plus(double* mainData, double* faceData, int local_count_x, int local_count_y, int local_count_z) {
    for (int i = 0; i < local_count_x; ++i) {
        for (int k = 0; k < local_count_z; ++k) {
            faceData[i * local_count_z + k] = mainData[IDX_LOCAL(i, local_count_y - 1, k, local_count_y, local_count_z)];
        }
    }
}

// 提取 Z 軸負方向邊界面 (k=0) 的數據
// Извлечение данных граничной плоскости в отрицательном направлении оси Z (k=0)
void Get_Z_Face_Minus(double* mainData, double* faceData, int local_count_x, int local_count_y, int local_count_z) {
    for (int i = 0; i < local_count_x; ++i) {
        for (int j = 0; j < local_count_y; ++j) {
            faceData[i * local_count_y + j] = mainData[IDX_LOCAL(i, j, 0, local_count_y, local_count_z)];
        }
    }
}

// 提取 Z 軸正方向邊界面 (k=local_count_z-1) 的數據
// Извлечение данных граничной плоскости в положительном направлении оси Z (k=local_count_z-1)
void Get_Z_Face_Plus(double* mainData, double* faceData, int local_count_x, int local_count_y, int local_count_z) {
    for (int i = 0; i < local_count_x; ++i) {
        for (int j = 0; j < local_count_y; ++j) {
            faceData[i * local_count_y + j] = mainData[IDX_LOCAL(i, j, local_count_z - 1, local_count_y, local_count_z)];
        }
    }
}


// 计算 t=1 (u_curr) 的单个点
// Вычисление одной точки для t=1 (u_curr)
inline double CalculateStencil_t1(
    int i, int j, int k, // 本地索引
                         // Локальные индексы
    const double* local_u_prev,
    const double* recv_x_minus, const double* recv_x_plus,
    const double* recv_y_minus, const double* recv_y_plus,
    const double* recv_z_minus, const double* recv_z_plus,
    int local_count_x, int local_count_y, int local_count_z,
    double delta_x, double delta_y, double delta_z,
    double a_squared, double delta_t
) {
    const double u0_center = local_u_prev[IDX_LOCAL(i, j, k, local_count_y, local_count_z)];

    // X-邻居 (周期性)
    // X-соседи (периодические)
    const double u0_x_minus = (i == 0) ? recv_x_minus[j * local_count_z + k] : local_u_prev[IDX_LOCAL(i - 1, j, k, local_count_y, local_count_z)];
    const double u0_x_plus  = (i == local_count_x - 1) ? recv_x_plus[j * local_count_z + k] : local_u_prev[IDX_LOCAL(i + 1, j, k, local_count_y, local_count_z)];

    // Y-邻居 (周期性)
    // Y-соседи (периодические)
    const double u0_y_minus = (j == 0) ? recv_y_minus[i * local_count_z + k] : local_u_prev[IDX_LOCAL(i, j - 1, k, local_count_y, local_count_z)];
    const double u0_y_plus  = (j == local_count_y - 1) ? recv_y_plus[i * local_count_z + k] : local_u_prev[IDX_LOCAL(i, j + 1, k, local_count_y, local_count_z)];

    // Z-邻居 (非周期性)
    // Z-соседи (непериодические)
    const double u0_z_minus = (k == 0) ? recv_z_minus[i * local_count_y + j] : local_u_prev[IDX_LOCAL(i, j, k - 1, local_count_y, local_count_z)];
    const double u0_z_plus  = (k == local_count_z - 1) ? recv_z_plus[i * local_count_y + j] : local_u_prev[IDX_LOCAL(i, j, k + 1, local_count_y, local_count_z)];

    double laplacian_u0 = (u0_x_plus - 2*u0_center + u0_x_minus) / (delta_x * delta_x) +
                          (u0_y_plus - 2*u0_center + u0_y_minus) / (delta_y * delta_y) +
                          (u0_z_plus - 2*u0_center + u0_z_minus) / (delta_z * delta_z);

    return u0_center + 0.5 * a_squared * delta_t * delta_t * laplacian_u0;
}


// 计算 t > 1 (u_next) 的单个点
// Вычисление одной точки для t > 1 (u_next)
inline double CalculateStencil_t_n(
    int i, int j, int k, // 本地索引
                         // Локальные индексы
    const double* local_u_prev, // t-1
    const double* local_u_curr, // t
    const double* recv_x_minus, const double* recv_x_plus,
    const double* recv_y_minus, const double* recv_y_plus,
    const double* recv_z_minus, const double* recv_z_plus,
    int local_count_x, int local_count_y, int local_count_z,
    double r_x, double r_y, double r_z
) {
    const double u_center = local_u_curr[IDX_LOCAL(i, j, k, local_count_y, local_count_z)];

    // X-邻居 (周期性)
    // X-соседи (периодические)
    const double u_x_minus = (i == 0) ? recv_x_minus[j * local_count_z + k] : local_u_curr[IDX_LOCAL(i - 1, j, k, local_count_y, local_count_z)];
    const double u_x_plus  = (i == local_count_x - 1) ? recv_x_plus[j * local_count_z + k] : local_u_curr[IDX_LOCAL(i + 1, j, k, local_count_y, local_count_z)];

    // Y-邻居 (周期性)
    // Y-соседи (периодические)
    const double u_y_minus = (j == 0) ? recv_y_minus[i * local_count_z + k] : local_u_curr[IDX_LOCAL(i, j - 1, k, local_count_y, local_count_z)];
    const double u_y_plus  = (j == local_count_y - 1) ? recv_y_plus[i * local_count_z + k] : local_u_curr[IDX_LOCAL(i, j + 1, k, local_count_y, local_count_z)];

    // Z-邻居 (非周期性)
    // Z-соседи (непериодические)
    const double u_z_minus = (k == 0) ? recv_z_minus[i * local_count_y + j] : local_u_curr[IDX_LOCAL(i, j, k - 1, local_count_y, local_count_z)];
    const double u_z_plus  = (k == local_count_z - 1) ? recv_z_plus[i * local_count_y + j] : local_u_curr[IDX_LOCAL(i, j, k + 1, local_count_y, local_count_z)];

    const double u_prev_center = local_u_prev[IDX_LOCAL(i, j, k, local_count_y, local_count_z)];

    return 2*u_center - u_prev_center +
           r_x * (u_x_plus - 2*u_center + u_x_minus) +
           r_y * (u_y_plus - 2*u_center + u_y_minus) +
           r_z * (u_z_plus - 2*u_center + u_z_minus);
}


inline void allocate_aligned(double** ptr, size_t count) {
    // POWER8 推荐 128 字节对齐 (Cache Line Size)
    // POWER8 рекомендует выравнивание 128 байт (размер строки кэша)
    size_t alignment = 128;
    size_t bytes = count * sizeof(double);

    // posix_memalign 返回 0 表示成功
    // posix_memalign возвращает 0 в случае успеха
    if (posix_memalign((void**)ptr, alignment, bytes) != 0) {
        std::cerr << "Memory allocation failed!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

void solve_serial(int count_x, int count_y, int count_z, int count_t,
                  double delta_x, double delta_y, double delta_z, double delta_t,
                  double Lx, double Ly, double Lz,
                  double a_squared, double a, double T)
{
    clock_t begin0, end0;
    #define IDX(i, j, k) ((size_t)(i) * count_y * count_z + (size_t)(j) * count_z + (k))
    size_t total_size = (size_t)count_x * count_y * count_z;
    double *u_prev, *u_curr, *u_next;
    allocate_aligned(&u_prev, total_size);
    allocate_aligned(&u_curr, total_size);
    allocate_aligned(&u_next, total_size);

    // 初始化 u^0 (t=0)
    // Инициализация u^0 (t=0)
    for(int i = 0; i < count_x; i++) {
        for(int j = 0; j < count_y; j++) {
            for(int k = 0; k < count_z; k++) {
                double x = i * delta_x;
                double y = j * delta_y;
                double z = k * delta_z;
                u_prev[IDX(i, j, k)] = phi(x, y, z, Lx, Ly, Lz);
            }
        }
    }

    double r_x = a_squared * delta_t * delta_t / (delta_x * delta_x);
    double r_y = a_squared * delta_t * delta_t / (delta_y * delta_y);
    double r_z = a_squared * delta_t * delta_t / (delta_z * delta_z);

    // 修正 u^1 (t=1) 的初始化
    // Исправленная инициализация u^1 (t=1)
    for(int i = 0; i < count_x; i++) {
        for(int j = 0; j < count_y; j++) {
            for(int k = 0; k < count_z; k++) {
                // *** 修改 1: 只对 Z 轴应用 Dirichlet 边界 (固定为0) ***
                // *** Изменение 1: Применение граничных условий Дирихле только к оси Z (фиксировано 0) ***
                if (k == 0 || k == count_z - 1) {
                    u_curr[IDX(i, j, k)] = 0.0;
                }
                else {
                    // *** 修改 2: X 和 Y 轴应用周期性环绕逻辑 ***
                    // *** Изменение 2: Применение периодической логики к осям X и Y ***

                    // X 轴周期性处理
                    // Периодическая обработка по оси X
                    int i_minus = (i == 0) ? count_x - 1 : i - 1; // 左边越界回到最右
                                                                   // Левая граница возвращается к правой
                    int i_plus  = (i == count_x - 1) ? 0 : i + 1; // 右边越界回到最左
                                                                   // Правая граница возвращается к левой

                    // Y 轴周期性处理
                    // Периодическая обработка по оси Y
                    int j_minus = (j == 0) ? count_y - 1 : j - 1; // 上边越界回到最下
                                                                   // Верхняя граница возвращается к нижней
                    int j_plus  = (j == count_y - 1) ? 0 : j + 1; // 下边越界回到最上
                                                                   // Нижняя граница возвращается к верхней

                    // Z 轴依然是普通邻居 (因为上面已经排除了边界 k)
                    // Ось Z остаётся обычными соседями (так как граница k исключена выше)
                    int k_minus = k - 1;
                    int k_plus  = k + 1;

                    const double u0_center = u_prev[IDX(i, j, k)];
                    const double u0_x_minus = u_prev[IDX(i_minus, j, k)];
                    const double u0_x_plus  = u_prev[IDX(i_plus, j, k)];
                    const double u0_y_minus = u_prev[IDX(i, j_minus, k)];
                    const double u0_y_plus  = u_prev[IDX(i, j_plus, k)];
                    const double u0_z_minus = u_prev[IDX(i, j, k_minus)];
                    const double u0_z_plus  = u_prev[IDX(i, j, k_plus)];

                    double laplacian_u0 = (u0_x_plus - 2*u0_center + u0_x_minus) / (delta_x * delta_x) +
                                          (u0_y_plus - 2*u0_center + u0_y_minus) / (delta_y * delta_y) +
                                          (u0_z_plus - 2*u0_center + u0_z_minus) / (delta_z * delta_z);
                    u_curr[IDX(i, j, k)] = u0_center + 0.5 * a_squared * delta_t * delta_t * laplacian_u0;
                }
            }
        }
    }

    begin0 = clock();
    for(int n = 1; n < count_t; n++) {
        for(int i = 0; i < count_x; i++) {
            for(int j = 0; j < count_y; j++) {
                for(int k = 0; k < count_z; k++) {
                    // *** 修改 3: 主循环同样逻辑，只固定 Z 轴 ***
                    // *** Изменение 3: Та же логика для основного цикла, фиксируется только ось Z ***
                    if (k == 0 || k == count_z - 1) {
                        u_next[IDX(i, j, k)] = 0.0;
                    }
                    else {
                        // X 轴周期性
                        // Периодичность по оси X
                        int i_minus = (i == 0) ? count_x - 1 : i - 1;
                        int i_plus  = (i == count_x - 1) ? 0 : i + 1;

                        // Y 轴周期性
                        // Периодичность по оси Y
                        int j_minus = (j == 0) ? count_y - 1 : j - 1;
                        int j_plus  = (j == count_y - 1) ? 0 : j + 1;

                        // Z 轴普通
                        // Обычная ось Z
                        int k_minus = k - 1;
                        int k_plus  = k + 1;

                        const double u_center = u_curr[IDX(i, j, k)];

                        const double u_x_minus = u_curr[IDX(i_minus, j, k)];
                        const double u_x_plus  = u_curr[IDX(i_plus, j, k)];

                        const double u_y_minus = u_curr[IDX(i, j_minus, k)];
                        const double u_y_plus  = u_curr[IDX(i, j_plus, k)];

                        const double u_z_minus = u_curr[IDX(i, j, k_minus)];
                        const double u_z_plus  = u_curr[IDX(i, j, k_plus)];

                        u_next[IDX(i, j, k)] = 2*u_center - u_prev[IDX(i, j, k)] +
                            r_x * (u_x_plus - 2*u_center + u_x_minus) +
                            r_y * (u_y_plus - 2*u_center + u_y_minus) +
                            r_z * (u_z_plus - 2*u_center + u_z_minus);
                    }
                }
            }
        }
        double *temp = u_prev;
        u_prev = u_curr;
        u_curr = u_next;
        u_next = temp;
    }
    end0 = clock();
    double final_error = calculate_error(u_curr, T, delta_x, delta_y, delta_z,
                                         count_x, count_y, count_z, Lx, Ly, Lz, a);
    cout << "\n=== Results (1 Process) ===" << endl;
    cout << "Calculation time: "  << (double)(end0-begin0)/CLOCKS_PER_SEC << " s" << endl;
    cout << "Final max error: " << scientific << final_error << endl;
    free(u_prev);
    free(u_curr);
    free(u_next);
    #undef IDX
}


/**
 * 并行求解器主函数 (当进程数 > 1 时调用)。
 * * 使用 MPI 笛卡尔拓扑进行 3D 域分解。
 * 1. 并行初始化：每个进程根据全局坐标独立计算其子区域的初始条件 u^0 (无需 Rank 0 分发)。
 * 2. 计算 t=1 (u_curr)，使用非阻塞通信 (Isend/Irecv) 交换 Halo。
 * 3. (性能优化) 使用高效的“盒子循环”(Box Method) 分别计算内部点和边界点，以最小化循环开销。
 * 4. 进入主时间循环 (t=1 到 count_t)，重复 Halo 交换和优化的计算。
 * 5. 分布式误差计算：每个进程计算本地误差，使用 MPI_Reduce 汇总全局最大误差 (无需收集全局网格)。
 * Основная функция параллельного решателя (вызывается, когда np > 1).
 * * Использует декартову топологию MPI для 3D-декомпозиции области.
 * 1. Параллельная инициализация: каждый процесс независимо вычисляет начальное условие u^0 для своей локальной подобласти (без Scatter).
 * 2. Вычисляет t=1 (u_curr), используя неблокирующие обмены (Isend/Irecv) для Halo-областей.
 * 3. (Оптимизация) Использует эффективный "Box" цикл для раздельного вычисления внутренних и граничных точек, минимизируя накладные расходы на циклы.
 * 4. Входит в основной временной цикл (t=1 до count_t), повторяя обмены Halo и оптимизированные вычисления.
 * 5. Распределенное вычисление ошибки: каждый процесс считает локальную ошибку, глобальный максимум собирается через MPI_Reduce (без Gather).
 */
void solve_parallel(int rank, int size, MPI_Comm comm, MPI_Comm cart_comm,
                    int count_x, int count_y, int count_z, int count_t,
                    double delta_x, double delta_y, double delta_z, double delta_t,
                    double Lx, double Ly, double Lz,
                    double a_squared, double a, double T,
                    int local_count_x, int local_count_y, int local_count_z,
                    size_t submatrix_elements,
                    int* outer_coords, int* dims)
{
    // 1. 定义指针
    // 1. Определение указателей
    double *local_u_prev, *local_u_curr, *local_u_next;
    double *tmp_x_minus, *tmp_x_plus, *tmp_y_minus, *tmp_y_plus, *tmp_z_minus, *tmp_z_plus;
    double *recv_x_minus, *recv_x_plus, *recv_y_minus, *recv_y_plus, *recv_z_minus, *recv_z_plus;

    // 2. 统一分配
    // 2. Единое выделение памяти
    allocate_aligned(&local_u_prev, submatrix_elements);
    allocate_aligned(&local_u_curr, submatrix_elements);
    allocate_aligned(&local_u_next, submatrix_elements);

    size_t size_x_face = local_count_y * local_count_z;
    size_t size_y_face = local_count_x * local_count_z;
    size_t size_z_face = local_count_x * local_count_y;

    allocate_aligned(&tmp_x_minus, size_x_face); allocate_aligned(&tmp_x_plus,  size_x_face);
    allocate_aligned(&recv_x_minus, size_x_face); allocate_aligned(&recv_x_plus,  size_x_face);

    allocate_aligned(&tmp_y_minus, size_y_face); allocate_aligned(&tmp_y_plus,  size_y_face);
    allocate_aligned(&recv_y_minus, size_y_face); allocate_aligned(&recv_y_plus,  size_y_face);

    allocate_aligned(&tmp_z_minus, size_z_face); allocate_aligned(&tmp_z_plus,  size_z_face);
    allocate_aligned(&recv_z_minus, size_z_face); allocate_aligned(&recv_z_plus,  size_z_face);

    // 获取邻居 (如果是周期性，boundary 处会返回由于循环连接的 rank，而不是 MPI_PROC_NULL)
    // Получение соседей (если периодичность включена, на границе будет возвращён ранг циклического соединения, а не MPI_PROC_NULL)
    int x_minus_neighbor, x_plus_neighbor, y_minus_neighbor, y_plus_neighbor, z_minus_neighbor, z_plus_neighbor;
    MPI_Cart_shift(cart_comm, 0, 1, &x_minus_neighbor, &x_plus_neighbor);
    MPI_Cart_shift(cart_comm, 1, 1, &y_minus_neighbor, &y_plus_neighbor);
    MPI_Cart_shift(cart_comm, 2, 1, &z_minus_neighbor, &z_plus_neighbor);

    MPI_Request requests[12];

    const int TAG_X_MINUS_DATA = 10; const int TAG_X_PLUS_DATA = 11;
    const int TAG_Y_MINUS_DATA = 20; const int TAG_Y_PLUS_DATA = 21;
    const int TAG_Z_MINUS_DATA = 30; const int TAG_Z_PLUS_DATA = 31;

    double r_x = a_squared * delta_t * delta_t / (delta_x * delta_x);
    double r_y = a_squared * delta_t * delta_t / (delta_y * delta_y);
    double r_z = a_squared * delta_t * delta_t / (delta_z * delta_z);

    double start_time, end_time, local_time_diff, global_time_diff;

    MPI_Barrier(comm);

    // 计算当前进程在全局网格中的起始偏移量 (Global Offsets)
    // Вычисление начальных смещений текущего процесса в глобальной сетке (Global Offsets)
    int start_i_global = 0;
    for (int c = 0; c < outer_coords[0]; ++c) {
        int size_c = count_x / dims[0];
        if (c < (count_x % dims[0])) size_c++;
        start_i_global += size_c;
    }
    int start_j_global = 0;
    for (int c = 0; c < outer_coords[1]; ++c) {
        int size_c = count_y / dims[1];
        if (c < (count_y % dims[1])) size_c++;
        start_j_global += size_c;
    }
    int start_k_global = 0;
    for (int c = 0; c < outer_coords[2]; ++c) {
        int size_c = count_z / dims[2];
        if (c < (count_z % dims[2])) size_c++;
        start_k_global += size_c;
    }

    // 每个进程独立计算自己这块区域的初始值
    // Каждый процесс независимо вычисляет начальные значения для своей области
    for(int i = 0; i < local_count_x; i++) {
        for(int j = 0; j < local_count_y; j++) {
            for(int k = 0; k < local_count_z; k++) {
                double x = (start_i_global + i) * delta_x;
                double y = (start_j_global + j) * delta_y;
                double z = (start_k_global + k) * delta_z;
                local_u_prev[IDX_LOCAL(i, j, k, local_count_y, local_count_z)] = phi(x, y, z, Lx, Ly, Lz);
            }
        }
    }

    // === 3. 并行计算 u^1 (u_curr) ===
    // === 3. Параллельное вычисление u^1 (u_curr) ===

    // 3a. 准备发送缓冲区
    // 3a. Подготовка буферов отправки
    memcpy(tmp_x_minus, local_u_prev, size_x_face * sizeof(double));
    memcpy(tmp_x_plus, local_u_prev + (local_count_x - 1) * size_x_face, size_x_face * sizeof(double));
    Get_Y_Face_Minus(local_u_prev, tmp_y_minus, local_count_x, local_count_y, local_count_z);
    Get_Y_Face_Plus(local_u_prev, tmp_y_plus, local_count_x, local_count_y, local_count_z);
    if (z_minus_neighbor != MPI_PROC_NULL) Get_Z_Face_Minus(local_u_prev, tmp_z_minus, local_count_x, local_count_y, local_count_z);
    if (z_plus_neighbor != MPI_PROC_NULL)  Get_Z_Face_Plus(local_u_prev, tmp_z_plus, local_count_x, local_count_y, local_count_z);

    // 3b. 非阻塞通信 (使用邻居 rank 判断是否需要通信)
    // 3b. Неблокирующая связь (используя ранг соседа для определения необходимости обмена)
    for(int i=0; i<12; ++i) requests[i] = MPI_REQUEST_NULL;


    MPI_Isend(tmp_x_minus, size_x_face, MPI_DOUBLE, x_minus_neighbor, TAG_X_MINUS_DATA, cart_comm, &requests[0]);
    MPI_Irecv(recv_x_minus, size_x_face, MPI_DOUBLE, x_minus_neighbor, TAG_X_PLUS_DATA, cart_comm, &requests[1]);
    MPI_Isend(tmp_x_plus, size_x_face, MPI_DOUBLE, x_plus_neighbor, TAG_X_PLUS_DATA, cart_comm, &requests[2]);
    MPI_Irecv(recv_x_plus, size_x_face, MPI_DOUBLE, x_plus_neighbor, TAG_X_MINUS_DATA, cart_comm, &requests[3]);
    MPI_Isend(tmp_y_minus, size_y_face, MPI_DOUBLE, y_minus_neighbor, TAG_Y_MINUS_DATA, cart_comm, &requests[4]);
    MPI_Irecv(recv_y_minus, size_y_face, MPI_DOUBLE, y_minus_neighbor, TAG_Y_PLUS_DATA, cart_comm, &requests[5]);
    MPI_Isend(tmp_y_plus, size_y_face, MPI_DOUBLE, y_plus_neighbor, TAG_Y_PLUS_DATA, cart_comm, &requests[6]);
    MPI_Irecv(recv_y_plus, size_y_face, MPI_DOUBLE, y_plus_neighbor, TAG_Y_MINUS_DATA, cart_comm, &requests[7]);
    if (z_minus_neighbor != MPI_PROC_NULL) {
        MPI_Isend(tmp_z_minus, size_z_face, MPI_DOUBLE, z_minus_neighbor, TAG_Z_MINUS_DATA, cart_comm, &requests[8]);
        MPI_Irecv(recv_z_minus, size_z_face, MPI_DOUBLE, z_minus_neighbor, TAG_Z_PLUS_DATA, cart_comm, &requests[9]);
    }
    if (z_plus_neighbor != MPI_PROC_NULL) {
        MPI_Isend(tmp_z_plus, size_z_face, MPI_DOUBLE, z_plus_neighbor, TAG_Z_PLUS_DATA, cart_comm, &requests[10]);
        MPI_Irecv(recv_z_plus, size_z_face, MPI_DOUBLE, z_plus_neighbor, TAG_Z_MINUS_DATA, cart_comm, &requests[11]);
    }

    // 3c. 计算内部点 (1 到 N-2)
    // 3c. Вычисление внутренних точек (от 1 до N-2)
    for (int i = 1; i < local_count_x - 1; ++i) {
        for (int j = 1; j < local_count_y - 1; ++j) {
            for (int k = 1; k < local_count_z - 1; ++k) {
                local_u_curr[IDX_LOCAL(i,j,k, local_count_y, local_count_z)] = CalculateStencil_t1(
                    i, j, k, local_u_prev,
                    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                    local_count_x, local_count_y, local_count_z,
                    delta_x, delta_y, delta_z, a_squared, delta_t
                );
            }
        }
    }

    // 3d. 等待通信
    // 3d. Ожидание завершения обмена
    MPI_Waitall(12, requests, MPI_STATUSES_IGNORE);

    // 3e. 计算边界层 ("Box" Method) - 仅当有邻居时才计算 (周期性或内部边界)
    // 3e. Вычисление граничного слоя (метод "Box") - вычисляется только при наличии соседей (периодичность или внутренняя граница)
    // X Faces
    int i = 0;
    for (int j = 0; j < local_count_y; ++j) for (int k = 0; k < local_count_z; ++k)
        local_u_curr[IDX_LOCAL(i,j,k, local_count_y, local_count_z)] = CalculateStencil_t1(
            i, j, k, local_u_prev, recv_x_minus, recv_x_plus, recv_y_minus, recv_y_plus, recv_z_minus, recv_z_plus,
            local_count_x, local_count_y, local_count_z, delta_x, delta_y, delta_z, a_squared, delta_t);

    i = local_count_x - 1;
    for (int j = 0; j < local_count_y; ++j) for (int k = 0; k < local_count_z; ++k)
        local_u_curr[IDX_LOCAL(i,j,k, local_count_y, local_count_z)] = CalculateStencil_t1(
            i, j, k, local_u_prev, recv_x_minus, recv_x_plus, recv_y_minus, recv_y_plus, recv_z_minus, recv_z_plus,
            local_count_x, local_count_y, local_count_z, delta_x, delta_y, delta_z, a_squared, delta_t);

    // Y Faces
    int j = 0;
    for (int i = 1; i < local_count_x - 1; ++i) for (int k = 0; k < local_count_z; ++k)
        local_u_curr[IDX_LOCAL(i,j,k, local_count_y, local_count_z)] = CalculateStencil_t1(
            i, j, k, local_u_prev, recv_x_minus, recv_x_plus, recv_y_minus, recv_y_plus, recv_z_minus, recv_z_plus,
            local_count_x, local_count_y, local_count_z, delta_x, delta_y, delta_z, a_squared, delta_t);

    j = local_count_y - 1;
    for (int i = 1; i < local_count_x - 1; ++i) for (int k = 0; k < local_count_z; ++k)
        local_u_curr[IDX_LOCAL(i,j,k, local_count_y, local_count_z)] = CalculateStencil_t1(
            i, j, k, local_u_prev, recv_x_minus, recv_x_plus, recv_y_minus, recv_y_plus, recv_z_minus, recv_z_plus,
            local_count_x, local_count_y, local_count_z, delta_x, delta_y, delta_z, a_squared, delta_t);

    // Z Faces
    if (z_minus_neighbor != MPI_PROC_NULL) {
        int k = 0;
        for (int i = 1; i < local_count_x - 1; ++i) for (int j = 1; j < local_count_y - 1; ++j)
            local_u_curr[IDX_LOCAL(i,j,k, local_count_y, local_count_z)] = CalculateStencil_t1(
                i, j, k, local_u_prev, recv_x_minus, recv_x_plus, recv_y_minus, recv_y_plus, recv_z_minus, recv_z_plus,
                local_count_x, local_count_y, local_count_z, delta_x, delta_y, delta_z, a_squared, delta_t);
    }
    if (z_plus_neighbor != MPI_PROC_NULL) {
        int k = local_count_z - 1;
        for (int i = 1; i < local_count_x - 1; ++i) for (int j = 1; j < local_count_y - 1; ++j)
            local_u_curr[IDX_LOCAL(i,j,k, local_count_y, local_count_z)] = CalculateStencil_t1(
                i, j, k, local_u_prev, recv_x_minus, recv_x_plus, recv_y_minus, recv_y_plus, recv_z_minus, recv_z_plus,
                local_count_x, local_count_y, local_count_z, delta_x, delta_y, delta_z, a_squared, delta_t);
    }

    // 3f. 应用 Dirichlet 边界条件 (仅当没有邻居时强制设为0，适用于物理墙壁)
    // 3f. Применение граничных условий Дирихле (принудительно 0 только при отсутствии соседей, для физических стенок)
    if (z_minus_neighbor == MPI_PROC_NULL) {
        for(int i=0; i<local_count_x; ++i) for(int j=0; j<local_count_y; ++j) local_u_curr[IDX_LOCAL(i,j,0, local_count_y, local_count_z)] = 0.0;
    }
    if (z_plus_neighbor == MPI_PROC_NULL) {
        for(int i=0; i<local_count_x; ++i) for(int j=0; j<local_count_y; ++j) local_u_curr[IDX_LOCAL(i,j,local_count_z-1, local_count_y, local_count_z)] = 0.0;
    }


    // === 4. MPI 主时间循环 ===
    // === 4. Основной временной цикл MPI ===
    MPI_Barrier(comm);
    start_time = MPI_Wtime();

    for (int n = 1; n < count_t; n++) {

        // 4a. 拷贝数据到 Buffer (仅当有邻居时)
        // 4a. Копирование данных в буфер (только при наличии соседей)
        memcpy(tmp_x_minus, local_u_curr, size_x_face * sizeof(double));
        memcpy(tmp_x_plus, local_u_curr + (local_count_x - 1) * size_x_face, size_x_face * sizeof(double));
        Get_Y_Face_Minus(local_u_curr, tmp_y_minus, local_count_x, local_count_y, local_count_z);
        Get_Y_Face_Plus(local_u_curr, tmp_y_plus, local_count_x, local_count_y, local_count_z);
        if (z_minus_neighbor != MPI_PROC_NULL) Get_Z_Face_Minus(local_u_curr, tmp_z_minus, local_count_x, local_count_y, local_count_z);
        if (z_plus_neighbor != MPI_PROC_NULL)  Get_Z_Face_Plus(local_u_curr, tmp_z_plus, local_count_x, local_count_y, local_count_z);

        // 4b. 开始通信
        // 4b. Начало обмена данными
        for(int i=0; i<12; ++i) requests[i] = MPI_REQUEST_NULL;

        MPI_Isend(tmp_x_minus, size_x_face, MPI_DOUBLE, x_minus_neighbor, TAG_X_MINUS_DATA, cart_comm, &requests[0]);
        MPI_Irecv(recv_x_minus, size_x_face, MPI_DOUBLE, x_minus_neighbor, TAG_X_PLUS_DATA, cart_comm, &requests[1]);
        MPI_Isend(tmp_x_plus, size_x_face, MPI_DOUBLE, x_plus_neighbor, TAG_X_PLUS_DATA, cart_comm, &requests[2]);
        MPI_Irecv(recv_x_plus, size_x_face, MPI_DOUBLE, x_plus_neighbor, TAG_X_MINUS_DATA, cart_comm, &requests[3]);
        MPI_Isend(tmp_y_minus, size_y_face, MPI_DOUBLE, y_minus_neighbor, TAG_Y_MINUS_DATA, cart_comm, &requests[4]);
        MPI_Irecv(recv_y_minus, size_y_face, MPI_DOUBLE, y_minus_neighbor, TAG_Y_PLUS_DATA, cart_comm, &requests[5]);
        MPI_Isend(tmp_y_plus, size_y_face, MPI_DOUBLE, y_plus_neighbor, TAG_Y_PLUS_DATA, cart_comm, &requests[6]);
        MPI_Irecv(recv_y_plus, size_y_face, MPI_DOUBLE, y_plus_neighbor, TAG_Y_MINUS_DATA, cart_comm, &requests[7]);

        if (z_minus_neighbor != MPI_PROC_NULL) {
            MPI_Isend(tmp_z_minus, size_z_face, MPI_DOUBLE, z_minus_neighbor, TAG_Z_MINUS_DATA, cart_comm, &requests[8]);
            MPI_Irecv(recv_z_minus, size_z_face, MPI_DOUBLE, z_minus_neighbor, TAG_Z_PLUS_DATA, cart_comm, &requests[9]);
        }
        if (z_plus_neighbor != MPI_PROC_NULL) {
            MPI_Isend(tmp_z_plus, size_z_face, MPI_DOUBLE, z_plus_neighbor, TAG_Z_PLUS_DATA, cart_comm, &requests[10]);
            MPI_Irecv(recv_z_plus, size_z_face, MPI_DOUBLE, z_plus_neighbor, TAG_Z_MINUS_DATA, cart_comm, &requests[11]);
        }

        // 4c. 计算内部点
        // 4c. Вычисление внутренних точек
        for (int i = 1; i < local_count_x - 1; ++i) {
            for (int j = 1; j < local_count_y - 1; ++j) {
                for (int k = 1; k < local_count_z - 1; ++k) {
                    local_u_next[IDX_LOCAL(i,j,k, local_count_y, local_count_z)] = CalculateStencil_t_n(
                        i, j, k, local_u_prev, local_u_curr,
                        nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                        local_count_x, local_count_y, local_count_z, r_x, r_y, r_z
                    );
                }
            }
        }

        // 4d. 等待
        // 4d. Ожидание
        MPI_Waitall(12, requests, MPI_STATUSES_IGNORE);

        // 4e. 计算边界层 (Box Method) - 使用邻居判断
        // 4e. Вычисление граничного слоя (метод "Box") - используя проверку соседей
        int i = 0;
        for (int j = 0; j < local_count_y; ++j) for (int k = 0; k < local_count_z; ++k)
            local_u_next[IDX_LOCAL(i,j,k, local_count_y, local_count_z)] = CalculateStencil_t_n(
                i, j, k, local_u_prev, local_u_curr, recv_x_minus, recv_x_plus, recv_y_minus, recv_y_plus, recv_z_minus, recv_z_plus,
                local_count_x, local_count_y, local_count_z, r_x, r_y, r_z);

        i = local_count_x - 1;
        for (int j = 0; j < local_count_y; ++j) for (int k = 0; k < local_count_z; ++k)
            local_u_next[IDX_LOCAL(i,j,k, local_count_y, local_count_z)] = CalculateStencil_t_n(
                i, j, k, local_u_prev, local_u_curr, recv_x_minus, recv_x_plus, recv_y_minus, recv_y_plus, recv_z_minus, recv_z_plus,
                local_count_x, local_count_y, local_count_z, r_x, r_y, r_z);

        int j = 0;
        for (int i = 1; i < local_count_x - 1; ++i) for (int k = 0; k < local_count_z; ++k)
            local_u_next[IDX_LOCAL(i,j,k, local_count_y, local_count_z)] = CalculateStencil_t_n(
                i, j, k, local_u_prev, local_u_curr, recv_x_minus, recv_x_plus, recv_y_minus, recv_y_plus, recv_z_minus, recv_z_plus,
                local_count_x, local_count_y, local_count_z, r_x, r_y, r_z);

        j = local_count_y - 1;
        for (int i = 1; i < local_count_x - 1; ++i) for (int k = 0; k < local_count_z; ++k)
            local_u_next[IDX_LOCAL(i,j,k, local_count_y, local_count_z)] = CalculateStencil_t_n(
                i, j, k, local_u_prev, local_u_curr, recv_x_minus, recv_x_plus, recv_y_minus, recv_y_plus, recv_z_minus, recv_z_plus,
                local_count_x, local_count_y, local_count_z, r_x, r_y, r_z);

        if (z_minus_neighbor != MPI_PROC_NULL) {
            int k = 0;
            for (int i = 1; i < local_count_x - 1; ++i) for (int j = 1; j < local_count_y - 1; ++j)
                local_u_next[IDX_LOCAL(i,j,k, local_count_y, local_count_z)] = CalculateStencil_t_n(
                    i, j, k, local_u_prev, local_u_curr, recv_x_minus, recv_x_plus, recv_y_minus, recv_y_plus, recv_z_minus, recv_z_plus,
                    local_count_x, local_count_y, local_count_z, r_x, r_y, r_z);
        }
        if (z_plus_neighbor != MPI_PROC_NULL) {
            int k = local_count_z - 1;
            for (int i = 1; i < local_count_x - 1; ++i) for (int j = 1; j < local_count_y - 1; ++j)
                local_u_next[IDX_LOCAL(i,j,k, local_count_y, local_count_z)] = CalculateStencil_t_n(
                    i, j, k, local_u_prev, local_u_curr, recv_x_minus, recv_x_plus, recv_y_minus, recv_y_plus, recv_z_minus, recv_z_plus,
                    local_count_x, local_count_y, local_count_z, r_x, r_y, r_z);
        }

        // 4f. 应用 Dirichlet 边界条件 (如果没有邻居)
        // 4f. Применение граничных условий Дирихле (если нет соседей)
        if (z_minus_neighbor == MPI_PROC_NULL) {
            for(int i=0; i<local_count_x; ++i) for(int j=0; j<local_count_y; ++j) local_u_next[IDX_LOCAL(i,j,0, local_count_y, local_count_z)] = 0.0;
        }
        if (z_plus_neighbor == MPI_PROC_NULL) {
            for(int i=0; i<local_count_x; ++i) for(int j=0; j<local_count_y; ++j) local_u_next[IDX_LOCAL(i,j,local_count_z-1, local_count_y, local_count_z)] = 0.0;
        }

        // 4g. 交换指针
        // 4g. Обмен указателями
        double *temp = local_u_prev;
        local_u_prev = local_u_curr;
        local_u_curr = local_u_next;
        local_u_next = temp;
    }

    end_time = MPI_Wtime();
    local_time_diff = end_time - start_time;
    MPI_Reduce(&local_time_diff, &global_time_diff, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    double local_max_error = calculate_local_error(
        local_u_curr, T, delta_x, delta_y, delta_z,
        local_count_x, local_count_y, local_count_z,
        start_i_global, start_j_global, start_k_global,
        Lx, Ly, Lz
    );

    double global_max_error = 0.0;
    MPI_Reduce(&local_max_error, &global_max_error, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

    if (rank == 0) {
        cout << "\n=== Results (" << size << " Processes) ===" << endl;
        cout << "Calculation time: "  << global_time_diff << " s" << endl;
        cout << "Final max error: " << scientific << global_max_error << endl;
    }

    free(local_u_prev); free(local_u_curr); free(local_u_next);
    free(tmp_x_minus); free(tmp_x_plus);
    free(tmp_y_minus); free(tmp_y_plus);
    free(tmp_z_minus); free(tmp_z_plus);
    free(recv_x_minus); free(recv_x_plus);
    free(recv_y_minus); free(recv_y_plus);
    free(recv_z_minus); free(recv_z_plus);
}



int main(int argc, char **argv) {

    int grid_size = 256;
    int time_steps = 20;
    if(argc >= 2) {
        grid_size = atoi(argv[1]);
    }
    if(argc >= 3) {
        time_steps = atoi(argv[2]);
    }

    int count_x = grid_size;
    int count_y = grid_size;
    int count_z = grid_size;
    int count_t = time_steps;

    double Lx = 1.0, Ly = 1.0, Lz = 1.0;
    double delta_x = Lx / (count_x);
    double delta_y = Ly / (count_y);
    double delta_z = Lz / (count_z - 1);
    double T = 0.01;
    double delta_t = T / count_t;
    double a_squared = calculate_a_squared(Lx, Ly, Lz);
    double a = sqrt(a_squared);

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm cart_comm;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(!check_stability(a, delta_x, delta_y, delta_z, delta_t, rank)) {
        if (rank == 0) {
            cout << "Warning: Stability condition may not be satisfied!" << endl;
        }
    }

    int dims[3] = {0, 0, 0};
    int periods[3] = {1, 1, 0};
    int outer_coords[3];
    MPI_Dims_create(size, 3, dims);
    MPI_Cart_create(comm, 3, dims, periods, 0, &cart_comm);
    MPI_Cart_coords(cart_comm, rank, 3, outer_coords);

    int local_count_x, local_count_y, local_count_z;

    int extra_count_x = count_x % dims[0];
    int extra_count_y = count_y % dims[1];
    int extra_count_z = count_z % dims[2];
    local_count_x = count_x / dims[0];
    local_count_y = count_y / dims[1];
    local_count_z = count_z / dims[2];
    if(outer_coords[0] < extra_count_x) local_count_x++;
    if(outer_coords[1] < extra_count_y) local_count_y++;
    if(outer_coords[2] < extra_count_z) local_count_z++;
    size_t submatrix_elements = (size_t)local_count_x * local_count_y * local_count_z;

    if (rank == 0) {
        cout << "Wave equation solver - Variant 7" << endl;
        cout << "Grid: " << count_x << "x" << count_y << "x" << count_z << endl;
        cout << "Time steps: " << count_t << endl;
        cout << "a² = " << a_squared << ", a = " << a << endl;
        cout << "MPI Processes: " << size << " (Grid: " << dims[0] << "x" << dims[1] << "x" << dims[2] << ")" << endl;
    }

    if (size == 1) {
        solve_serial(count_x, count_y, count_z, count_t,
                     delta_x, delta_y, delta_z, delta_t,
                     Lx, Ly, Lz, a_squared, a, T);
    }
    else {
        solve_parallel(rank, size, comm, cart_comm,
                       count_x, count_y, count_z, count_t,
                       delta_x, delta_y, delta_z, delta_t,
                       Lx, Ly, Lz, a_squared, a, T,
                       local_count_x, local_count_y, local_count_z,
                       submatrix_elements,
                       outer_coords, dims);
    }

    MPI_Finalize();
    return 0;
}