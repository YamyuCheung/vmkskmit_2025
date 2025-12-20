#include <iostream>
#include <cmath>
#include <ctime>
#include <omp.h>
#include <iomanip>
using namespace std;

const double PI = 3.14159265358979323846;

// 解析解函数
// Функция аналитического решения
double u_analytical(double t, double x, double y, double z,
                    double Lx, double Ly, double Lz) {
    double at = PI;  // 变体7中 at = π
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

// 计算误差
// Вычисление ошибки
// --- 更改点: 数组类型从 double*** u 改为 double* u ---
double calculate_error(double *u, double t,
                       double delta_x, double delta_y, double delta_z,
                       int count_x, int count_y, int count_z,
                       double Lx, double Ly, double Lz, double a) {
    double max_error = 0.0;

    // --- 更改点: 定义 1D 索引宏 ---
    // (i, j, k) -> i * Ny*Nz + j * Nz + k
    #define IDX(i, j, k) ((size_t)(i) * count_y * count_z + (size_t)(j) * count_z + (k))

    #pragma omp parallel for reduction(max:max_error)
    for(int i = 0; i < count_x; i++) {
        for(int j = 0; j < count_y; j++) {
            for(int k = 0; k < count_z; k++) {
                double x = i * delta_x;
                double y = j * delta_y;
                double z = k * delta_z;
                double exact = u_analytical(t, x, y, z, Lx, Ly, Lz);
                // --- 更改点: 使用 1D 索引访问 ---
                double error = fabs(u[IDX(i, j, k)] - exact);
                if(error > max_error) max_error = error;
            }
        }
    }
    // --- 更改点: 释放宏定义 ---
    #undef IDX
    return max_error;
}

// 检查稳定性条件（CFL条件）
// Проверка условия устойчивости (условие CFL)
bool check_stability(double a, double delta_x, double delta_y, double delta_z, double delta_t) {

    // 计算 3D CFL 条件
    // Вычисление 3D условия CFL
    double cfl_3D = a * delta_t * sqrt( 1.0/(delta_x*delta_x) +
                                       1.0/(delta_y*delta_y) +
                                       1.0/(delta_z*delta_z) );

    cout << "3D CFL number = " << cfl_3D << " (should be <= 1 for stability)" << endl;

    return cfl_3D <= 1.0;
}

// 辅助函数：获取 X 方向的周期性索引
// Вспомогательная функция: получение периодического индекса по направлению X
inline void get_periodic_indices_x(int i, int count_x, int& i_minus, int& i_plus) {
    // 0 的左边是 N-1
    i_minus = (i == 0) ? count_x - 1 : i - 1;
    // N-1 的右边是 0
    i_plus  = (i == count_x - 1) ? 0 : i + 1;
}

// 辅助函数：获取 Y 方向的周期性索引
// Вспомогательная функция: получение периодического индекса по направлению Y
inline void get_periodic_indices_y(int j, int count_y, int& j_minus, int& j_plus) {
    j_minus = (j == 0) ? count_y - 1 : j - 1;
    j_plus  = (j == count_y - 1) ? 0 : j + 1;
}


int main(int argc, const char **argv) {

    // 解析命令行参数
    // Разбор параметров командной строки
    int grid_size = 128;
    int time_steps = 20;
    int num_threads = 4;
    if(argc >= 2) {
        grid_size = atoi(argv[1]);
    }
    if(argc >= 3) {
        time_steps = atoi(argv[2]);
    }
    if(argc >= 4) {
        num_threads = atoi(argv[3]);
    }

    // 网格参数
    // Параметры сетки
    int count_x = grid_size;
    int count_y = grid_size;
    int count_z = grid_size;
    int count_t = time_steps;

    // 域的大小
    // Размер области
    double Lx = 1.0, Ly = 1.0, Lz = 1.0;

    // 计算网格步长
    // Вычисление шагов сетки
    double delta_x = Lx / (count_x);
    double delta_y = Ly / (count_y);
    double delta_z = Lz / (count_z - 1);

    // 时间参数
    // Параметры времени
    double T = 0.01;  // 总时间 // Общее время
    double delta_t = T / count_t;

    // 计算a²
    // Вычисление a²
    double a_squared = calculate_a_squared(Lx, Ly, Lz);
    double a = sqrt(a_squared);

    cout << "Wave equation solver - Variant 7" << endl;
    cout << "Grid: " << count_x << "x" << count_y << "x" << count_z << endl;
    cout << "Time steps: " << count_t << endl;
    cout << "a² = " << a_squared << ", a = " << a << endl;

    // 检查稳定性
    // Проверка устойчивости
    if(!check_stability(a, delta_x, delta_y, delta_z, delta_t)) {
        cout << "Warning: Stability condition may not be satisfied!" << endl;
    }

    omp_set_num_threads(num_threads);

    // (i, j, k) -> i * Ny*Nz + j * Nz + k
    // 使用 size_t 来防止大网格时索引溢出
    #define IDX(i, j, k) ((size_t)(i) * count_y * count_z + (size_t)(j) * count_z + (k))

    // 计算总大小
    size_t total_size = (size_t)count_x * count_y * count_z;

    // 声明并分配一维数组
    double *u_prev = new double[total_size];
    double *u_curr = new double[total_size];
    double *u_next = new double[total_size];

    // 初始化 u^0 (t=0)
    // Инициализация u^0 (t=0)
    #pragma omp parallel for
    for(int i = 0; i < count_x; i++) {
        for(int j = 0; j < count_y; j++) {
            for(int k = 0; k < count_z; k++) {
                double x = i * delta_x;
                double y = j * delta_y;
                double z = k * delta_z;
                // --- 更改点: 使用 1D 索引 ---
                u_prev[IDX(i, j, k)] = phi(x, y, z, Lx, Ly, Lz);
            }
        }
    }

    // 计算 u^1 使用 Taylor 展开
    double r_x = a_squared * delta_t * delta_t / (delta_x * delta_x);
    double r_y = a_squared * delta_t * delta_t / (delta_y * delta_y);
    double r_z = a_squared * delta_t * delta_t / (delta_z * delta_z);

    // 修正 u^1 (t=1) 的初始化
    // Инициализация u^1 (t=1)
    #pragma omp parallel for collapse(2) schedule(static)
    for(int i = 0; i < count_x; i++) {
        for(int j = 0; j < count_y; j++) {
            for(int k = 0; k < count_z; k++) {

                // 應用 Dirichlet (z) 邊界 (условие первого рода)
                // Применение граничных условий первого рода (z)
                if (k == 0 || k == count_z - 1) {
                    // --- 更改点: 使用 1D 索引 ---
                    u_curr[IDX(i, j, k)] = 0.0;
                }
                else {
                    // 內部節點 (z方向)，應用 Periodic (x, y) 邊界

                    // --- X 週期邊界 ---
                    int i1, i2;
                    get_periodic_indices_x(i, count_x, i1, i2);

                    // --- Y 週期邊界 ---
                    int j1, j2;
                    get_periodic_indices_y(j, count_y, j1, j2);

                    // --- Z 內部節點 (非邊界) ---
                    int k1 = k - 1;
                    int k2 = k + 1;

                    // 獲取鄰居值 (u_prev 即 u^0)
                    // --- 更改点: 全部使用 1D 索引 ---
                    const double u0_center = u_prev[IDX(i, j, k)];   // (i, j, k) at t=0
                    const double u0_x_minus = u_prev[IDX(i1, j, k)]; // (i-1, j, k) at t=0
                    const double u0_x_plus = u_prev[IDX(i2, j, k)];  // (i+1, j, k) at t=0
                    const double u0_y_minus = u_prev[IDX(i, j1, k)]; // (i, j-1, k) at t=0
                    const double u0_y_plus = u_prev[IDX(i, j2, k)];  // (i, j+1, k) at t=0
                    const double u0_z_minus = u_prev[IDX(i, j, k1)]; // (i, j, k-1) at t=0
                    const double u0_z_plus = u_prev[IDX(i, j, k2)];  // (i, j, k+1) at t=0

                    // 計算 u^0 的拉普拉斯算子
                    double laplacian_u0 = (u0_x_plus - 2*u0_center + u0_x_minus) / (delta_x * delta_x) +
                                          (u0_y_plus - 2*u0_center + u0_y_minus) / (delta_y * delta_y) +
                                          (u0_z_plus - 2*u0_center + u0_z_minus) / (delta_z * delta_z);

                    // 根據 Taylor 展開計算 u^1
                    // --- 更改点: 使用 1D 索引 ---
                    u_curr[IDX(i, j, k)] = u0_center + 0.5 * a_squared * delta_t * delta_t * laplacian_u0;
                }
            }
        }
    }


    // 开始计时
    // Начало измерения времени
    double begin = omp_get_wtime();

    // 时间推进（从 n=1 开始）
    // Шаг во времени (начиная с n=1)
    for(int n = 1; n < count_t; n++) {

#pragma omp parallel for collapse(2) schedule(static)
        for(int i = 0; i < count_x; i++) {       // 注意：X 从 0 到 N-1 全算
            for(int j = 0; j < count_y; j++) {   // 注意：Y 从 0 到 N-1 全算
                for(int k = 0; k < count_z; k++) {

                    // 1. Z 方向：Dirichlet 边界 (直接置 0)
                    if (k == 0 || k == count_z - 1) {
                        u_next[IDX(i, j, k)] = 0.0;
                        continue; // 跳过后续计算
                    }

                    // 2. X 方向：周期性索引处理
                    int i_minus, i_plus;
                    // 你的辅助函数逻辑是正确的，这里展开写更直观：
                    i_minus = (i == 0) ? count_x - 1 : i - 1;
                    i_plus  = (i == count_x - 1) ? 0 : i + 1;

                    // 3. Y 方向：周期性索引处理
                    int j_minus, j_plus;
                    j_minus = (j == 0) ? count_y - 1 : j - 1;
                    j_plus  = (j == count_y - 1) ? 0 : j + 1;

                    // 4. Z 方向：普通内部点
                    int k_minus = k - 1;
                    int k_plus  = k + 1;

                    // --- 获取值 ---
                    double u_center = u_curr[IDX(i, j, k)];

                    // 使用上面计算出的周期性索引
                    double u_x_minus = u_curr[IDX(i_minus, j, k)];
                    double u_x_plus  = u_curr[IDX(i_plus, j, k)];

                    double u_y_minus = u_curr[IDX(i, j_minus, k)];
                    double u_y_plus  = u_curr[IDX(i, j_plus, k)];

                    double u_z_minus = u_curr[IDX(i, j, k_minus)];
                    double u_z_plus  = u_curr[IDX(i, j, k_plus)];

                    // --- 计算更新 (公式不变) ---
                    u_next[IDX(i, j, k)] = 2*u_center - u_prev[IDX(i, j, k)] +
                        r_x * (u_x_plus - 2*u_center + u_x_minus) +
                        r_y * (u_y_plus - 2*u_center + u_y_minus) +
                        r_z * (u_z_plus - 2*u_center + u_z_minus);
                }
            }
        }


        // 更新时间层
        // (你注释掉的循环更新方式也可以用，但指针交换更高效)
        // --- 更改点: 临时指针类型改为 double* ---
        double *temp = u_prev;
        u_prev = u_curr;
        u_curr = u_next;
        u_next = temp;

    }

    // 结束计时
    // Конец измерения времени
    double end = omp_get_wtime();

    // 计算最终误差
    // (调用不变, 因为 u_curr 已经是 double* 类型)
    double final_error = calculate_error(u_curr, T, delta_x, delta_y, delta_z,
                                         count_x, count_y, count_z, Lx, Ly, Lz, a);

    cout << "\n=== Results ===" << endl;
    cout << "Number of threads: " << omp_get_max_threads() << endl;
    cout << "Calculation time: " << fixed  << end - begin << " s" << endl;
    cout << "Final max error: " << scientific << final_error << endl;

    // 释放 1D 数组 ---
    delete[] u_prev;
    delete[] u_curr;
    delete[] u_next;

    // 释放宏定义 ---
    #undef IDX

    return 0;
}