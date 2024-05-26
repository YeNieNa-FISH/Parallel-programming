#include<iostream>
#include<algorithm>
#include<random>
#include<cstdlib>
#include<vector>
#include<xmmintrin.h>
#include<chrono>
#include<ctime>
#include<immintrin.h>

using namespace std;
using namespace chrono;

static const int n_ = 3000;

static const int row_ = n_;
static const int column_ = n_;

float m[row_][column_];
float m_copy[row_][column_];

void reset(float(*m)[column_], float(*m_copy)[column_]) {
    for (int row = 0; row < row_; row++) {
        for (int col = 0; col < column_; col++) {
            m_copy[row][col] = m[row][col];
        }
    }
}
/*
* 串行算法
*/
void simple_serial_elimination(float(*m)[column_]) {
    auto t1 = system_clock::now();//计时开始
    for (int i = 0; i < column_ - 1; ++i) {
        for (int j = i + 1; j < row_; ++j) {
            float temp = m[j][i] / m[i][i];
            for (int k = i + 1; k < column_; ++k) {
                m[j][k] -= m[i][k] * temp;
            }
            m[j][i] = 0.00f;
        }
    }
    auto t2 = system_clock::now();//计时结束
    cout << "串行用时\t" << duration_cast<nanoseconds>(t2 - t1).count() << '\n';
}

void align_simple_serial_elimination(float(*m_)[column_]) {
    float** m = static_cast<float**>(_aligned_malloc(row_ * sizeof(float*), 16));
    for (size_t i = 0; i < row_; ++i) {
        m[i] = static_cast<float*>(_aligned_malloc(column_ * sizeof(float), 16));
    }
    for (int i = 0; i < row_; i++)
    {
        for (int j = 0; j < column_; j++)
        {
            m[i][j] = m_[i][j];
        }
    }
    auto t1 = system_clock::now();//计时开始
    for (int i = 0; i < column_ - 1; ++i) {
        for (int j = i + 1; j < row_; ++j) {
            float temp = m[j][i] / m[i][i];
            for (int k = i + 1; k < column_; ++k) {
                m[j][k] -= m[i][k] * temp;
            }
            m[j][i] = 0.00f;
        }
    }
    auto t2 = system_clock::now();//计时结束
    cout << "对齐串行用时\t" << duration_cast<nanoseconds>(t2 - t1).count() << '\n';
    for (size_t i = 0; i < row_; ++i) {
        _aligned_free(m[i]); // 释放每一行的内存
    }
    _aligned_free(m); // 释放行指针数组的内存
}

/*
* SSE并行算法
*/
void SSE_parallel_elimination(float(*m)[column_]) {
    auto t1 = system_clock::now();//计时开始
    for (int i = 0; i < column_ - 1; ++i) {
        for (int j = i + 1; j < row_; ++j) {
            float temp = m[j][i] / m[i][i];
            __m128 temp_vec = _mm_set1_ps(temp); // 将temp广播到四个元素中

            int k;
            // 使用SSE处理能被4整除的部分
            for (k = i + 1; k <= column_ - 4; k += 4) {
                __m128 vec_i = _mm_loadu_ps(&m[i][k]); // 加载m[i][k]到k+3之间的四个元素
                __m128 vec_j = _mm_loadu_ps(&m[j][k]); // 加载m[j][k]到k+3之间的四个元素
                __m128 result = _mm_sub_ps(vec_j, _mm_mul_ps(vec_i, temp_vec));
                _mm_storeu_ps(&m[j][k], result);
            }
            // 处理剩余的元素
            for (; k < column_; ++k) {
                m[j][k] -= m[i][k] * temp;
            }
            m[j][i] = 0.0f;
        }
    }
    auto t2 = system_clock::now();//计时结束
    cout << "SSE用时\t\t" << duration_cast<nanoseconds>(t2 - t1).count() << '\n';
}

void align_SSE_parallel_elimination(float(*m_)[column_]) {
    float** m = static_cast<float**>(_aligned_malloc(row_ * sizeof(float*), 16));
    for (size_t i = 0; i < row_; ++i) {
        m[i] = static_cast<float*>(_aligned_malloc(column_ * sizeof(float), 16));
    }
    for (int i = 0; i < row_; i++)
    {
        for (int j = 0; j < column_; j++)
        {
            m[i][j] = m_[i][j];

        }
        //cout << &m[i][0] << "\t\t" << &m[i][column_ - 1] <<"\t\t";
    }
    auto t1 = system_clock::now();//计时开始
    for (int i = 0; i < column_ - 1; ++i) {
        for (int j = i + 1; j < row_; ++j) {
            float temp = m[j][i] / m[i][i];
            __m128 temp_vec = _mm_set1_ps(temp); // 将temp广播到四个元素中

            int k;
            // 使用SSE处理能被4整除的部分
            for (k = i + 1; k <= column_ - 4; k += 4) {
                __m128 vec_i = _mm_loadu_ps(&m[i][k]); // 加载m[i][k]到k+3之间的四个元素
                __m128 vec_j = _mm_loadu_ps(&m[j][k]); // 加载m[j][k]到k+3之间的四个元素
                __m128 result = _mm_sub_ps(vec_j, _mm_mul_ps(vec_i, temp_vec));
                _mm_storeu_ps(&m[j][k], result);
            }
            // 处理剩余的元素
            for (; k < column_; ++k) {
                m[j][k] -= m[i][k] * temp;
            }
            m[j][i] = 0.0f;
        }
    }
    auto t2 = system_clock::now();//计时结束
    cout << "对齐SSE用时\t" << duration_cast<nanoseconds>(t2 - t1).count() << '\n';
    for (size_t i = 0; i < row_; ++i) {
        _aligned_free(m[i]); // 释放每一行的内存
    }
    _aligned_free(m); // 释放行指针数组的内存
}

/*
* AVX2
*/

void AVX_parallel_elimination(float(*m)[column_]) {
    auto t1 = system_clock::now();
    for (int i = 0; i < column_ - 1; ++i) {
        for (int j = i + 1; j < row_; ++j) {
            float temp = m[j][i] / m[i][i];
            __m256 temp_v = _mm256_set1_ps(temp); // 创建一个包含 8 个 temp 的向量
            for (int k = i + 1; k < column_; k += 8) { // 每次迭代处理 8 个元素
                __m256 m_i_k = _mm256_loadu_ps(&m[i][k]); // 加载 m[i][k] 到 m[i][k+7] 的元素
                __m256 m_j_k = _mm256_loadu_ps(&m[j][k]); // 加载 m[j][k] 到 m[j][k+7] 的元素
                __m256 result = _mm256_fnmadd_ps(m_i_k, temp_v, m_j_k); // 计算 -m[i][k]*temp + m[j][k]
                _mm256_storeu_ps(&m[j][k], result); // 存储结果到 m[j][k] 到 m[j][k+7]
            }
            m[j][i] = 0.00f;
        }
    }
    auto t2 = system_clock::now();
    cout << "AVX用时\t\t" << duration_cast<nanoseconds>(t2 - t1).count() << '\n';
}

void align_AVX_parallel_elimination(float(*m_)[column_]) {
    float** m = static_cast<float**>(_aligned_malloc(row_ * sizeof(float*), 32));
    for (size_t i = 0; i < row_; ++i) {
        m[i] = static_cast<float*>(_aligned_malloc(column_ * sizeof(float), 32));
    }
    for (int i = 0; i < row_; i++)
    {
        for (int j = 0; j < column_; j++)
        {
            m[i][j] = m_[i][j];
        }
    }
    auto t1 = system_clock::now();//计时开始
    for (int i = 0; i < column_ - 1; ++i) {
        for (int j = i + 1; j < row_; ++j) {
            float temp = m[j][i] / m[i][i];
            __m256 temp_v = _mm256_set1_ps(temp); // 创建一个包含 8 个 temp 的向量
            for (int k = i + 1; k + 7 < column_; k += 8) { // 确保至少还有8个元素需要处理
                __m256 m_i_k = _mm256_loadu_ps(&m[i][k]);
                __m256 m_j_k = _mm256_loadu_ps(&m[j][k]);
                __m256 result = _mm256_fnmadd_ps(m_i_k, temp_v, m_j_k);
                _mm256_storeu_ps(&m[j][k], result);
            }
            // 处理剩余的元素
            for (int k = (column_ & ~7) + i + 1; k < column_; k++) {
                m[j][k] -= m[i][k] * temp;
            }
            m[j][i] = 0.00f;
        }
    }
    auto t2 = system_clock::now();//计时结束
    cout << "对齐AVX用时\t" << duration_cast<nanoseconds>(t2 - t1).count() << '\n';
    for (size_t i = 0; i < row_; ++i) {
        _aligned_free(m[i]); // 释放每一行的内存
    }
    _aligned_free(m); // 释放行指针数组的内存
}

int main() {
    srand(time(0));
    for (int row = 0; row < row_; row++) {
        for (int col = 0; col < column_; col++) {
            m[row][col] = float(rand());
        }
    }

//****************************************************************
    reset(m, m_copy);
    simple_serial_elimination(m_copy);
    reset(m, m_copy);
    SSE_parallel_elimination(m_copy);
    reset(m, m_copy);
    align_SSE_parallel_elimination(m_copy);
    reset(m, m_copy);
    AVX_parallel_elimination(m_copy);
    reset(m, m_copy);
    align_AVX_parallel_elimination(m_copy);

    //NEON_parallel_elimination(m);

//****************************************************************
    
    
}
