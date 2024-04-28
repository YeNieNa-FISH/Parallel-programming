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
* �����㷨
*/
void simple_serial_elimination(float(*m)[column_]) {
    auto t1 = system_clock::now();//��ʱ��ʼ
    for (int i = 0; i < column_ - 1; ++i) {
        for (int j = i + 1; j < row_; ++j) {
            float temp = m[j][i] / m[i][i];
            for (int k = i + 1; k < column_; ++k) {
                m[j][k] -= m[i][k] * temp;
            }
            m[j][i] = 0.00f;
        }
    }
    auto t2 = system_clock::now();//��ʱ����
    cout << "������ʱ\t" << duration_cast<nanoseconds>(t2 - t1).count() << '\n';
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
    auto t1 = system_clock::now();//��ʱ��ʼ
    for (int i = 0; i < column_ - 1; ++i) {
        for (int j = i + 1; j < row_; ++j) {
            float temp = m[j][i] / m[i][i];
            for (int k = i + 1; k < column_; ++k) {
                m[j][k] -= m[i][k] * temp;
            }
            m[j][i] = 0.00f;
        }
    }
    auto t2 = system_clock::now();//��ʱ����
    cout << "���봮����ʱ\t" << duration_cast<nanoseconds>(t2 - t1).count() << '\n';
    for (size_t i = 0; i < row_; ++i) {
        _aligned_free(m[i]); // �ͷ�ÿһ�е��ڴ�
    }
    _aligned_free(m); // �ͷ���ָ��������ڴ�
}

/*
* SSE�����㷨
*/
void SSE_parallel_elimination(float(*m)[column_]) {
    auto t1 = system_clock::now();//��ʱ��ʼ
    for (int i = 0; i < column_ - 1; ++i) {
        for (int j = i + 1; j < row_; ++j) {
            float temp = m[j][i] / m[i][i];
            __m128 temp_vec = _mm_set1_ps(temp); // ��temp�㲥���ĸ�Ԫ����

            int k;
            // ʹ��SSE�����ܱ�4�����Ĳ���
            for (k = i + 1; k <= column_ - 4; k += 4) {
                __m128 vec_i = _mm_loadu_ps(&m[i][k]); // ����m[i][k]��k+3֮����ĸ�Ԫ��
                __m128 vec_j = _mm_loadu_ps(&m[j][k]); // ����m[j][k]��k+3֮����ĸ�Ԫ��
                __m128 result = _mm_sub_ps(vec_j, _mm_mul_ps(vec_i, temp_vec));
                _mm_storeu_ps(&m[j][k], result);
            }
            // ����ʣ���Ԫ��
            for (; k < column_; ++k) {
                m[j][k] -= m[i][k] * temp;
            }
            m[j][i] = 0.0f;
        }
    }
    auto t2 = system_clock::now();//��ʱ����
    cout << "SSE��ʱ\t\t" << duration_cast<nanoseconds>(t2 - t1).count() << '\n';
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
    auto t1 = system_clock::now();//��ʱ��ʼ
    for (int i = 0; i < column_ - 1; ++i) {
        for (int j = i + 1; j < row_; ++j) {
            float temp = m[j][i] / m[i][i];
            __m128 temp_vec = _mm_set1_ps(temp); // ��temp�㲥���ĸ�Ԫ����

            int k;
            // ʹ��SSE�����ܱ�4�����Ĳ���
            for (k = i + 1; k <= column_ - 4; k += 4) {
                __m128 vec_i = _mm_loadu_ps(&m[i][k]); // ����m[i][k]��k+3֮����ĸ�Ԫ��
                __m128 vec_j = _mm_loadu_ps(&m[j][k]); // ����m[j][k]��k+3֮����ĸ�Ԫ��
                __m128 result = _mm_sub_ps(vec_j, _mm_mul_ps(vec_i, temp_vec));
                _mm_storeu_ps(&m[j][k], result);
            }
            // ����ʣ���Ԫ��
            for (; k < column_; ++k) {
                m[j][k] -= m[i][k] * temp;
            }
            m[j][i] = 0.0f;
        }
    }
    auto t2 = system_clock::now();//��ʱ����
    cout << "����SSE��ʱ\t" << duration_cast<nanoseconds>(t2 - t1).count() << '\n';
    for (size_t i = 0; i < row_; ++i) {
        _aligned_free(m[i]); // �ͷ�ÿһ�е��ڴ�
    }
    _aligned_free(m); // �ͷ���ָ��������ڴ�
}

/*
* AVX2
*/

void AVX_parallel_elimination(float(*m)[column_]) {
    auto t1 = system_clock::now();
    for (int i = 0; i < column_ - 1; ++i) {
        for (int j = i + 1; j < row_; ++j) {
            float temp = m[j][i] / m[i][i];
            __m256 temp_v = _mm256_set1_ps(temp); // ����һ������ 8 �� temp ������
            for (int k = i + 1; k < column_; k += 8) { // ÿ�ε������� 8 ��Ԫ��
                __m256 m_i_k = _mm256_loadu_ps(&m[i][k]); // ���� m[i][k] �� m[i][k+7] ��Ԫ��
                __m256 m_j_k = _mm256_loadu_ps(&m[j][k]); // ���� m[j][k] �� m[j][k+7] ��Ԫ��
                __m256 result = _mm256_fnmadd_ps(m_i_k, temp_v, m_j_k); // ���� -m[i][k]*temp + m[j][k]
                _mm256_storeu_ps(&m[j][k], result); // �洢����� m[j][k] �� m[j][k+7]
            }
            m[j][i] = 0.00f;
        }
    }
    auto t2 = system_clock::now();
    cout << "AVX��ʱ\t\t" << duration_cast<nanoseconds>(t2 - t1).count() << '\n';
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
    auto t1 = system_clock::now();//��ʱ��ʼ
    for (int i = 0; i < column_ - 1; ++i) {
        for (int j = i + 1; j < row_; ++j) {
            float temp = m[j][i] / m[i][i];
            __m256 temp_v = _mm256_set1_ps(temp); // ����һ������ 8 �� temp ������
            for (int k = i + 1; k + 7 < column_; k += 8) { // ȷ�����ٻ���8��Ԫ����Ҫ����
                __m256 m_i_k = _mm256_loadu_ps(&m[i][k]);
                __m256 m_j_k = _mm256_loadu_ps(&m[j][k]);
                __m256 result = _mm256_fnmadd_ps(m_i_k, temp_v, m_j_k);
                _mm256_storeu_ps(&m[j][k], result);
            }
            // ����ʣ���Ԫ��
            for (int k = (column_ & ~7) + i + 1; k < column_; k++) {
                m[j][k] -= m[i][k] * temp;
            }
            m[j][i] = 0.00f;
        }
    }
    auto t2 = system_clock::now();//��ʱ����
    cout << "����AVX��ʱ\t" << duration_cast<nanoseconds>(t2 - t1).count() << '\n';
    for (size_t i = 0; i < row_; ++i) {
        _aligned_free(m[i]); // �ͷ�ÿһ�е��ڴ�
    }
    _aligned_free(m); // �ͷ���ָ��������ڴ�
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