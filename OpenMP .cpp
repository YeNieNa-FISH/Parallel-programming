#include<omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <iostream>
#include <chrono>
#include <ctime>
#include <xmmintrin.h> //SSE
#include <emmintrin.h> //SSE2
#include <pmmintrin.h> //SSE3
#include <tmmintrin.h> //SSSE3
#include <smmintrin.h> //SSE4.1
#include <nmmintrin.h> //SSSE4.2
using namespace std;
using namespace chrono;


int n=2000;//���ݹ�ģ
float** m_data ;//��������
const int thread_count=7;//�߳���

//��ʼ����������
void init() {
	for (int i = 0; i < n; i++) {
		m_data[i] = new float[n];
	}
	srand((int)time(0));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			m_data[i][j] = rand() % 100;
		}
	}
}

//��ӡ����
void printMatrix(float**m_data) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << m_data[i][j] << " ";
		}
		cout << endl;
	}
}

//��ͨ��˹��ȥ�����㷨
float** normal_Gauss(float** m_data) {
	for (int k = 0; k < n; k++) {
		for (int j = k + 1; j < n; j++) {
			m_data[k][j] = m_data[k][j] / m_data[k][k];
		}
		m_data[k][k] = 1.0;
		for (int i = k + 1; i < n; i++) {
			for (int j = k + 1; j < n; j++) {
				m_data[i][j] = m_data[i][j] - m_data[i][k] * m_data[k][j];
			}
			m_data[i][k] = 0;
		}
	}
	return m_data;
}



//test��dynamicģʽ ���л��֣�(����黮��)
float** dynamic_row_block_test(float** m_data)
{
	int k;
	int i;
	int j;
	float temp;
	//����ѭ��֮�ⴴ���̣߳������̷߳����������٣�ע�⹲��������˽�б���������
	#pragma omp parallel , num_threads(thread_count), private(i, j, k, temp)
	for ( k = 0; k < n; k++)
	{
		//���в��֣�Ҳ���Կ��ǲ��л�
		#pragma omp single
		temp = m_data[k][k];
		for ( j = k + 1; j < n; j++)
		{
			m_data[k][j] = m_data[k][j] / temp;
		}
		m_data[k][k] = 1.0;
		//���в��֣�ʹ���л���
		#pragma omp for schedule(dynamic,n/thread_count)
		for ( i = k +1 ; i < n; i++)
		{
			temp = m_data[i][k];
			for ( j = k + 1; j < n; ++j)
			{
				m_data[i][j] = m_data[i][j] -temp * m_data[k][j];
			}
			m_data[i][k] = 0.0;
		}
		//�뿪forѭ��ʱ�������߳�Ĭ��ͬ����������һ�еĴ���
	}
	return m_data;
}


//test:dynamicģʽ + SSE ���л��֣�������黮�֣�
float** dynamic_sse_row_block_test(float** m_data)
{
	int k;
	int i;
	int j;
	float temp;
	//����ѭ��֮�ⴴ���̣߳������̷߳����������٣�ע�⹲��������˽�б���������
#pragma omp parallel , num_threads(thread_count), private(i, j, k, temp)
	for (k = 0; k < n; k++)
	{
		//���в��֣�Ҳ���Կ��ǲ��л�
#pragma omp single
		temp = m_data[k][k];
		__m128 t1, t0;
		t0 = _mm_set_ps(temp,temp,temp,temp);
		for ( j = k + 1; j + 4 < n; j = j + 4)
		{
			t1 = _mm_loadu_ps(m_data[k] + j);
			t1 = _mm_div_ps(t1, t0);
			_mm_store_ps(m_data[k] + j, t1);
		}
		for (int j = (n / 4) * 4; j < n; j++)
		{
			m_data[k][j] = m_data[k][j] / temp;
		}
		m_data[k][k] = float(1.0);
		//���в��֣�ʹ���л���
#pragma omp for schedule(dynamic,n/thread_count)
		__m128 t2, t3, t4, vx;
		for ( i = k+1; i < n; i++)
		{
			temp = m_data[i][k];
			t2 = _mm_set_ps(temp,temp,temp,temp);
			for (int j = k + 1; j + 4 < n; j = j + 4)
			{
				t3 = _mm_loadu_ps(m_data[k] + j);
				t4 = _mm_loadu_ps(m_data[i] + j);
				vx = _mm_mul_ps(t2, t3);
				t4 = _mm_sub_ps(t4, vx);
				_mm_store_ps(m_data[i] + j, t4);
			}
			for (int j = (n / 4) * 4; j < n; j++)
			{
				m_data[i][j] = m_data[i][j] - temp * m_data[k][j];
			}
			m_data[i][k] = 0;
		}
	}
	return m_data;
}


//test: staticģʽ ���л��֣�������黮�֣�
float** static_row_block_test(float** m_data)
{
	int k;
	int i;
	int j;
	float temp;
	//����ѭ��֮�ⴴ���̣߳������̷߳����������٣�ע�⹲��������˽�б���������
#pragma omp parallel , num_threads(thread_count), private(i, j, k, temp)
	for (k = 0; k < n; k++)
	{
		//���в��֣�Ҳ���Կ��ǲ��л�
#pragma omp single
		temp = m_data[k][k];
		for (j = k + 1; j < n; j++)
		{
			m_data[k][j] = m_data[k][j] / temp;
		}
		m_data[k][k] = 1.0;
		//���в��֣�ʹ���л���
#pragma omp for schedule(static,n/thread_count)
		for (i = k + 1; i < n; i++)
		{
			temp = m_data[i][k];
			for (j = k + 1; j < n; ++j)
			{
				m_data[i][j] = m_data[i][j] - temp * m_data[k][j];
			}
			m_data[i][k] = 0.0;
		}
		//�뿪forѭ��ʱ�������߳�Ĭ��ͬ����������һ�еĴ���
	}
	return m_data;
}


//test��staticģʽ + SSE ���л��֣�������黮�֣�
float** static_sse_row_block_test(float** m_data)
{
	int k;
	int i;
	int j;
	float temp;
	//����ѭ��֮�ⴴ���̣߳������̷߳����������٣�ע�⹲��������˽�б���������
#pragma omp parallel , num_threads(thread_count), private(i, j, k, temp)
	for (k = 0; k < n; k++)
	{
		//���в��֣�Ҳ���Կ��ǲ��л�
#pragma omp single
		temp = m_data[k][k];
		__m128 t1, t0;
		t0 = _mm_set_ps(temp, temp, temp, temp);
		for (j = k + 1; j + 4 < n; j = j + 4)
		{
			t1 = _mm_loadu_ps(m_data[k] + j);
			t1 = _mm_div_ps(t1, t0);
			_mm_store_ps(m_data[k] + j, t1);
		}
		for (int j = (n / 4) * 4; j < n; j++)
		{
			m_data[k][j] = m_data[k][j] / temp;
		}
		m_data[k][k] = float(1.0);
		//���в��֣�ʹ���л���
#pragma omp for schedule(static,n/thread_count)
		__m128 t2, t3, t4, vx;
		for (i = k + 1; i < n; i++)
		{
			temp = m_data[i][k];
			t2 = _mm_set_ps(temp, temp, temp, temp);
			for (int j = k + 1; j + 4 < n; j = j + 4)
			{
				t3 = _mm_loadu_ps(m_data[k] + j);
				t4 = _mm_loadu_ps(m_data[i] + j);
				vx = _mm_mul_ps(t2, t3);
				t4 = _mm_sub_ps(t4, vx);
				_mm_store_ps(m_data[i] + j, t4);
			}
			for (int j = (n / 4) * 4; j < n; j++)
			{
				m_data[i][j] = m_data[i][j] - temp * m_data[k][j];
			}
			m_data[i][k] = 0;
		}
	}
	return m_data;
}


//test: staticģʽ ���л��֣���ѭ�����֣�
float** static_row_loop_test(float**m_data)
{
	int k;
	int i;
	int j;
	float temp;
	//����ѭ��֮�ⴴ���̣߳������̷߳����������٣�ע�⹲��������˽�б���������
#pragma omp parallel , num_threads(thread_count), private(i, j, k, temp)
	for (k = 0; k < n; k++)
	{
		//���в��֣�Ҳ���Կ��ǲ��л�
#pragma omp single
		temp = m_data[k][k];
		for (j = k + 1; j < n; j++)
		{
			m_data[k][j] = m_data[k][j] / temp;
		}
		m_data[k][k] = 1.0;
		//���в��֣�ʹ���л���
#pragma omp for schedule(static,1)
		for (i = k + 1; i < n; i++)
		{
			temp = m_data[i][k];
			for (j = k + 1; j < n; ++j)
			{
				m_data[i][j] = m_data[i][j] - temp * m_data[k][j];
			}
			m_data[i][k] = 0.0;
		}
		//�뿪forѭ��ʱ�������߳�Ĭ��ͬ����������һ�еĴ���
	}
	return m_data;
}


//test: staticģʽ (�л���) ������黮�֣�
float** static_col_block_test(float** m_data)
{
	int k;
	int i;
	int j;
	float temp;
	//����ѭ��֮�ⴴ���̣߳������̷߳����������٣�ע�⹲��������˽�б���������
#pragma omp parallel , num_threads(thread_count), private(i, j, k, temp)
	for (k = 0; k < n; k++)
	{
		//���в��֣�Ҳ���Կ��ǲ��л�
#pragma omp single
		temp = m_data[k][k];
		for (j = k + 1; j < n; j++)
		{
			m_data[k][j] = m_data[k][j] / temp;
		}
		m_data[k][k] = 1.0;
		//���в��֣�ʹ���л���

		for (i = k + 1; i < n; i++)
		{
#pragma omp for schedule(static,n/thread_count)
			temp = m_data[i][k];
			for (j = k + 1; j < n; ++j)
			{
				m_data[i][j] = m_data[i][j] - temp * m_data[k][j];
			}
			m_data[i][k] = 0.0;
		}
		//�뿪forѭ��ʱ�������߳�Ĭ��ͬ����������һ�еĴ���
	}
	return m_data;
}


//test��guidedģʽ  ���л��֣� ������黮��--��С��Ϊn/thread_count��
float** guided_row_block_test(float** m_data)
{
	int k;
	int i;
	int j;
	float temp;
	//����ѭ��֮�ⴴ���̣߳������̷߳����������٣�ע�⹲��������˽�б���������
#pragma omp parallel , num_threads(thread_count), private(i, j, k, temp)
	for (k = 0; k < n; k++)
	{
		//���в��֣�Ҳ���Կ��ǲ��л�
#pragma omp single
		temp = m_data[k][k];
		for (j = k + 1; j < n; j++)
		{
			m_data[k][j] = m_data[k][j] / temp;
		}
		m_data[k][k] = 1.0;
		//���в��֣�ʹ���л���
#pragma omp for schedule(guided,n/thread_count)
		for (i = k + 1; i < n; i++)
		{
			temp = m_data[i][k];
			for (j = k + 1; j < n; ++j)
			{
				m_data[i][j] = m_data[i][j] - temp * m_data[k][j];
			}
			m_data[i][k] = 0.0;
		}
		//�뿪forѭ��ʱ�������߳�Ĭ��ͬ����������һ�еĴ���
	}
	return m_data;
}



int main() {

	int m[9] = { 100,200,300,400,500,1000,2000,3000,4000 };
	for (int i = 0; i < 9; i++) {
		n = m[i];
        cout << n << ": ";
        long long sum = 0;
		for(int j =0;j<10;j++){
            m_data = new float* [n];//��������

            //������ȷ��

            init();
            auto t1 = system_clock::now();//��ʱ��ʼ

            //normal_Gauss(m_data);//�����㷨
            //dynamic_row_block_test(m_data);			//��̬  �л���  �黮��
            //dynamic_sse_row_block_test(m_data);		//��̬  SSE  �л���  �黮��
            //static_row_block_test(m_data);			//��̬  �л���  �黮��
            //static_sse_row_block_test(m_data);		//��̬  SSE  �л���  �黮��
            //static_row_loop_test(m_data);				//��̬  �л���  ѭ������
            //static_col_block_test(m_data);			//��̬  �л���  �黮��
            guided_row_block_test(m_data);			//guided �л���  �黮�֣���С��Ϊn/thread_count)--------���ز���

            auto t2 = system_clock::now();//��ʱ����
            //cout << n << "��ʱ\t\t" << duration_cast<nanoseconds>(t2 - t1).count() << '\n';
            sum += duration_cast<nanoseconds>(t2 - t1).count();

            cout << duration_cast<nanoseconds>(t2 - t1).count() << ' ';
		}
		cout<<sum / 10<<endl;

		//�ͷŶ�̬�ڴ�
		for (int i = 0; i < n; i++) {
			delete[] m_data[i];
		}
		delete[] m_data;
	}

	return 0;
}