/*
#include<iostream>
#include<vector>
#include<chrono>

using namespace std;
using namespace chrono;

int main() {
	vector<vector<int>> matrix;
	vector<int> vec;
	int n;
	//int temp;
	cout << "����N" << endl;
	cin >> n;
	cout << "�������" << endl;
	for (int i = 0; i < n; i++) {
		vector<int> New;
		matrix.push_back(New);
		for (int j = 0; j < n; j++) {
			//cin >> temp;
			matrix[i].push_back(1);
		}
	}
	cout << "��������" << endl;
	for (int i = 0; i < n; i++) {
		//cin >> temp;
		vec.push_back(1);
	}

	typedef high_resolution_clock Clock;
	auto t1 = Clock::now();//��ʱ��ʼ

	int* sum1 = new int[n];
	int* sum2 = new int[n];
	int flag1 = 0;
	int flag2 = 0;

	//�㷨һ
	for (int x = 0; x < 100; x++) {
		for (int i = 0; i < n; i++) {
			sum1[i] = 0;
			for (int j = 0; j < n; j++) {
				sum1[i] += matrix[j][i] * vec[j];
				//flag1++;
			}
		}
	}
	
	auto t2 = Clock::now();//��ʱ����
	cout << "�㷨һ��ʱ��" << duration_cast<nanoseconds>(t2 - t1).count() << '\n';

	auto t3 = Clock::now();//��ʱ��ʼ

	for (int x = 0; x < 100; x++) {
		for (int i = 0; i < n; i++) sum2[i] = 0;
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < n; i++) {
				sum2[i] += matrix[j][i] * vec[j];
				//flag2++;
			}
		}
	}
	
	auto t4 = Clock::now();//��ʱ����
	cout << "�㷨����ʱ��" << duration_cast<nanoseconds>(t4 - t3).count() << '\n';
	//cout << "1:" << flag1 << " ;2:" << flag2 << endl;
	
}*/