
#include<iostream>
#include<vector>
#include<chrono>

using namespace std;
using namespace chrono;

const int counts = 100;

void su(int n, vector<int> v) {
	if (n == 1) return;
	else {
		for (int i = 0; i < n / 2; i++)
			v[i] += v[n - i - 1];
		n = n / 2;
		su(n, v);
	}
}

int main() {
	vector<int> vec;
	
	int N = 0;
	int n = 2;
	cout << "输入N：";
	cin >> N;
	for (int i = 0; i < N; i++) {
		n = n * 2;
	}
	//int* vec = new int[n];
	for (int i = 0; i < n; i++) {
		vec.push_back(1);
		//vec[i] = 1;
	}
	
	int sum = 0;
	typedef high_resolution_clock Clock;
	auto t1 = system_clock::now();//计时开始
	for (int x = 0; x < counts; x++) {
		for (int i = 0; i < n; i++) {
			sum += vec[i];
		}
	}
	auto t2 = system_clock::now();//计时结束
	cout << "算法一用时：" << duration_cast<nanoseconds>(t2 - t1).count() << '\n';

	//****************************************************************
	
	auto t3 = system_clock::now();//计时开始
	for (int x = 0; x < counts; x++) {
		int sum1 = 0;
		int sum2 = 0;
		for (int i = 0; i < n; i += 2) {
			sum1 += vec[i];
			sum2 += vec[i + 1];
		}
	}
	auto t4 = system_clock::now();//计时结束
	cout << "算法二用时：" << duration_cast<nanoseconds>(t4 - t3).count() << '\n';

	

	auto t5 = system_clock::now();
	for (int x = 0; x < counts; x++) {
		for (int m = n; m > 1; m /= 2) {
			for (int i = 0; i < m / 2; i++) {
				vec[i] = vec[i * 2] + vec[i * 2 + 1];
			}
		}
	}
	auto t6 = system_clock::now();//计时结束
	cout << "算法三用时：" << duration_cast<nanoseconds>(t6 - t5).count() << '\n';

	
	
	auto t7 = system_clock::now();
	for (int x = 0; x < counts; x++) {
		su(n, vec);
	}
	auto t8 = system_clock::now();//计时结束
	cout << "算法四用时：" << duration_cast<nanoseconds>(t8 - t7).count() << '\n';
	
}