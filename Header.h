#pragma once
#include<iostream>
#include<vector>
#include <fstream>
using namespace std;
template<typename T>
class vector1
{
protected:
	int n;
public:
	vector<T> cont;
	vector1(int n=1);
	vector1(const vector1& a);
	vector1(ifstream & stream);
	int get_size() const;
	void operator=(T a);
	vector1 operator+(vector1& a);
	vector1 operator -(vector1 &a);
	vector1 operator *(T a);
	friend vector1 operator *(vector1&b, T& a);
	T& operator [](int m);
	double scalar_mul(vector1&a);
	double norm();
	vector<vector<T>> out_mul(vector1&a);
	double angle(vector1&a);
	T max_norm();
	void file_bout(ofstream & stream);
	vector1 file_bin(ifstream& stream);
//	ofstream& operator<<(ofstream & stream, vector1<T> s);

	~vector1();

};


template<typename T>
class matrix
{
protected:
	int m, n;
public:
	vector<vector1<T>> cont;
	matrix();
	matrix(int m, int n);
	matrix(const matrix &a);
	matrix(ifstream & stream);
	int get_width() const;
	int get_height() const;
	void operator =(T a);
	matrix operator +(matrix &a);
	matrix operator -(matrix &a);
	matrix operator *(matrix& a);
    matrix operator *(T a);
	T& operator ()(const int& m,const int& n);
	vector1<T>& operator[](int m);
	matrix col(int m);
	matrix Hadamar_mul(matrix&a);
	T trace();
	T det();
	double norm();
	int rk();
	matrix<T> inverse_matrix();
	matrix<T> transpose();
	void file_bout(ofstream & stream);
	matrix file_bin(ifstream& stream);
	void append_row(matrix<T> a, int k);
	void append_col(matrix<T> a, int k);
	virtual ~matrix();
};



template<typename T>
class single_matrix :public matrix<T>
{
public:
	single_matrix(int m):matrix<T>(m,m) {
		for (int i = 0; i < n; i++) {
			cont[i] = 0;
			cont[i][i] = 1;
		}
	}
};


template<typename T>
class diag_matrix :public matrix<T>
{
public:
	diag_matrix(int m): matrix<T>(m, m) {}
	T& operator ()(int& m, int& n) {
		if (m != n)
			return 0;
		else
			return cont[m][n];
	}
};

template<typename T>
class vt_matrix :public matrix<T>
{
	
public:
	vt_matrix(int m): matrix<T>(m, m) {}
	T& operator ()(int& m, int& n) {
		if (m > n)
			return 0;
		else
			return cont[m][n];
	}
};

template<typename T>
class nt_matrix :public matrix<T>
{
public:
	nt_matrix(int m) : matrix<T>(m, m) {}
	T& operator ()(int& m, int& n) {
		if (m >= n)
			return cont[m][n];
		else
			return 0;
	}
};

template<typename T>
class sym_matrix :public matrix<T>
{
public:
	sym_matrix(int m) : matrix<T>(m, m) {}
	T& operator ()(int& m, int& n) {
		if (m > n)
			return cont[n][m];
		else
			return cont[m][n];
	}
};
