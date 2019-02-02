#pragma once
#include"Header.h"
#include<math.h>
#include <fstream>
#include<string>
using namespace std;

template<typename T>
matrix<T>::matrix()
{
	cont = vector<vector1<T>>(1);
	for (int i = 0; i < m; i++)
		cont[i] = vector1<T>(1); 
}
template<typename T>
matrix<T>::matrix(int m, int n) :m(m), n(n)
{
    vector1<T> mat;
	cont = vector<vector1<T>>(m);
	for (int i = 0; i < m; i++)
		cont[i]= vector1<T>(n);
	
}

template<typename T>
matrix<T>::matrix(const matrix<T> & a) :m(a.get_height()), n(a.get_width())
{
	cont = vector<vector1<T>>(m);
		for (int i = 0; i < m; i++)
			cont[i] = vector1<T>(n);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			cont[i][j] = a.cont[i].cont[j];
}

template<typename T>
matrix<T>::matrix(ifstream & stream)
{
	if (!stream.is_open())
		throw std::invalid_argument("sorry");

	while (!stream.eof()) {
		vector1<T> a(stream);
		if (a.get_size() == 1)
			break;
		cont.push_back(a);
		m++;
	}
	n = cont[0].get_size();
}

template<typename T>
int matrix<T>::get_width() const
{
	return n;
}

template<typename T>
int matrix<T>::get_height() const
{
	return m;
}

template<typename T>
void matrix<T>::operator=(T a)
{
	for (int i = 0; i < this->get_height(); i++)
		for (int j = 0; j < this->get_width(); j++)
			(*this)(i,j) = a;

}

template<typename T>
matrix<T> matrix<T>::operator+(matrix<T> & a)
{
	if (this->get_width() != a.get_width() || this->get_height() != a.get_height())
		throw std::invalid_argument("sorry");
	matrix<T> m(this->get_height(), this->get_width());
	for (int i = 0; i < m.get_height(); i++)
		for (int j = 0; j < m.get_width(); j++)
			m(i, j) = a(i, j) + (*this)(i, j);
	return m;
}

template<typename T>
matrix<T> matrix<T>::operator-(matrix<T> & a)
{
	if (this->get_width() != a.get_width() || this->get_height() != a.get_height())
		throw std::invalid_argument("sorry");
	matrix<T> m(this->get_height(), this->get_width());
	for (int i = 0; i < m.get_height(); i++)
		for (int j = 0; j < m.get_width(); j++)
			m(i, j) = a(i, j) - (*this)(i, j);
	return m;
}

template<typename T>
matrix<T> matrix<T>::operator*(matrix<T> & a)
{
	if (this->get_width() != a.get_height())
		throw std::invalid_argument("sorry");
	matrix<T> m(this->get_height(), a.get_width());
	m = 0;
	for (int i = 0; i < m.get_height(); i++)
		for (int j = 0; j < m.get_width(); j++)
			for (int k = 0; k < this->get_width(); k++)
				m(i, j) += (*this)(i, k) * a(k, j);
	return m;
}
template<typename T>
matrix<T> matrix<T>::operator*(T  a)
{
	matrix<T> m((*this));
	for (int i = 0; i < this->get_height(); i++)
		for (int j = 0; j < this->get_width(); j++)
			m(i, j)*=a;
	return m;
}
template<typename T>
matrix<T> operator*(T  a, matrix<T> & b)
{
	return b*a;
}
template<typename T>
T& matrix<T>::operator()(const int& m,const int& n)
{
	if(matrix<T>::m<=m|| matrix<T>::n <= n)
		throw std::invalid_argument("sorry");
	return this->cont[m][n];
}
template<typename T>
vector1<T> & matrix<T>::operator[](int m)
{
	if (matrix<T>::m <= m)
		throw std::invalid_argument("sorry");
	return cont[m];
}
template<typename T>
matrix<T> matrix<T>::col(int m)
{
	if (matrix<T>::n <= m)
		throw std::invalid_argument("sorry");
    matrix<T> new_v(matrix<T>::m,1);
	for (int i = 0; i <matrix<T>::m; i++) {
		new_v(i, 0) = cont[i][m];
	}
	return new_v;
}
template<typename T>
matrix<T> matrix<T>::Hadamar_mul(matrix & a)
{
	if (m != a.get_width() || n != a.get_height())
		throw std::invalid_argument("sorry");
	matrix<T> arr(m, n);
	arr = 0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			arr(i, j) = (*this)(i, j)*a(i, j);

			return arr;
}
template<typename T>
T matrix<T>::trace()
{
	T sum = 0;
	for (int i = 0; i < (m > n ? n : m); i++)
		sum += (*this)(i, i);
	return sum;

}
template<typename T>
T matrix<T>::det()
{
	if (m != n)
		throw std::invalid_argument("sorry");
	matrix<T> arr(*this);
	for (int i = 0; i < n; i++) {
		if (arr(i, i) == 0) {
			for (int k = i + 1; k < n; k++)
				if (arr(k, i) != 0)
					{
					arr[i] = arr[i] + arr[k];
					break;
					}
			if (arr(i, i) == 0)
				return 0;
		}
		for (int j = i + 1; j < n; j++) 
			if (arr(i,i) != 0)
				arr[j] = arr[j] - (arr[i] * (arr(j, i) / arr(i, i)));
	}
	T sum = 1;
	for (int i = 0; i < n; i++)
		sum *= arr(i, i);
	return sum;
}
template<typename T>
double matrix<T>::norm()
{
	double sum = 0;
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			sum += pow(cont[i][j],2);
	return sqrt(sum);
}
template<typename T>
int matrix<T>::rk()
{
	int rk = m < n ? m : n;
	int mx = m > n ? n : m;
	int numb = 0;
	matrix<T> arr(*this);
	for (int i = 0; i < n&&numb<mx; i++, numb++)
	{
		if (arr(numb, i) == 0){
			for (int k = numb + 1; k < m; k++) 
				if (arr(k, i) != 0)
				{
					arr[i] = arr[i] + arr[k];
					break;
				}
			if (arr(numb, i) == 0)
			{
				rk--; numb--; continue;
			}
			}
		for (int j = i + 1; j < m; j++)
			if (arr(numb, i) != 0)
				arr[j] = arr[j] - (arr[numb] * (arr(j, i) / arr(numb, i)));
	}
	return rk;
}
template<typename T>
matrix<T> matrix<T>::inverse_matrix()
{
	if(m!=n)
		throw std::invalid_argument("sorry");
	matrix<T> arr(*this);
	single_matrix<T> ed(n);
	for (int i = 0; i < n; i++)
	{
		if (arr(i, i) == 0) 
		{
			for (int k = i + 1; k < n; k++)
				if (arr(k, i) != 0)
				{
					arr[i] = arr[i] + arr[k];
					ed[i] = ed[i] + ed[k];
					break;
				}
			if (arr(i, i) == 0)
				throw std::invalid_argument("sorry");
		}
		for (int j = i + 1; j < n; j++)
			 {
				T a = arr(j, i) / arr(i, i);
				arr[j] = arr[j] - ((vector1<T>)arr[i] * (T)a);
				ed[j] = ed[j] - ((vector1<T>)ed[i] * (T)a);
			}
	}

	for (int i = n-1; i >= 0; i--) 
		for (int j = i -1; j >=0; j--) {
			T a = arr(j, i) / arr(i, i);
			arr[j] = arr[j] - ((vector1<T>)arr[i] * (T)a);
			ed[j] = ed[j] - ((vector1<T>)ed[i] * (T)a);
			}
	for (int i = 0; i < n; i++)
		ed[i] = ed[i] * (1 / arr(i, i));

	return ed;
}
template<typename T>
matrix<T> matrix<T>::transpose()
{
	matrix<T> trans(n, m);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			trans(i, j) = (*this)(j, i);
	return trans;
}
template<typename T>
void matrix<T>::file_bout(ofstream & stream)
{
	if (!stream.is_open())
		throw std::invalid_argument("sorry");
	stream.write(reinterpret_cast<char*>(&m), sizeof(int));
	stream.write(reinterpret_cast<char*>(&n), sizeof(int));
	//stream.write(m<< n;
	for (int i = 0; i < m; i++)
		for(int j=0;j<n;j++)
		stream.write(reinterpret_cast<char*>(&cont[i][j]),sizeof(float));
}
template<typename T>
matrix<T> matrix<T>::file_bin(ifstream & stream)
{
	int a,b;
	if (!stream.is_open())
		throw std::invalid_argument("sorry");
	stream.read(reinterpret_cast<char*>(&a), sizeof(int));
	stream.read(reinterpret_cast<char*>(&b), sizeof(int));
	matrix<T> s(a, b);
	for (int i = 0; i < a; i++)
		for(int j=0;j<b;j++)
			stream.read(reinterpret_cast<char*>(&s[i][j]), sizeof(float));
	return s;
}
template<typename T>
void matrix<T>::append_row(matrix<T> a, int k)
{
	for (int i = 0; i < n; i++)
		cont[k][i] = a(0, i);
}
template<typename T>
void matrix<T>::append_col(matrix<T> a, int k)
{
	for (int i = 0; i < m; i++)
		cont[i][k] = a(i, 0);
}
template<typename T>
matrix<T>::~matrix()
{
	//destuctor
}






template<typename T>
vector1<T>::vector1(int n) :n(n)
{
	cont = vector<T>(n);
}

template<typename T>
vector1<T>::vector1(const vector1 & a) :n(a.get_size())
{
	cont= vector<T>(n);
		for (int i = 0; i < n; i++)
			(*this)[i] = a.cont[i];

}

template<typename T>
vector1<T>::vector1(ifstream & stream)
{
	string c;
	char c1='\0' ;
	int i = 0;
	n = 0;
	if (!stream.is_open())
		throw std::invalid_argument("sorry");
	//ios::pos_type pos = stream.tellg();
	ios::pos_type pos = stream.tellg();
	while(c==""&&!stream.eof())
	getline(stream, c);
	while ((i = c.find('\t'))!=string::npos) {
		c[i] = 'a';
		n++;
	}
	n++;
	cont = vector<T>(n);
		stream.seekg(pos, ios::beg);
		for (int i = 0; i < n; i++)
			stream >> cont[i];
}

template<typename T>
int vector1<T>::get_size() const
{
	return n;
}

template<typename T>
void vector1<T>::operator=(T a)
{
	for (int i = 0; i < n; i++)
		cont[i] = a;
}

/*template<typename T>
void vector1<T>::operator=(vector1 & a)
{
	for (int i = 0; i < n; i++)
		cont[i] = a[i];
}
*/
template<typename T>
vector1<T> vector1<T>::operator+(vector1 & a)
{
	if (n != a.get_size())
		throw std::invalid_argument("sorry");
	vector1<T> arr(n);
	for (int i = 0; i < n; i++)
		arr[i] = (*this)[i] + a[i];
		return arr;
}

template<typename T>
vector1<T> vector1<T>::operator-(vector1 & a)
{
	if (n != a.get_size())
		throw std::invalid_argument("sorry");
	vector1<T> arr(n);
	for (int i = 0; i < n; i++)
		arr[i] = (*this)[i] - a[i];
		return arr;
}
template<typename T>
T & vector1<T>::operator[](int m)
{
	return cont[m];
}
template<typename T>
double vector1<T>::scalar_mul(vector1 & a)
{
	if(n!=a.get_size())
		throw std::invalid_argument("sorry");
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += cont[i] * a[i];
	return sum;
}
template<typename T>
double vector1<T>::norm()
{
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += cont[i] * cont[i];
	return sqrt(sum);
}
template<typename T>
vector<vector<T>> vector1<T>::out_mul(vector1 & a)
{
	vector<vector<T>>  matr(n, vector<T>(a.get_size()));
	for (int i = 0; i < n; i++)
		for (int j = 0; j < a.get_size(); j++)
			matr[i][j] = cont[i] * a[j];
	return matr;
}
template<typename T>
double vector1<T>::angle(vector1 & a)
{
	double c = this->norm(), b = a.norm();
	if(c==0||b==0)
		throw std::invalid_argument("sorry");
	double sum=this->scalar_mul(a)/(c*b);
	return acos(sum);
}
template<typename T>
T vector1<T>::max_norm()
{
	T max = (*this)[0];
	for (int i = 0; i < n; i++)
	{
		if ((*this)[i] >  max)
		{
			 max = (*this)[i];
		}
	}
	return max;
}
template<typename T>
void vector1<T>::file_bout(ofstream & stream)
{
	if (!stream.is_open())
		throw std::invalid_argument("sorry");
	stream.write(reinterpret_cast<char*>(&n), sizeof(int));
	for (int i = 0; i < n; i++)
		//stream.write((char*)cont, n*sizeof(float));
		stream.write(reinterpret_cast<char*>(&cont[i]), sizeof(float));
}
template<typename T>
vector1<T> vector1<T>::file_bin(ifstream & stream)
{
	int a=0;
	if (!stream.is_open())
		throw std::invalid_argument("sorry");
	stream.read(reinterpret_cast<char*>(&a), sizeof(int));
	//stream >> a;
	vector1<T> s(a);
	for (int i = 0; i < a; i++)
		stream.read(reinterpret_cast<char*>(&s[i]), sizeof(float));
	return s;
}
/*template<typename T>
ofstream & vector1<T>::operator<<(ofstream& stream, vector1<T> s)
{
	if(!stream.is_open())
		throw std::invalid_argument("sorry");
	stream << this->n << endl;
	for (int i = 0; i < this->n; i++)
		stream << (this*)[i]<<" ";
	stream << endl;
}*/
template<typename T>
vector1<T>::~vector1()
{
	//destructor
}
template<typename T>
vector1<T> vector1<T>:: operator*(T a)
{
	vector1<T> arr(n);
	for (int i = 0; i < n; i++)
		arr[i] = (*this)[i] * a;
		return arr;
}
template<typename T>
vector1<T> operator*(vector1<T> & b, T & a)
{
	return a*b;
}

/*
PCA::PCA(matrix<double>& a) :data(a){}

void PCA::centering()
{
	vector<double> m(data.get_width());
	for (int i = 0; i < m.size(); i++) {
		for (int j = 0; j < data.get_height(); j++)
			m[i] += data(i, j);
		m[i] /= data.get_height();
	}
	for (int i = 0; i < m.size(); i++)
		for (int j = 0; j < data.get_height(); j++)
			data(i, j) -= m[i];
}

void PCA::rationing()
{
	vector<double> m(data.get_width());
	for (int i = 0; i < m.size(); i++) {
		for (int j = 0; j < data.get_height(); j++)
			m[i] += pow(data(i, j),2);
		m[i] /= data.get_height();
		m[i] = sqrt(m[i]);
	}
	for (int i = 0; i < m.size(); i++) 
		for (int j = 0; j < data.get_height(); j++)
			data(i, j) /= m[i];
}

*/