#include "Header.h"
#include"def.cpp"
#include"PCA.h"
#include <fstream>
using namespace std;



template<typename T>
ofstream & operator<<(ofstream& stream, vector1<T>& s)
{
	if (!stream.is_open())
		throw std::invalid_argument("sorry");
	stream << s.get_size() << endl;
	for (int i = 0; i < s.get_size(); i++)
		stream << s[i] << " ";
	stream << endl;
}



template<typename T>
vector1<T> operator>>(ifstream& stream, vector1<T>& s)
{
	int a;
	if (!stream.is_open())
		throw std::invalid_argument("sorry");
	stream >> a;
	//s = vector1<T>(a);
	vector1<T> b(a);
	for (int i = 0; i < a; i++)
		stream >> b[i] ;
	return b;
}



template<typename T>
ofstream & operator<<(ofstream& stream, matrix<T> s)
{
	if (!stream.is_open())
		throw std::invalid_argument("sorry");
	stream << s.get_height()<<" "<<s.get_width() << endl;
	for (int i = 0; i < s.get_height(); i++) {
		for (int j = 0; j < s.get_width(); j++)
			stream << s(i,j) << " ";
		stream << endl;
	}
}



template<typename T>
matrix<T> operator>>(ifstream& stream, matrix<T>& s)
{
	int a, b;
	if (!stream.is_open())
		throw std::invalid_argument("sorry");
	stream >> a>>b;
	matrix<T> c(a,b);
	for (int i = 0; i < a; i++)
		for(int j=0;j<b;j++)
		stream >> c(i,j);
	return c;
}



/*template<typename T>
ostream& operator<<(ostream & stream, matrix<T> s)
{
	for (int i = 0; i < s.get_height(); i++)
	{
		stream << "|";
		for (int j = 0; j < s.get_width(); j++)
			stream << s(i, j) << " ";
		stream << "|";
		stream << endl;
	}
	return stream;
}*/





template<typename T>
ostream& operator<<(ostream & stream, vector<vector<T>> s)
{
	for (int i = 0; i < s.size(); i++)
	{
		for (int j = 0; j < s[0].size(); j++)
			stream << s[i][j] << " ";
		stream << endl;
	}
	return stream;
}




template<typename T>
ostream& operator<<(ostream & stream, vector1<T> s)
{
	stream << "[";
	for (int i = 0; i < s.get_size(); i++)
		stream << s[i] << " ";
	stream << "]";
	stream << endl;
	return stream;
}



void test_vect()
{
	int m, n;
	cout << "Razmer Vect1" << endl;
	cin >> m;
	cout << "Razmer Vect2" << endl;
	cin >> n;
	vector1<double> m1(m);
	vector1<double> n1(n);
	for (int i = 0; i < m; i++) 
		cin >> m1[i];
	for (int i = 0; i < n; i++)
		cin >> n1[i];
	cout << endl;
	cout << m1;
	cout << endl;
	cout << n1;
	cout << " vnesh proizv" << endl << m1.out_mul(n1) << endl << endl;
	cout << " angle" << endl << m1.angle(n1) << endl << endl;
	cout << " skalar" << endl << m1.scalar_mul(n1) << endl << endl;
	cout << "frab_norm" << endl << m1.norm() << endl << endl;
	cout << "max_norm" << endl << m1.max_norm() << endl << endl;
}
void test_matr()
{
	cout << "Razmer matr1" << endl;
	int m, n,m1,n1;
	cin >> m >> n;
	cout << "Razmer matr2" << endl;
	cin>> m1 >> n1;
	matrix<double> m2(m, n);
	matrix<double> n2(m1,n1);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
		cin >> m2(i,j);
	for (int i = 0; i < m1; i++)
		for (int j = 0; j < n1; j++)
			cin >> n2(i, j);
	cout << m2;
	cout << endl;
	cout << n2;
	cout << " norm" << endl << m2.norm() << endl << endl;
	cout << " rk" << endl << m2.rk() << endl << endl;
	cout << " trans" << endl << m2.transpose() << endl << endl;
	try
	{
		cout << " inverse matr" << endl << m2.inverse_matrix() << endl << m2.inverse_matrix()*m2 << endl << endl;
		cout << " det" << endl << m2.det() << endl << endl;
		cout << " hamadar" << endl << m2.Hadamar_mul(n2) << endl << endl;
	}
	catch (exception& some) {
		cerr << some.what();
	}
}
void test_text1() {
	int n,m;
	//vector1<float> a();
	ofstream strm("data.txt");
	cin >> n;
	for (int i = 0; i < n; i++) {
		cout << "size then vctr" << endl;
		cin >> m;
		vector1<float> a(m);
		for (int j = 0; j < m; j++)
			cin >> a[j];
		strm << a;
	}
	strm.close();

	ifstream str("data.txt");
	vector<vector1<float>> mass;
	for (int i = 0; i < n; i++) {
		vector1<float> a(2);
		vector1<float> b (str >> a);
		mass.push_back(b);
	}
	system("pause");
}

void test_text2() {
	int n, m,k;
	//vector1<float> a();
	ofstream strm("data_matr.txt");
	cin >> n;
	for (int i = 0; i < n; i++) {
		cout << "size then mtrx" << endl;
		cin >> m>>k;
		matrix<float> a(m,k);
		for (int j = 0; j < m; j++)
			for(int g = 0;g<k;g++)
               cin >> a(j,g);
		strm << a;
	}
	strm.close();

	ifstream str("data_matr.txt");
	vector<matrix<float>> mass;
	for (int i = 0; i < n; i++) {
		matrix<float> a(2,2);
		matrix<float> b(str >> a);
		mass.push_back(b);
	}
	system("pause");
}

void test_text1bin() {
	int n, m;
	//vector1<float> a();
	ofstream strm("data1.dat",ios::binary);
	cin >> n;
	for (int i = 0; i < n; i++) {
		cout << "size then vctr" << endl;
		cin >> m;
		vector1<float> a(m);
		for (int j = 0; j < m; j++)
			cin >> a[j];
		a.file_bout(strm);
	}
	strm.close();

	ifstream str("data1.dat",ios::binary);
	vector<vector1<float>> mass;
	for (int i = 0; i < n; i++) {
		vector1<float> a(2);
		vector1<float> b =a.file_bin(str);
		mass.push_back(b);
	}
	system("pause");
}

void test_text2bin() {
	int n, m, k;
	//vector1<float> a();
	ofstream strm("data_matr1.dat", ios::binary);
	cin >> n;
	for (int i = 0; i < n; i++) {
		cout << "size then mtrx" << endl;
		cin >> m >> k;
		matrix<float> a(m, k);
		for (int j = 0; j < m; j++)
			for (int g = 0; g<k; g++)
				cin >> a(j, g);
		a.file_bout(strm);
	}
	strm.close();

	ifstream str("data_matr1.dat", ios::binary);
	vector<matrix<float>> mass;
	for (int i = 0; i < n; i++) {
		matrix<float> a(2, 2);
		matrix<float> b(a.file_bin(str));
		mass.push_back(b);
	}
	system("pause");
}

int main()
{
	ifstream stream("task.txt");
	matrix<double> matr(stream);
	PCA test(matr);
	vector<matrix<double>> a=test.NIPALS(12);
	for (int i = 0; i < 3; i++)
		cout << a[i]<<endl;
	vector<double> b = test.TRVC_ERVC(12);
	for (int i = 0; i < 2; i++)
		cout <<b[i] << endl;
	vector<vector<double>> c = test.lever_and_deviation(12);
	for (int i = 0; i < 2; i++)
		for(int j=0;j<c[i].size();j++)
		cout << c[i][j] << endl;
	system("pause");


/*
test_text2bin();
	int m;
	cout << "size then vctr" << endl;
	cin >> m;
	vector1<float> a(m);
	for (int j = 0; j < m; j++)
		cin >> a[j];
	vector1<float> b(a);

	*/


/*	ofstream stream("data1.txt",ios::binary);
	if (!stream.is_open())
		throw std::invalid_argument("sorry");
	int b;
	float a = 10000000000;
	for (int i = 0; i < 100; i++)
	{
		stream << a;
		//matrix<float> tst;
	//	stream >> tst;
	}
	/*for (int i = 0; i < 1; i++)
	{
		cin >> a;
		vector1<float> tst(a);
		for (int j = 0; j < a; j++)
			cin >> tst[j];
		stream << tst;
	}*/
	//test_matr();
	//test_vect();
	
	system("pause");
}