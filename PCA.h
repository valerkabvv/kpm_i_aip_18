#pragma once
#define eps  1e-8 
#include"Header.h"
#include"def.cpp"
#include<vector>


template<typename T>
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
}


class PCA {
private:
	matrix<double> data;
public:
	PCA(matrix<double> &a) :data(a) {};
	matrix<double>  centering(matrix<double> new_data)
	{
		//matrix<double> new_data(data);
		vector<double> m(new_data.get_width());
		for (int i = 0; i < m.size(); i++) {
			for (int j = 0; j < new_data.get_height(); j++)
				m[i] += new_data(j, i);
			m[i] /= new_data.get_height();
		}
		for (int i = 0; i < m.size(); i++)
			for (int j = 0; j < new_data.get_height(); j++)
				new_data(j, i) -= m[i];
		return new_data;
	};
	matrix<double> rationing(matrix<double> new_data)
	{
		//matrix<double> new_data(data);
		vector<double> m(new_data.get_width());
		for (int i = 0; i < m.size(); i++) {
			for (int j = 0; j < new_data.get_height(); j++)
				m[i] += pow(new_data(j, i), 2);
			m[i] /= (new_data.get_height()-1);
			m[i] = sqrt(m[i]);
		}
		for (int i = 0; i < m.size(); i++)
			for (int j = 0; j < new_data.get_height(); j++)
				new_data(j, i) /= m[i];
		return new_data;
	};
	vector<matrix<double>> NIPALS(int dim)
	{
		matrix<double> D(this->rationing(this->centering(data)));
		matrix<double> E(D);
		matrix<double> T(E.get_height(), dim);
		matrix<double> P(dim, E.get_width());
		matrix<double> d(E.get_height(), 1);
		matrix<double> p(E.get_width(), 1);
		for (int i = 0; i < dim; i++) {
			matrix<double> t(E.col(i));
			do {
				 p = ((t.transpose()*E)*(1 / (t.transpose()*t)(0, 0)));
				 p = p.transpose();
				p = p*(1 / p.norm());
				matrix<double> t_old(t);
				t = (E*p)* (1/(p.transpose()*p)(0,0));
				d = t_old - t;
			} while (d.norm() > eps);
			E = E - t*p.transpose();
			P.append_row(p.transpose(), i);
			T.append_col(t, i);
		}
		vector<matrix<double>> res;
		res.push_back(E);
		res.push_back(P*(-1));
		res.push_back(T*(-1));
		return res;
	}
	vector<vector<double>> lever_and_deviation(int dim) {
		vector<matrix<double>> a = NIPALS(dim);
		matrix<double> eigenval = (a[2].transpose()*a[2]).inverse_matrix();
		vector<double> lever;
		matrix<double> help(a[2].transpose());
		for (int i = 0; i < a[2].get_height(); i++) {
			matrix<double> b(help.col(i));
			lever.push_back((b.transpose()*eigenval*b)(0,0));
		}
		vector<double> dev;
		for (int i = 0; i < a[0].get_height(); i++)
			dev.push_back(a[0][i].scalar_mul(a[0][i]));
		vector<vector<double>> res;
		res.push_back(lever);
		res.push_back(dev);
		return res;
	}
	vector<double> TRVC_ERVC(int dim) {
		vector<matrix<double>> a = NIPALS(dim);
		vector<double> res(2);
		for (int i = 0; i < a[0].get_height(); i++)
			for (int j = 0; j < a[0].get_width(); j++)
				res[0] += pow(a[0](i, j),2);
		for (int i = 0; i < data.get_height(); i++)
			for (int j = 0; j < data.get_width(); j++)
				res[1] += pow(data(i, j),2);
		res[1] = 1 - (res[0] / res[1]);
		res[0] /= a[0].get_height()*a[0].get_width();
		return res;
	}
};