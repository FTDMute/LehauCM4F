#include <iostream>
#include <math.h>
#include <functional> 
#include <string> 
#include <time.h> 
#include <random> 
#include <vector> 
#include <iostream> 
#include "DopFunctions.h"
using namespace std;

void BisectionCalculate(double start, double end) {
	double mid;
	cout << "Метод бисекции  " << endl;
	for (int i = 0; i < 2; i++) {
		if (i == 0) {
			end = (start + end) / 2.;
		}
		else {
			start = 0.;
			end = 2.;
			start = (start + end) / 2.;
		}
		while (abs(start - end) > 1e-5) {
			mid = (start + end) / 2.;
			if (Function(mid) == 0) {
				cout << mid << endl;
				break;
			}
			if (Function(start) * Function(mid) < 0) {
				end = mid;
			}
			else {
				start = mid;
			}
			if (abs(	start - end) < 1e-5) {
				cout << "ответ:  " << (start + end) / 2 << endl;
			}
		}
	}
}

void ChordCalculate(double start, double end) {
	double mid;
	double X1, Xn, tru;
	cout << "Метод хорд  " << endl;
	for (int i = 0; i < 2; i++) {
		if (i == 0) {
			tru = 0;
			Xn = (start + end) / 2.;
			X1 = 10;
			while (tru != 1){
				X1 = Xn - ((Function(Xn) / (Function(Xn) - Function(start))) * (Xn - start));
				if (abs(X1-Xn) < 1e-3) {
					cout << "ответ:  " << X1 << endl;
					tru = 1;
				}
				else {
					Xn = X1;
				}
			}
		}
		else {
			tru = 0;
			Xn = (start + end) / 2.;
			X1 = 10;
			while (tru != 1) {
				X1 = Xn - ((Function(Xn) / (Function(Xn) - Function(end))) * (Xn - end));
				if (abs(X1 - Xn) < 1e-3) {
					cout << "ответ:  " << X1 << endl;
					tru = 1;
				}
				else {
					Xn = X1;
				}
			}
		}
	}
}

void IterationCalculate(double start, double end) {
	double mid, approach;
	cout << "Метод итераций  " << endl;
	for (int i = 0; i < 2; i++) {
		if (i == 0) {
			end = (start + end) / 2.;
			double x0 = (start + end) / 2.,  // Начальное приближение
				x1 = end;
			for (;;)
			{
				x1 = 1/(2*cos(x0));
				if (fabs(x1 - x0) < 1e-5) break;
				x0 = x1;
			}
			cout << "ответ:  " << x1 << endl;
		}
		else {
			start = 0, end = 2;
			start = (start + end) / 2.;
			double x0 = (start + end) / 2,  // Начальное приближение (start + end) / 2.
				x1 = end;
			for (;;)
			{
				x1 = acos(1/(2*x0));
				if (fabs(x1 - x0) < 1e-5) break;
				x0 = x1;
			}
			cout << "ответ:  " << x1 << endl;
		}
	}
}

void CasCalculate(double start, double end) {
	double mid;
	double X1, Xn, tru;
	cout << "Метод касательных  " << endl;
	for (int i = 0; i < 2; i++) {
		if (i == 0) {
			tru = 0;
			end = (start + end) / 2.;
			Xn = (start + end) / 2.;
			X1 = 10;
			while (tru != 1) {
				X1 = Xn - (Function(Xn) / DiffFunction(Xn));
				if (abs(X1 - Xn) < 1e-3) {
					cout << "ответ:  " << X1 << endl;
					tru = 1;
				}
				else {
					Xn = X1;
				}
			}
		}
		else {
			tru = 0;
			start = (start + end) / 2.;
			end = 2.;
			Xn = (start + end) / 2.;
			X1 = 10;
			while (tru != 1) {
				X1 = Xn - (Function(Xn) / DiffFunction(Xn));
				if (abs(X1 - Xn) < 1e-3) {
					cout << "ответ:  " << X1 << endl;
					tru = 1;
				}
				else {
					Xn = X1;
				}
			}
		}
	}
}

void SecCalculate(double start, double end) {
	double Xl;
	double X1, Xn, tru;
	cout << "Метод секущих  " << endl;
	for (int i = 0; i < 2; i++) {
		if (i == 0) {
			tru = 0;
			end = (start + end) / 2.;
			Xn = (start + end) / 2.;
			Xl = Xn / 2.;
			X1 = 10;
			while (tru != 1) {
				X1 = Xn - ((Xn - Xl) * Function(Xn)) / (Function(Xn) - Function(Xl));
				if (abs(X1 - Xn) < 1e-3) {
					cout << "ответ:  " << X1 << endl;
					tru = 1;
				}
				else {
					Xl = Xn;
					Xn = X1;
				}
			}
		}
		else {
			tru = 0;
			start = 0, end = 2;
			start = (start + end) / 2.;
			end = 2.;
			Xn = (start + end) / 2.;
			Xl = Xn / 2.;
			X1 = 10;
			while (tru != 1) {
				X1 = Xn - ((Xn - Xl) * Function(Xn)) / (Function(Xn) - Function(Xl));
				if (abs(X1 - Xn) < 1e-3) {
					cout << "ответ:  " << X1 << endl;
					tru = 1;
				}
				else {
					Xl = Xn;
					Xn = X1;
				}
			}
		}
	}
}

vector<double> sim_sne(double eps, double x_0, double y_0, function<double(double, double)> f_1, function<double(double, double)> f_2)

{
	double x_k = x_0;
	double y_k = y_0;
	double ro = 0.;
	int it = 0;
	do
	{
		if (fabs(f_1(x_k, y_k) - x_k) <= fabs(f_2(x_k, y_k) - y_k))
			ro = fabs(f_2(x_k, y_k) - y_k);
		else
			ro = fabs(f_1(x_k, y_k) - x_k);
		x_k = f_1(x_k, y_k);
		y_k = f_2(x_k, y_k);
		it++;
	} while (ro > eps);
	std::vector<double> roots;
	roots.push_back(x_k);
	roots.push_back(y_k);
	roots.push_back((double)it);
	return roots;
}

vector<double> nm_sne_exact(double eps, double x_0, double y_0,function<double(double, double)> f_1, function<double(double, double)> df_1_x, function<double(double, double)> df_1_y, function<double(double, double)> f_2, function<double(double, double)> df_2_x, function<double(double, double)> df_2_y)
{
	double x_k = x_0;
	double y_k = y_0;
	double ro = 0.;
	int it = 0;
	do
	{
		double temp_x = x_k;
		double temp_y = y_k;
		double det = df_1_x(x_k, y_k) * df_2_y(x_k, y_k) - df_1_y(x_k, y_k) * df_2_x(x_k, y_k);
		x_k = x_k - (df_2_y(x_k, y_k) * f_1(x_k, y_k) - df_1_y(x_k, y_k) * f_2(x_k, y_k)) / det;
		y_k = y_k - (-df_2_x(x_k, y_k) * f_1(x_k, y_k) + df_1_x(x_k, y_k) * f_2(x_k, y_k)) / det;
		if (fabs(temp_x - x_k) <= fabs(temp_y - y_k))
			ro = fabs(temp_y - y_k);
		else
			ro = fabs(temp_x - x_k);
		it++;
	} while (ro > eps);
	vector<double> roots;
	roots.push_back(x_k);
	roots.push_back(y_k);
	roots.push_back((double)it);
	return roots;
}

vector<double> nm_sne_approx(double eps, double x_0, double y_0, function<double(double, double)> f_1, function<double(double, double)> f_2)
{
	double x_k = x_0;
	double y_k = y_0;
	double ro = 0;
	int it = 0;
	double h = 0.001;
	do
	{
		double temp_x = x_k;
		double temp_y = y_k;
		double df_1_x = (f_1(x_k + h, y_k) - f_1(x_k - h, y_k)) / (2. * h);
		double df_1_y = (f_1(x_k, y_k + h) - f_1(x_k, y_k - h)) / (2. * h);
		double df_2_x = (f_2(x_k + h, y_k) - f_2(x_k - h, y_k)) / (2. * h);
		double df_2_y = (f_2(x_k, y_k + h) - f_2(x_k, y_k - h)) / (2. * h);
		double det = df_1_x * df_2_y - df_1_y * df_2_x;
		x_k = x_k - (df_2_y * f_1(x_k, y_k) - df_1_y * f_2(x_k, y_k)) / det;
		y_k = y_k - (-df_2_x * f_1(x_k, y_k) + df_1_x * f_2(x_k, y_k)) / det;
		if (fabs(temp_x - x_k) <= fabs(temp_y - y_k))
			ro = fabs(temp_y - y_k);
		else
			ro = fabs(temp_x - x_k);
		it++;
	} while (ro > eps);
	std::vector<double> roots;
	roots.push_back(x_k);
	roots.push_back(y_k);
	roots.push_back((double)it);
	return roots;
}

int main() {
	setlocale(LC_ALL, "Russian");

	double start, end, startt, endd;
	cout << "Введите начало и конец для (1)" << endl;
	cin >> startt;
	cin >> endd;

	//1

	BisectionCalculate(startt, endd);
	ChordCalculate(startt, endd);
	IterationCalculate(startt, endd);
	CasCalculate(startt, endd);
	SecCalculate(startt, endd);

	cout << "Введите начало и конец для (2) и (3)" << endl;
	cin >> start;
	cin >> end;

	//2

	cout << "Задача 2" << endl;
	vector<double> roots = sim_sne(1e-5, 0, 0, f_1, f_2);
	double x_0 = roots[0];
	double y_0 = roots[1];
	int it = (int)roots[2];
	double F_xx = exp(x_0 * x_0 + y_0 * y_0) * (2 + 4. * x_0 * x_0) - sin(x_0 + y_0);
	double F_xy = exp(x_0 * x_0 + y_0 * y_0) * 4. * x_0 * y_0 - sin(x_0 + y_0);
	double F_yy = exp(x_0 * x_0 + y_0 * y_0) * (2 + 4. * y_0 * y_0) - sin(x_0 + y_0);
	cout << "Решение системы двух нелинейных уравнений: x = " << x_0 << " , y = " << y_0 << endl;
	if (F_xx * F_yy - F_xy * F_xy > 0)
	{
		cout << "точка (" << x_0 << "," << y_0 << ") является точкой экстремума" << endl;
		if (F_xx > 0)
		{
			cout << "точка (" << x_0 << "," << y_0 << ") является точкой минимума" << endl;
			cout << "значение функции z в этой точке: " << exp(x_0 * x_0 + y_0 * y_0) + sin(x_0 + y_0) << endl;
		}
	}
	cout << "------------------------------------------------------" << endl;

	// 3 

	cout << "Задача 3" << endl;
	cout << "Решение уравнения Ньютона с точными производными" << endl << endl;
	roots = nm_sne_exact(1e-5, 0.5, 0.5, F_1, dF_1_x, dF_1_y, F_2, dF_2_x, dF_2_y);
	x_0 = roots[0];
	y_0 = roots[1];
	it = (int)roots[2];
	F_xx = exp(x_0 * x_0 + y_0 * y_0) * (2 + 4. * x_0 * x_0) - sin(x_0 + y_0);
	F_xy = exp(x_0 * x_0 + y_0 * y_0) * 4. * x_0 * y_0 - sin(x_0 + y_0);
	F_yy = exp(x_0 * x_0 + y_0 * y_0) * (2 + 4. * y_0 * y_0) - sin(x_0 + y_0);
	cout << "Решение системы двух нелинейных уравнений: x = " << x_0 << " , y = " << y_0 << ". количество итераций: " << it << endl;
	if (F_xx * F_yy - F_xy * F_xy > 0)
	{
		cout << "точка (" << x_0 << "," << y_0 << ") является точкой экстремума" << endl;
		if (F_xx > 0)
		{
			cout << "точка (" << x_0 << "," << y_0 << ") является точкой минимума" << endl;
			cout << "значение функции z в этой точке: " << exp(x_0 * x_0 + y_0 * y_0) + sin(x_0 + y_0) << endl << endl;
		}
	}
	cout << "Решение уравнения Ньютона с приближенным вычислением производных" << endl << endl;
	roots = nm_sne_approx(1e-5, 0.5, 0.5, F_1, F_2);
	x_0 = roots[0];
	y_0 = roots[1];
	it = (int)roots[2];
	F_xx = exp(x_0 * x_0 + y_0 * y_0) * (2 + 4. * x_0 * x_0) - sin(x_0 + y_0);
	F_xy = exp(x_0 * x_0 + y_0 * y_0) * 4. * x_0 * y_0 - sin(x_0 + y_0);
	F_yy = exp(x_0 * x_0 + y_0 * y_0) * (2 + 4. * y_0 * y_0) - sin(x_0 + y_0);
	cout << "Решение системы двух нелинейных уравнений: x = " << x_0 << " , y = " << y_0 << ". количество итераций: " << it << endl;
	if (F_xx * F_yy - F_xy * F_xy > 0)
	{
		cout << "точка (" << x_0 << "," << y_0 << ") является точкой экстремума" << endl;
		if (F_xx > 0)
		{
			cout << "точка (" << x_0 << "," << y_0 << ") является точкой минимума" << endl;
			cout << "значение функции z в этой точке: " << exp(x_0 * x_0 + y_0 * y_0) + sin(x_0 + y_0) << endl << endl;
		}
	}
	return 0;
}