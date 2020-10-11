#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>

using std::cout;
using std::endl;

#define epc 0.1

double F(const double& x) {
	return (cos(x)* tanh(x));
}

std::vector<std::pair<double, double>> optpas(const double& a, const double& b) {
	std::vector<std::pair<double, double>> val;
	size_t N = 1;
	double del = (b - a) / (N + 1);
	double miny;
	while (del > epc / 2) {
		std::vector<double> vecty;
		del = (b - a) / (N + 1);
		for (size_t k = 1; k <= N; ++k) {
			vecty.push_back(F((b - a) / (N + 1) * k + a));
		}
		size_t ymik = std::min_element(vecty.begin(), vecty.end()) - vecty.begin() + 1;
		miny = (b - a) / (N + 1) * ymik + a;
		val.push_back({ miny, del });
		N++;
	}
	return val;
}

unsigned int fibbon(const size_t& n) {
	if (n < 1)
		return 0;
	unsigned int n1 = 0, n2 = 1, nn = 0;
	for (size_t i = 1; i < n; ++i) {
		nn = n1 + n2;
		n1 = n2;
		n2 = nn;
	}
	return nn;
}


std::vector<double> fibbon(size_t N, double& a, double& b) {
	std::vector<double> val;
	double  x1 = a + (b - a) * fibbon(N) / fibbon(N + 2);
	double  x2 = a + b - x1;
	double  y1 = F(x1);
	double  y2 = F(x2);
	while (N--) {
		if (y1 > y2) {
			a = x1;
			x1 = x2;
			x2 = b - (x1 - a);
			y1 = y2;
			y2 = F(x2);
			val.push_back(x2);
		}
		else {
			b = x2;
			x2 = x1;
			x1 = a + (b - x2);
			y2 = y1;
			y1 = F(x1);
			val.push_back(x1);
		}
	}
	return val;
}

void PrintValues(const std::vector<double>& val) {
	cout << std::string(40, '-') << endl;
	cout << std::setw(18) << std::left << "| Количество точек N " << "|";
	cout << std::setw(19) << std::left << " Значение x в точке |" << endl;
	cout << std::string(40, '-') << endl;

	for (size_t i = 0; i < val.size(); i++) {
		cout << "|";
		cout << std::setw(10) << std::right << i + 1;
		cout << std::setw(10) << "|";
		cout << std::setw(15) << std::right << std::setprecision(11) << val[i];
		cout << std::setw(2) << "|" << endl;
	}
	cout << std::string(40, '-') << endl;
	cout << val[7] << endl;
}

void PrintValues(const std::vector<std::pair<double, double>>& val) {
	cout << std::string(60, '-') << endl;
	cout << std::setw(20) << std::left << "| Количество точек N " << "|";
	cout << std::setw(20) << std::left << " Значение x в минимуме |";
	cout << std::setw(10) << std::left << " Погрешность |" << endl;
	cout << std::string(60, '-') << endl;

	int i = 0;
	double an1, an2;
	for (auto const& val : val)
	{
		cout << "|";
		cout << std::setw(10) << std::right << i + 1;
		cout << std::setw(10) << "|";
		cout << std::setw(15) << std::setprecision(11) << val.first << std::setw(6) << "|";
		cout << std::setw(10) << std::setprecision(3) << val.second;
		cout << std::setw(5) << "|" << endl;
		i++;
		an1 = val.first;
		an2 = val.second;
	}
	cout << std::string(60, '-') << endl;
	cout << an1 << " +- " << an2 << endl;
}
int main() {
	setlocale(LC_ALL, "Russian");

	size_t N = 30;
	double a = 7.0, b = 11.0;

	cout << "Функция cos(x)* tanh(x) для интервала [7, 11]" << endl << endl;
	cout << "Метод оптимального пассивного поиска" << endl;
	cout << "Для погрешности 0,1" << endl;
	PrintValues(optpas(a, b));
	cout << endl;
	cout << "Метод Фибоначчи" << endl;
	cout << "При N=30 для точности 0,1 достаточно 8 итераций" << endl;
	PrintValues(fibbon(N, a, b));
	return 0;
}
