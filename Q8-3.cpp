#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
using namespace std;

const int N = 4;    // 階數
const int m = 16;   // 離散點數
const double L = 1; // 區間長度

// 原始函數 f(x) = x^2 * sin(x)
double f(double x) {
    return x * x * sin(x);
}

// S_N(x) 近似三角多項式
double S_N(double x, const vector<double>& a, const vector<double>& b, int N, double L) {
    double result = a[0] / 2.0;
    for (int k = 1; k <= N; ++k) {
        result += a[k] * cos(2 * M_PI * k * x / L);
        result += b[k] * sin(2 * M_PI * k * x / L);
    }
    return result;
}

int main() {
    cout << fixed << setprecision(6);

    // 1. 離散點與函數值
    vector<double> x_j(m), y_j(m);
    for (int j = 0; j < m; ++j) {
        x_j[j] = j * (L / m);        // 均勻點
        y_j[j] = f(x_j[j]);          // f(x_j)
    }

    // 2. 計算傅立葉係數
    vector<double> a(N + 1, 0.0);
    vector<double> b(N + 1, 0.0);

    // a_0 特別處理
    for (int j = 0; j < m; ++j)
        a[0] += y_j[j];
    a[0] *= 2.0 / m;

    for (int k = 1; k <= N; ++k) {
        for (int j = 0; j < m; ++j) {
            double angle = 2 * M_PI * k * x_j[j] / L;
            a[k] += y_j[j] * cos(angle);
            b[k] += y_j[j] * sin(angle);
        }
        a[k] *= 2.0 / m;
        b[k] *= 2.0 / m;
    }

    // 3. 印出 S_4(x) 係數
    cout << "--- S_4(x) Fourier coefficients ---\n";
    cout << "a_0 = " << a[0] << "\n";
    for (int k = 1; k <= N; ++k) {
        cout << "a_" << k << " = " << a[k] << "\n";
        cout << "b_" << k << " = " << b[k] << "\n";
    }

    // 4. S_4(x) 的積分（只需 a_0）
    double integral_s4 = (a[0] / 2.0) * L;
    cout << "\nIntegral of S_4(x) over [0,1] = " << integral_s4 << endl;

    // 5. 誤差 E(S_4)
    double error = 0.0;
    for (int j = 0; j < m; ++j) {
        double approx = S_N(x_j[j], a, b, N, L);
        error += pow(y_j[j] - approx, 2);
    }
    cout << "Squared Error E(S_4) = " << error << endl;

    return 0;
}
