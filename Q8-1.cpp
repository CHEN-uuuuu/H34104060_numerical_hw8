#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

// 高斯消去法求解 3x3 線性方程組
void solve3x3(double A[3][4], double result[3]) {
    // 前向消去
    for (int i = 0; i < 3; ++i) {
        // 將主對角線歸一
        double pivot = A[i][i];
        for (int j = 0; j < 4; ++j)
            A[i][j] /= pivot;

        // 將其他行的第 i 列清為 0
        for (int k = 0; k < 3; ++k) {
            if (k == i) continue;
            double factor = A[k][i];
            for (int j = 0; j < 4; ++j)
                A[k][j] -= factor * A[i][j];
        }
    }

    for (int i = 0; i < 3; ++i)
        result[i] = A[i][3];
}

int main() {
    vector<double> x = {4.0, 4.2, 4.5, 4.7, 5.1, 5.5, 5.9, 6.3};
    vector<double> y = {102.6, 113.2, 130.1, 142.1, 167.5, 195.1, 224.9, 256.8};
    int n = x.size();

    cout << fixed << setprecision(6);

    // ---------- Part (a): y = ax^2 + bx + c ----------
    double Sx = 0, Sx2 = 0, Sx3 = 0, Sx4 = 0;
    double Sy = 0, Sxy = 0, Sx2y = 0;
    for (int i = 0; i < n; ++i) {
        double xi = x[i], yi = y[i];
        double xi2 = xi * xi;
        double xi3 = xi2 * xi;
        double xi4 = xi3 * xi;

        Sx += xi;
        Sx2 += xi2;
        Sx3 += xi3;
        Sx4 += xi4;
        Sy += yi;
        Sxy += xi * yi;
        Sx2y += xi2 * yi;
    }

    double A[3][4] = {
        {Sx4, Sx3, Sx2, Sx2y},
        {Sx3, Sx2, Sx,  Sxy},
        {Sx2, Sx,  (double)n, Sy}
    };

    double coef[3]; // a, b, c
    solve3x3(A, coef);
    double a = coef[0], b = coef[1], c = coef[2];

    cout << "--- Part a: y = ax^2 + bx + c ---\n";
    cout << "Coefficients: a = " << a << ", b = " << b << ", c = " << c << endl;

    double error_a = 0.0;
    for (int i = 0; i < n; ++i) {
        double y_hat = a * x[i] * x[i] + b * x[i] + c;
        error_a += pow(y[i] - y_hat, 2);
    }
    cout << "Sum of squared error: " << error_a << "\n\n";

    // ---------- Part (b): y = b * e^(a * x) ----------
    vector<double> ln_y(n);
    for (int i = 0; i < n; ++i)
        ln_y[i] = log(y[i]);

    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
    for (int i = 0; i < n; ++i) {
        sumX += x[i];
        sumY += ln_y[i];
        sumXY += x[i] * ln_y[i];
        sumX2 += x[i] * x[i];
    }

    double denom = n * sumX2 - sumX * sumX;
    double a_b = (n * sumXY - sumX * sumY) / denom;
    double ln_b_b = (sumY * sumX2 - sumX * sumXY) / denom;
    double b_b = exp(ln_b_b);

    cout << "--- Part b: y = b * e^(a * x) ---\n";
    cout << "Parameters: b = " << b_b << ", a = " << a_b << endl;

    double error_b = 0.0;
    for (int i = 0; i < n; ++i) {
        double y_hat = b_b * exp(a_b * x[i]);
        error_b += pow(y[i] - y_hat, 2);
    }
    cout << "Sum of squared error: " << error_b << "\n\n";

    // ---------- Part (c): y = b * x^a ----------
    vector<double> ln_x(n);
    for (int i = 0; i < n; ++i)
        ln_x[i] = log(x[i]);

    sumX = sumY = sumXY = sumX2 = 0;
    for (int i = 0; i < n; ++i) {
        sumX += ln_x[i];
        sumY += ln_y[i];
        sumXY += ln_x[i] * ln_y[i];
        sumX2 += ln_x[i] * ln_x[i];
    }

    denom = n * sumX2 - sumX * sumX;
    double a_c = (n * sumXY - sumX * sumY) / denom;
    double ln_b_c = (sumY * sumX2 - sumX * sumXY) / denom;
    double b_c = exp(ln_b_c);

    cout << "--- Part c: y = b * x^a ---\n";
    cout << "Parameters: b = " << b_c << ", a = " << a_c << endl;

    double error_c = 0.0;
    for (int i = 0; i < n; ++i) {
        double y_hat = b_c * pow(x[i], a_c);
        error_c += pow(y[i] - y_hat, 2);
    }
    cout << "Sum of squared error: " << error_c << endl;

    return 0;
}
