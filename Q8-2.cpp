#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

const double a = -1.0;
const double b = 1.0;
const int N = 1000; // 分段數（必須為偶數）

// 被近似的函數 f(x)
double f(double x) {
    return 0.5 * cos(x) - 0.5 * sin(x) + 0.25 * x;
}

// 基底函數 φ_i(x)
double phi(int i, double x) {
    if (i == 0) return 1.0;
    if (i == 1) return x;
    if (i == 2) return x * x;
    return 0.0;
}

// ✅ 改為 template，可接受 lambda、functor 或普通函數
template <typename Func>
double simpson_integral(Func func, double a, double b, int n) {
    double h = (b - a) / n;
    double sum = func(a) + func(b);
    for (int i = 1; i < n; i += 2)
        sum += 4 * func(a + i * h);
    for (int i = 2; i < n; i += 2)
        sum += 2 * func(a + i * h);
    return sum * h / 3.0;
}

// 用來傳遞 λx → φ_i(x) * φ_j(x)
struct BasisProduct {
    int i, j;
    double operator()(double x) const {
        return phi(i, x) * phi(j, x);
    }
};

// 傳遞 λx → f(x) * φ_i(x)
struct FTimesBasis {
    int i;
    double operator()(double x) const {
        return f(x) * phi(i, x);
    }
};

int main() {
    const int degree = 2;
    double G[3][3];  // Gram 矩陣
    double b_vec[3]; // 右側向量

    // 建立 Gram 矩陣 G[i][j] = ∫ φ_i(x) φ_j(x) dx
    for (int i = 0; i <= degree; ++i) {
        for (int j = 0; j <= degree; ++j) {
            BasisProduct bp{i, j};
            G[i][j] = simpson_integral(bp, a, b, N);
        }
        FTimesBasis fb{i};
        b_vec[i] = simpson_integral(fb, a, b, N);
    }

    // 高斯消去法解線性系統 G * x = b
    double A[3][4]; // 擴增矩陣
    for (int i = 0; i <= degree; ++i) {
        for (int j = 0; j <= degree; ++j)
            A[i][j] = G[i][j];
        A[i][3] = b_vec[i];
    }

    // 高斯消去法
    for (int i = 0; i < 3; ++i) {
        double pivot = A[i][i];
        for (int j = 0; j < 4; ++j)
            A[i][j] /= pivot;
        for (int k = 0; k < 3; ++k) {
            if (k == i) continue;
            double factor = A[k][i];
            for (int j = 0; j < 4; ++j)
                A[k][j] -= factor * A[i][j];
        }
    }

    double coef[3];
    for (int i = 0; i < 3; ++i)
        coef[i] = A[i][3];

    cout << fixed << setprecision(6);
    cout << "Degree 2 least squares polynomial approximation:\n";
    cout << "P2(x) = " << coef[0]
         << " + " << coef[1] << " * x"
         << " + " << coef[2] << " * x^2" << endl;

    return 0;
}
