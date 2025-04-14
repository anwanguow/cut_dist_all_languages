#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <limits>
#include <algorithm>
#include <utility>

using namespace std;

using Vector = vector<double>;
using Matrix = vector<Vector>;

size_t numRows(const Matrix& M) {
    return M.size();
}

size_t numCols(const Matrix& M) {
    return M.empty() ? 0 : M[0].size();
}

Matrix subtractMatrix(const Matrix& A, const Matrix& B) {
    size_t r = numRows(A), c = numCols(A);
    Matrix C(r, Vector(c, 0.0));
    for (size_t i = 0; i < r; i++){
        for (size_t j = 0; j < c; j++){
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return C;
}

Matrix divideMatrix(const Matrix& A, double scalar) {
    size_t r = numRows(A), c = numCols(A);
    Matrix C(r, Vector(c, 0.0));
    for (size_t i = 0; i < r; i++){
        for (size_t j = 0; j < c; j++){
            C[i][j] = A[i][j] / scalar;
        }
    }
    return C;
}

double matrixSum(const Matrix& A) {
    double s = 0.0;
    for (const auto &row : A)
        for (double val : row)
            s += val;
    return s;
}

Vector colSum(const Matrix& A) {
    size_t r = numRows(A), c = numCols(A);
    Vector sum(c, 0.0);
    for (size_t j = 0; j < c; j++){
        for (size_t i = 0; i < r; i++){
            sum[j] += A[i][j];
        }
    }
    return sum;
}

Vector rowSum(const Matrix& A) {
    size_t r = numRows(A), c = numCols(A);
    Vector sum(r, 0.0);
    for (size_t i = 0; i < r; i++){
        for (size_t j = 0; j < c; j++){
            sum[i] += A[i][j];
        }
    }
    return sum;
}

void appendColumn(Matrix& A, const Vector& col) {
    size_t r = numRows(A);
    for (size_t i = 0; i < r; i++){
        A[i].push_back(col[i]);
    }
}

void appendRow(Matrix& A, const Vector& row) {
    A.push_back(row);
}

Matrix transpose(const Matrix& A) {
    size_t r = numRows(A), c = numCols(A);
    Matrix T(c, Vector(r, 0.0));
    for (size_t i = 0; i < r; i++){
        for (size_t j = 0; j < c; j++){
            T[j][i] = A[i][j];
        }
    }
    return T;
}

Matrix matMultiply(const Matrix& A, const Matrix& B) {
    size_t r = numRows(A), c = numCols(B), inner = numCols(A);
    Matrix C(r, Vector(c, 0.0));
    for (size_t i = 0; i < r; i++){
        for (size_t j = 0; j < c; j++){
            for (size_t k = 0; k < inner; k++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

double normVector(const Vector& v) {
    double s = 0.0;
    for (double x : v)
        s += x * x;
    return sqrt(s);
}

void normalizeVector(Vector& v) {
    double n = normVector(v);
    if (n > 1e-8) {
        for (auto &x : v)
            x /= n;
    }
}

void normalizeColumns(Matrix& X) {
    size_t rows = numRows(X), cols = numCols(X);
    for (size_t j = 0; j < cols; j++){
        double sum_sq = 0.0;
        for (size_t i = 0; i < rows; i++){
            sum_sq += X[i][j] * X[i][j];
        }
        double nrm = sqrt(sum_sq);
        if (nrm > 1e-8) {
            for (size_t i = 0; i < rows; i++){
                X[i][j] /= nrm;
            }
        }
    }
}

Matrix randomNormalMatrix(size_t rows, size_t cols) {
    Matrix M(rows, Vector(cols, 0.0));
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> d(0, 1);
    for (size_t i = 0; i < rows; i++){
        for (size_t j = 0; j < cols; j++){
            M[i][j] = d(gen);
        }
    }
    return M;
}

Vector colSquareSum(const Matrix& M) {
    size_t rows = numRows(M), cols = numCols(M);
    Vector res(cols, 0.0);
    for (size_t j = 0; j < cols; j++){
        for (size_t i = 0; i < rows; i++){
            res[j] += M[i][j] * M[i][j];
        }
    }
    return res;
}

Vector sqrtVector(const Vector& v) {
    Vector res(v.size());
    for (size_t i = 0; i < v.size(); i++){
        res[i] = sqrt(v[i]);
    }
    return res;
}

Matrix broadcastDivide(const Matrix& M, const Vector& v) {
    size_t rows = numRows(M), cols = numCols(M);
    Matrix res = M;
    for (size_t j = 0; j < cols; j++){
        for (size_t i = 0; i < rows; i++){
            res[i][j] /= v[j];
        }
    }
    return res;
}

Vector colElementwiseProductSum(const Matrix& A, const Matrix& B) {
    size_t rows = numRows(A), cols = numCols(A);
    Vector res(cols, 0.0);
    for (size_t j = 0; j < cols; j++){
        for (size_t i = 0; i < rows; i++){
            res[j] += A[i][j] * B[i][j];
        }
    }
    return res;
}

Vector elementwiseMultiply(const Vector& a, const Vector& b) {
    Vector res(a.size(), 0.0);
    for (size_t i = 0; i < a.size(); i++){
        res[i] = a[i] * b[i];
    }
    return res;
}

Matrix elementwiseMultiply(const Matrix& A, const Matrix& B) {
    size_t rows = numRows(A), cols = numCols(A);
    Matrix res(rows, Vector(cols, 0.0));
    for (size_t i = 0; i < rows; i++){
        for (size_t j = 0; j < cols; j++){
            res[i][j] = A[i][j] * B[i][j];
        }
    }
    return res;
}

Matrix broadcastMultiply(const Matrix& M, const Vector& v) {
    size_t rows = numRows(M), cols = numCols(M);
    Matrix res = M;
    for (size_t j = 0; j < cols; j++){
        for (size_t i = 0; i < rows; i++){
            res[i][j] *= v[j];
        }
    }
    return res;
}

Matrix elementwiseSubtract(const Matrix& A, const Matrix& B) {
    size_t rows = numRows(A), cols = numCols(A);
    Matrix res(rows, Vector(cols, 0.0));
    for (size_t i = 0; i < rows; i++){
        for (size_t j = 0; j < cols; j++){
            res[i][j] = A[i][j] - B[i][j];
        }
    }
    return res;
}

Matrix add(const Matrix& A, const Matrix& B) {
    size_t rows = numRows(A), cols = numCols(A);
    Matrix res(rows, Vector(cols, 0.0));
    for (size_t i = 0; i < rows; i++){
        for (size_t j = 0; j < cols; j++){
            res[i][j] = A[i][j] + B[i][j];
        }
    }
    return res;
}

double frobeniusNorm(const Matrix& M) {
    double sum = 0.0;
    for (const auto &row : M)
        for (double val : row)
            sum += val * val;
    return sqrt(sum);
}

double matrixElementwiseSum(const Matrix& M) {
    double s = 0.0;
    for (const auto &row : M)
        for (double val : row)
            s += val;
    return s;
}

double clamp_val(double x, double lower, double upper) {
    return max(lower, min(x, upper));
}

pair<double, Matrix> cut_norm_quad(const Matrix& V, const Matrix& A) {
    size_t n = numRows(A);
    size_t p = numRows(V);
    size_t n2 = numCols(V);
    size_t half = n2 / 2;
    
    Matrix U(p, Vector(half, 0.0));
    Matrix V2(p, Vector(half, 0.0));
    for (size_t i = 0; i < p; i++){
        for (size_t j = 0; j < half; j++){
            U[i][j] = V[i][j];
            V2[i][j] = V[i][j + half];
        }
    }
    
    Matrix A_T = transpose(A);
    Matrix part1 = matMultiply(V2, A_T);
    Matrix part2 = matMultiply(U, A);
    
    Matrix g(p, Vector(2 * n, 0.0));
    for (size_t i = 0; i < p; i++){
        for (size_t j = 0; j < n; j++){
            g[i][j] = 2 * part1[i][j];
            g[i][j + n] = 2 * part2[i][j];
        }
    }
    
    double f = 0.0;
    for (size_t i = 0; i < p; i++){
        for (size_t j = 0; j < n; j++){
            f += g[i][j] * U[i][j] + g[i][j + n] * V2[i][j];
        }
    }
    f /= 2.0;
    
    return {f, g};
}

pair<Matrix, Matrix> optimize(Matrix x, const Matrix& A, double xtol=1e-8, double ftol=1e-10, double gtol=1e-8, double rho=1e-4,
             double eta=0.1, double gamma=0.85, double tau=1e-3,
             int nt=5, int mxitr=600) {
    Matrix crit(mxitr, Vector(3, 1.0));
    size_t n = numRows(x);
    size_t p = numCols(x);
    Vector nrmx = colSquareSum(x);
    if (normVector(nrmx) > 1e-8) {
        Vector sqrt_nrmx = sqrtVector(nrmx);
        x = broadcastDivide(x, sqrt_nrmx);
    }
    double f;
    Matrix g;
    tie(f, g) = cut_norm_quad(x, A);
    Vector xtg = colElementwiseProductSum(x, g);
    Vector gg = colSquareSum(g);
    Vector xx = colSquareSum(x);
    Vector xxgg = elementwiseMultiply(xx, gg);
    Matrix dtX = broadcastMultiply(x, xtg);
    dtX = elementwiseSubtract(dtX, g);
    double nrmG = frobeniusNorm(dtX);
    double Q = 1.0;
    double Cval = f;
    double tau_orig = tau;
    tau = tau_orig;
    for (size_t itr = 0; itr < static_cast<size_t>(mxitr); itr++){
        Matrix xp = x;
        double fp = f;
        Matrix gp = g;
        Matrix dtXP = dtX;
        int nls = 1;
        double deriv = rho * nrmG * nrmG;
        while (true) {
            double tau2 = tau / 2.0;
            Vector beta(p, 0.0);
            for (size_t j = 0; j < p; j++){
                beta[j] = 1.0 + tau2 * tau2 * (- xtg[j] * xtg[j] + xxgg[j]);
            }
            Vector a1(p, 0.0);
            for (size_t j = 0; j < p; j++){
                double temp = 1.0 + tau2 * xtg[j];
                a1[j] = ((temp * temp) - tau2 * tau2 * xxgg[j]) / beta[j];
            }
            Vector a2(p, 0.0);
            for (size_t j = 0; j < p; j++){
                a2[j] = -tau * xx[j] / beta[j];
            }
            Matrix x_new = add(broadcastMultiply(xp, a1), broadcastMultiply(gp, a2));
            x = x_new;
            tie(f, g) = cut_norm_quad(x, A);
            if (f <= Cval - tau * deriv || nls >= 5) {
                break;
            }
            tau = eta * tau;
            nls = nls + 1;
        }
        xtg = colElementwiseProductSum(x, g);
        gg = colSquareSum(g);
        xx = colSquareSum(x);
        xxgg = elementwiseMultiply(xx, gg);
        dtX = broadcastMultiply(x, xtg);
        dtX = elementwiseSubtract(dtX, g);
        nrmG = frobeniusNorm(dtX);
        Matrix s = elementwiseSubtract(x, xp);
        double XDiff = frobeniusNorm(s) / sqrt(n);
        double FDiff = fabs(fp - f) / (fabs(fp) + 1.0);
        crit[itr][0] = nrmG;
        crit[itr][1] = XDiff;
        crit[itr][2] = FDiff;
        int start = max(0, static_cast<int>(itr) - min(nt, static_cast<int>(itr)));
        Vector mcrit(3, 0.0);
        int count = 0;
        for (int k = start; k <= static_cast<int>(itr); k++){
            for (size_t j = 0; j < 3; j++){
                mcrit[j] += crit[k][j];
            }
            count++;
        }
        for (size_t j = 0; j < 3; j++){
            mcrit[j] /= count;
        }
        if ((XDiff < xtol && FDiff < ftol) || nrmG < gtol || (mcrit[1] < 10 * xtol && mcrit[2] < 10 * ftol)) {
            break;
        }
        Matrix y = elementwiseSubtract(dtX, dtXP);
        Matrix sy_mat = elementwiseMultiply(s, y);
        double sy = fabs(matrixSum(sy_mat));
        tau = tau_orig;
        if (sy > 0) {
            double sum_ss = matrixElementwiseSum(elementwiseMultiply(s, s));
            double sum_yy = matrixElementwiseSum(elementwiseMultiply(y, y));
            if (itr % 2 == 0) {
                tau = sum_ss / sy;
            } else {
                tau = sy / sum_yy;
            }
            tau = clamp_val(tau, 1e-20, 1e20);
        }
        double Qp = Q;
        Q = gamma * Qp + 1.0;
        Cval = (gamma * Qp * Cval + f) / Q;
    }
    return {x, g};
}

double cut_norm(const Matrix& A_input) {
    Matrix A = A_input;
    size_t n1 = numRows(A);
    Vector A_col = colSum(A);
    Vector A_row = rowSum(A);
    double A_tot = matrixSum(A);
    Vector negA_row(n1, 0.0);
    for (size_t i = 0; i < n1; i++){
        negA_row[i] = -A_row[i];
    }
    appendColumn(A, negA_row);
    Vector newRow;
    for (double val : A_col)
        newRow.push_back(-val);
    newRow.push_back(A_tot);
    appendRow(A, newRow);
    size_t p_val = max(min(static_cast<size_t>(round(sqrt(2 * n1) / 2.0)), static_cast<size_t>(100)), static_cast<size_t>(1));
    size_t n2 = 2 * n1 + 2;
    Matrix x0 = randomNormalMatrix(p_val, n2);
    normalizeColumns(x0);
    Matrix x, dummyG;
    tie(x, dummyG) = optimize(x0, A);
    size_t half = n2 / 2;
    Matrix U(p_val, Vector(half, 0.0));
    Matrix V(p_val, Vector(half, 0.0));
    for (size_t i = 0; i < p_val; i++){
        for (size_t j = 0; j < half; j++){
            U[i][j] = x[i][j];
            V[i][j] = x[i][j + half];
        }
    }
    Matrix U_trans = transpose(U);
    Matrix prod = matMultiply(U_trans, V);
    double sdp = 0.0;
    size_t r = numRows(A), c = numCols(A);
    for (size_t i = 0; i < r && i < prod.size(); i++){
        for (size_t j = 0; j < c && j < prod[0].size(); j++){
            sdp += A[i][j] * prod[i][j];
        }
    }
    double sdp_result = fabs(sdp) / 4.0;
    return sdp_result;
}

double cut_distance(const Matrix& G, const Matrix& H) {
    size_t n = numRows(G);
    Matrix diff = subtractMatrix(G, H);
    Matrix A = divideMatrix(diff, double(n * n));
    return cut_norm(A);
}

Matrix erdosRenyiGraph(int n, double p) {
    Matrix A(n, Vector(n, 0.0));
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);
    for (int i = 0; i < n; i++){
        for (int j = i + 1; j < n; j++){
            if (dis(gen) < p) {
                A[i][j] = 1.0;
                A[j][i] = 1.0;
            }
        }
    }
    return A;
}

int main() {
    int n = 400;
    double p1 = 0.2, p2 = 0.5;
    Matrix A = erdosRenyiGraph(n, p1);
    Matrix B = erdosRenyiGraph(n, p2);
    
    double s = cut_distance(A, B);
    cout << "The cut distance between A and B: " << s << endl;
    
    return 0;
}

