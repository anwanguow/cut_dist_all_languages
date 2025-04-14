#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef struct {
    int rows;
    int cols;
    double *data;
} Matrix;

void normalize_columns(Matrix *X);

Matrix create_matrix(int rows, int cols) {
    Matrix m;
    m.rows = rows;
    m.cols = cols;
    m.data = (double*)calloc(rows * cols, sizeof(double));
    return m;
}

void free_matrix(Matrix *m) {
    if (m->data) {
        free(m->data);
    }
    m->data = NULL;
    m->rows = 0;
    m->cols = 0;
}

Matrix copy_matrix(const Matrix *m) {
    Matrix copy = create_matrix(m->rows, m->cols);
    for (int i = 0; i < m->rows * m->cols; i++) {
        copy.data[i] = m->data[i];
    }
    return copy;
}

Matrix subtract_matrix(const Matrix *A, const Matrix *B) {
    Matrix C = create_matrix(A->rows, A->cols);
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            C.data[i * A->cols + j] = A->data[i * A->cols + j] - B->data[i * A->cols + j];
        }
    }
    return C;
}

Matrix divide_matrix(const Matrix *A, double scalar) {
    Matrix C = create_matrix(A->rows, A->cols);
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            C.data[i * A->cols + j] = A->data[i * A->cols + j] / scalar;
        }
    }
    return C;
}

double matrix_sum(const Matrix *A) {
    double s = 0.0;
    for (int i = 0; i < A->rows * A->cols; i++) {
        s += A->data[i];
    }
    return s;
}

double* col_sum(const Matrix *A) {
    double *sum = (double*)calloc(A->cols, sizeof(double));
    for (int j = 0; j < A->cols; j++) {
        for (int i = 0; i < A->rows; i++) {
            sum[j] += A->data[i * A->cols + j];
        }
    }
    return sum;
}

double* row_sum(const Matrix *A) {
    double *sum = (double*)calloc(A->rows, sizeof(double));
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            sum[i] += A->data[i * A->cols + j];
        }
    }
    return sum;
}

Matrix append_column(Matrix A, const double* col) {
    Matrix B = create_matrix(A.rows, A.cols + 1);
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < A.cols; j++) {
            B.data[i * (A.cols + 1) + j] = A.data[i * A.cols + j];
        }
        B.data[i * (A.cols + 1) + A.cols] = col[i];
    }
    free_matrix(&A);
    return B;
}

Matrix append_row(Matrix A, const double* row) {
    Matrix B = create_matrix(A.rows + 1, A.cols);
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < A.cols; j++) {
            B.data[i * A.cols + j] = A.data[i * A.cols + j];
        }
    }
    for (int j = 0; j < A.cols; j++) {
        B.data[A.rows * A.cols + j] = row[j];
    }
    free_matrix(&A);
    return B;
}

Matrix transpose_matrix(const Matrix *A) {
    Matrix T = create_matrix(A->cols, A->rows);
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            T.data[j * A->rows + i] = A->data[i * A->cols + j];
        }
    }
    return T;
}

Matrix mat_multiply(const Matrix *A, const Matrix *B) {
    Matrix C = create_matrix(A->rows, B->cols);
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < B->cols; j++) {
            double sum = 0.0;
            for (int k = 0; k < A->cols; k++) {
                sum += A->data[i * A->cols + k] * B->data[k * B->cols + j];
            }
            C.data[i * B->cols + j] = sum;
        }
    }
    return C;
}

double norm_vector(const double* v, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) {
        s += v[i] * v[i];
    }
    return sqrt(s);
}

void normalize_vector(double* v, int n) {
    double nrm = norm_vector(v, n);
    if (nrm > 1e-8) {
        for (int i = 0; i < n; i++) {
            v[i] /= nrm;
        }
    }
}

double* sqrt_vector(const double* v, int n) {
    double *res = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        res[i] = sqrt(v[i]);
    }
    return res;
}

double* elementwise_multiply_vectors(const double* a, const double* b, int n) {
    double *res = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        res[i] = a[i] * b[i];
    }
    return res;
}

Matrix broadcast_divide(const Matrix *M, const double* v) {
    Matrix res = create_matrix(M->rows, M->cols);
    for (int j = 0; j < M->cols; j++) {
        for (int i = 0; i < M->rows; i++) {
            res.data[i * M->cols + j] = M->data[i * M->cols + j] / v[j];
        }
    }
    return res;
}

double* col_elementwise_product_sum(const Matrix *A, const Matrix *B) {
    double *res = (double*)calloc(A->cols, sizeof(double));
    for (int j = 0; j < A->cols; j++) {
        for (int i = 0; i < A->rows; i++) {
            res[j] += A->data[i * A->cols + j] * B->data[i * A->cols + j];
        }
    }
    return res;
}

Matrix broadcast_multiply(const Matrix *M, const double* v) {
    Matrix res = create_matrix(M->rows, M->cols);
    for (int j = 0; j < M->cols; j++) {
        for (int i = 0; i < M->rows; i++) {
            res.data[i * M->cols + j] = M->data[i * M->cols + j] * v[j];
        }
    }
    return res;
}

Matrix elementwise_multiply_matrix(const Matrix *A, const Matrix *B) {
    Matrix res = create_matrix(A->rows, A->cols);
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            res.data[i * A->cols + j] = A->data[i * A->cols + j] * B->data[i * A->cols + j];
        }
    }
    return res;
}

Matrix elementwise_subtract(const Matrix *A, const Matrix *B) {
    Matrix res = create_matrix(A->rows, A->cols);
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            res.data[i * A->cols + j] = A->data[i * A->cols + j] - B->data[i * A->cols + j];
        }
    }
    return res;
}

Matrix add_matrix(const Matrix *A, const Matrix *B) {
    Matrix res = create_matrix(A->rows, A->cols);
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            res.data[i * A->cols + j] = A->data[i * A->cols + j] + B->data[i * A->cols + j];
        }
    }
    return res;
}

double frobenius_norm(const Matrix *M) {
    double sum = 0.0;
    for (int i = 0; i < M->rows * M->cols; i++) {
        sum += M->data[i] * M->data[i];
    }
    return sqrt(sum);
}

double matrix_elementwise_sum(const Matrix *M) {
    double s = 0.0;
    for (int i = 0; i < M->rows * M->cols; i++) {
        s += M->data[i];
    }
    return s;
}

double clamp_val(double x, double lower, double upper) {
    if (x < lower) return lower;
    if (x > upper) return upper;
    return x;
}

double rand_normal() {
    static int haveSpare = 0;
    static double spare;
    if (haveSpare) {
        haveSpare = 0;
        return spare;
    }
    haveSpare = 1;
    double u, v, s;
    do {
        u = (rand() / ((double)RAND_MAX)) * 2.0 - 1.0;
        v = (rand() / ((double)RAND_MAX)) * 2.0 - 1.0;
        s = u * u + v * v;
    } while (s >= 1 || s == 0);
    s = sqrt(-2.0 * log(s) / s);
    spare = v * s;
    return u * s;
}

Matrix random_normal_matrix(int rows, int cols) {
    Matrix M = create_matrix(rows, cols);
    for (int i = 0; i < rows * cols; i++) {
        M.data[i] = rand_normal();
    }
    return M;
}

double* col_square_sum(const Matrix *M) {
    double *res = (double*)calloc(M->cols, sizeof(double));
    for (int j = 0; j < M->cols; j++) {
        for (int i = 0; i < M->rows; i++) {
            res[j] += M->data[i * M->cols + j] * M->data[i * M->cols + j];
        }
    }
    return res;
}

typedef struct {
    double f;
    Matrix g;
} CutNormQuadResult;

CutNormQuadResult cut_norm_quad(const Matrix *V, const Matrix *A) {
    int n = A->rows;
    int p = V->rows;
    int n2 = V->cols;
    int half = n2 / 2;
    Matrix U = create_matrix(p, half);
    Matrix V2 = create_matrix(p, half);
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < half; j++) {
            U.data[i * half + j] = V->data[i * n2 + j];
            V2.data[i * half + j] = V->data[i * n2 + (j + half)];
        }
    }
    Matrix A_T = transpose_matrix(A);
    Matrix part1 = mat_multiply(&V2, &A_T);
    Matrix part2 = mat_multiply(&U, A);
    Matrix g = create_matrix(p, 2 * n);
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < n; j++) {
            g.data[i * (2 * n) + j] = 2.0 * part1.data[i * n + j];
            g.data[i * (2 * n) + j + n] = 2.0 * part2.data[i * n + j];
        }
    }
    double f_val = 0.0;
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < n; j++) {
            f_val += g.data[i * (2 * n) + j] * U.data[i * half + j] +
                     g.data[i * (2 * n) + j + n] * V2.data[i * half + j];
        }
    }
    f_val /= 2.0;
    free_matrix(&A_T);
    free_matrix(&part1);
    free_matrix(&part2);
    free_matrix(&U);
    free_matrix(&V2);
    CutNormQuadResult result;
    result.f = f_val;
    result.g = g;
    return result;
}

typedef struct {
    Matrix x;
    Matrix g;
} OptimizeResult;

OptimizeResult optimize(Matrix x, const Matrix *A, double xtol, double ftol, double gtol, double rho,
                          double eta, double gamma, double tau, int nt, int mxitr) {
    double* crit = (double*)calloc(mxitr * 3, sizeof(double));
    int n = x.rows;
    int p = x.cols;
    double* nrmx = col_square_sum(&x);
    double norm_nrmx = norm_vector(nrmx, x.cols);
    if (norm_nrmx > 1e-8) {
        double* sqrt_nrmx = sqrt_vector(nrmx, x.cols);
        Matrix new_x = broadcast_divide(&x, sqrt_nrmx);
        free_matrix(&x);
        x = new_x;
        free(sqrt_nrmx);
    }
    free(nrmx);
    double f;
    CutNormQuadResult cnq = cut_norm_quad(&x, A);
    f = cnq.f;
    Matrix g = cnq.g;
    double* xtg = col_elementwise_product_sum(&x, &g);
    double* gg = col_square_sum(&g);
    double* xx = col_square_sum(&x);
    double* xxgg = elementwise_multiply_vectors(xx, gg, x.cols);
    Matrix dtX = broadcast_multiply(&x, xtg);
    Matrix temp = elementwise_subtract(&dtX, &g);
    free_matrix(&dtX);
    dtX = temp;
    double nrmG = frobenius_norm(&dtX);
    double Q = 1.0;
    double Cval = f;
    double tau_orig = tau;
    tau = tau_orig;
    for (int itr = 0; itr < mxitr; itr++) {
        Matrix xp = copy_matrix(&x);
        double fp = f;
        Matrix gp = copy_matrix(&g);
        Matrix dtXP = copy_matrix(&dtX);
        int nls = 1;
        double deriv = rho * nrmG * nrmG;
        while (1) {
            double tau2 = tau / 2.0;
            double* beta = (double*)malloc(p * sizeof(double));
            for (int j = 0; j < p; j++) {
                beta[j] = 1.0 + tau2 * tau2 * (- xtg[j] * xtg[j] + xxgg[j]);
            }
            double* a1 = (double*)malloc(p * sizeof(double));
            for (int j = 0; j < p; j++) {
                double temp_val = 1.0 + tau2 * xtg[j];
                a1[j] = ((temp_val * temp_val) - tau2 * tau2 * xxgg[j]) / beta[j];
            }
            double* a2 = (double*)malloc(p * sizeof(double));
            for (int j = 0; j < p; j++) {
                a2[j] = -tau * xx[j] / beta[j];
            }
            Matrix bp1 = broadcast_multiply(&xp, a1);
            Matrix bp2 = broadcast_multiply(&gp, a2);
            Matrix x_new = add_matrix(&bp1, &bp2);
            free_matrix(&bp1);
            free_matrix(&bp2);
            free(beta);
            free(a1);
            free(a2);
            free_matrix(&x);
            x = x_new;
            CutNormQuadResult cnq_new = cut_norm_quad(&x, A);
            f = cnq_new.f;
            free_matrix(&g);
            g = cnq_new.g;
            if (f <= Cval - tau * deriv || nls >= 5)
                break;
            tau = eta * tau;
            nls = nls + 1;
        }
        free(xtg);
        xtg = col_elementwise_product_sum(&x, &g);
        free(gg);
        gg = col_square_sum(&g);
        free(xx);
        xx = col_square_sum(&x);
        free(xxgg);
        xxgg = elementwise_multiply_vectors(xx, gg, x.cols);
        free_matrix(&dtX);
        dtX = broadcast_multiply(&x, xtg);
        Matrix sub = elementwise_subtract(&dtX, &g);
        free_matrix(&dtX);
        dtX = sub;
        nrmG = frobenius_norm(&dtX);
        Matrix s = elementwise_subtract(&x, &xp);
        double XDiff = frobenius_norm(&s) / sqrt(n);
        double FDiff = fabs(fp - f) / (fabs(fp) + 1.0);
        crit[itr * 3 + 0] = nrmG;
        crit[itr * 3 + 1] = XDiff;
        crit[itr * 3 + 2] = FDiff;
        int start = (itr - ((nt < itr) ? nt : itr));
        double mcrit0 = 0.0, mcrit1 = 0.0, mcrit2 = 0.0;
        int count = 0;
        for (int k = start; k <= itr; k++) {
            mcrit0 += crit[k * 3 + 0];
            mcrit1 += crit[k * 3 + 1];
            mcrit2 += crit[k * 3 + 2];
            count++;
        }
        mcrit0 /= count;
        mcrit1 /= count;
        mcrit2 /= count;
        if ((XDiff < xtol && FDiff < ftol) || nrmG < gtol ||
            (mcrit1 < 10 * xtol && mcrit2 < 10 * ftol)) {
            free_matrix(&xp);
            free_matrix(&gp);
            free_matrix(&dtXP);
            free_matrix(&s);
            break;
        }
        Matrix y = elementwise_subtract(&dtX, &dtXP);
        Matrix s_elem = elementwise_multiply_matrix(&s, &y);
        double sy = fabs(matrix_sum(&s_elem));
        free_matrix(&s_elem);
        tau = tau_orig;
        if (sy > 0) {
            Matrix s_sq = elementwise_multiply_matrix(&s, &s);
            double sum_ss = matrix_elementwise_sum(&s_sq);
            free_matrix(&s_sq);
            Matrix y_sq = elementwise_multiply_matrix(&y, &y);
            double sum_yy = matrix_elementwise_sum(&y_sq);
            free_matrix(&y_sq);
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
        free_matrix(&xp);
        free_matrix(&gp);
        free_matrix(&dtXP);
        free_matrix(&s);
    }
    free(xtg);
    free(gg);
    free(xx);
    free(xxgg);
    free(crit);
    OptimizeResult res;
    res.x = x;
    res.g = g;
    free_matrix(&dtX);
    return res;
}

double cut_norm(const Matrix *A_input) {
    Matrix A = copy_matrix(A_input);
    int n1 = A.rows;
    double* A_col = col_sum(&A);
    double* A_row = row_sum(&A);
    double A_tot = matrix_sum(&A);
    double* negA_row = (double*)malloc(n1 * sizeof(double));
    for (int i = 0; i < n1; i++) {
        negA_row[i] = -A_row[i];
    }
    free(A_row);
    A = append_column(A, negA_row);
    free(negA_row);
    double* newRow = (double*)malloc((n1 + 1) * sizeof(double));
    for (int i = 0; i < n1; i++) {
        newRow[i] = -A_col[i];
    }
    newRow[n1] = A_tot;
    free(A_col);
    A = append_row(A, newRow);
    free(newRow);
    double temp_val = sqrt(2.0 * n1) / 2.0;
    int p_val = (int)round(temp_val);
    if (p_val > 100) p_val = 100;
    if (p_val < 1) p_val = 1;
    int n2 = 2 * n1 + 2;
    Matrix x0 = random_normal_matrix(p_val, n2);
    normalize_columns(&x0);
    OptimizeResult optRes = optimize(x0, &A, 1e-8, 1e-10, 1e-8, 1e-4, 0.1, 0.85, 1e-3, 5, 600);
    Matrix x = optRes.x;
    free_matrix(&optRes.g);
    int half = n2 / 2;
    Matrix U = create_matrix(p_val, half);
    Matrix V = create_matrix(p_val, half);
    for (int i = 0; i < p_val; i++) {
        for (int j = 0; j < half; j++) {
            U.data[i * half + j] = x.data[i * n2 + j];
            V.data[i * half + j] = x.data[i * n2 + (j + half)];
        }
    }
    Matrix U_trans = transpose_matrix(&U);
    Matrix prod = mat_multiply(&U_trans, &V);
    double sdp = 0.0;
    int r = A.rows;
    int c = A.cols;
    int pr = prod.rows;
    int pc = prod.cols;
    for (int i = 0; i < r && i < pr; i++) {
        for (int j = 0; j < c && j < pc; j++) {
            sdp += A.data[i * A.cols + j] * prod.data[i * pc + j];
        }
    }
    double sdp_result = fabs(sdp) / 4.0;
    free_matrix(&A);
    free_matrix(&x);
    free_matrix(&U);
    free_matrix(&V);
    free_matrix(&U_trans);
    free_matrix(&prod);
    return sdp_result;
}

double cut_distance(const Matrix *G, const Matrix *H) {
    Matrix diff = subtract_matrix(G, H);
    Matrix A = divide_matrix(&diff, (double)(G->rows * G->rows));
    free_matrix(&diff);
    double norm_val = cut_norm(&A);
    free_matrix(&A);
    return norm_val;
}

Matrix erdos_renyi_graph(int n, double p) {
    Matrix A = create_matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double r = ((double)rand()) / RAND_MAX;
            if (r < p) {
                A.data[i * n + j] = 1.0;
                A.data[j * n + i] = 1.0;
            }
        }
    }
    return A;
}

void normalize_columns(Matrix *X) {
    int rows = X->rows;
    int cols = X->cols;
    for (int j = 0; j < cols; j++){
        double sum_sq = 0.0;
        for (int i = 0; i < rows; i++){
            sum_sq += X->data[i * cols + j] * X->data[i * cols + j];
        }
        double nrm = sqrt(sum_sq);
        if (nrm > 1e-8) {
            for (int i = 0; i < rows; i++){
                X->data[i * cols + j] /= nrm;
            }
        }
    }
}

int main(void) {
    srand((unsigned)time(NULL));
    int n = 400;
    double p1 = 0.2, p2 = 0.5;
    Matrix A = erdos_renyi_graph(n, p1);
    Matrix B = erdos_renyi_graph(n, p2);
    double s = cut_distance(&A, &B);
    printf("The cut distance between A and B: %f\n", s);
    free_matrix(&A);
    free_matrix(&B);
    return 0;
}

