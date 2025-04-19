#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef struct {
    int rows;
    int cols;
    double *data;
} Matrix;

typedef struct {
    double f;
    Matrix g;
} CutNormQuadResult;

typedef struct {
    Matrix x;
    Matrix g;
} OptimizeResult;

Matrix create_matrix(int rows, int cols);
void free_matrix(Matrix *m);
Matrix copy_matrix(const Matrix *m);
Matrix subtract_matrix(const Matrix *A, const Matrix *B);
Matrix divide_matrix(const Matrix *A, double scalar);
double matrix_sum(const Matrix *A);
double* col_sum(const Matrix *A);
double* row_sum(const Matrix *A);
Matrix append_column(Matrix A, const double* col);
Matrix append_row(Matrix A, const double* row);
Matrix transpose_matrix(const Matrix *A);
Matrix mat_multiply(const Matrix *A, const Matrix *B);
double norm_vector(const double* v, int n);
void normalize_vector(double* v, int n);
double* sqrt_vector(const double* v, int n);
double* elementwise_multiply_vectors(const double* a, const double* b, int n);
Matrix broadcast_divide(const Matrix *M, const double* v);
double* col_elementwise_product_sum(const Matrix *A, const Matrix *B);
Matrix broadcast_multiply(const Matrix *M, const double* v);
Matrix elementwise_multiply_matrix(const Matrix *A, const Matrix *B);
Matrix elementwise_subtract(const Matrix *A, const Matrix *B);
Matrix add_matrix(const Matrix *A, const Matrix *B);
double frobenius_norm(const Matrix *M);
double matrix_elementwise_sum(const Matrix *M);
double clamp_val(double x, double lower, double upper);
double rand_normal(void);
Matrix random_normal_matrix(int rows, int cols);
double* col_square_sum(const Matrix *M);
CutNormQuadResult cut_norm_quad(const Matrix *V, const Matrix *A);
OptimizeResult optimize(Matrix x, const Matrix *A,
                        double xtol, double ftol, double gtol,
                        double rho, double eta, double gamma,
                        double tau, int nt, int mxitr);
double cut_norm(const Matrix *A_input);
double cut_distance(const Matrix *G, const Matrix *H);
Matrix erdos_renyi_graph(int n, double p);
void normalize_columns(Matrix *X);

Matrix create_matrix(int rows, int cols) {
    Matrix m;
    m.rows = rows;
    m.cols = cols;
    m.data = (double*)calloc(rows * cols, sizeof(double));
    if (!m.data) { perror("calloc"); exit(EXIT_FAILURE); }
    return m;
}

void free_matrix(Matrix *m) {
    if (m->data) free(m->data);
    m->data = NULL;
    m->rows = m->cols = 0;
}

Matrix copy_matrix(const Matrix *m) {
    Matrix c = create_matrix(m->rows, m->cols);
    for (int i = 0; i < m->rows * m->cols; i++) c.data[i] = m->data[i];
    return c;
}

Matrix subtract_matrix(const Matrix *A, const Matrix *B) {
    Matrix C = create_matrix(A->rows, A->cols);
    for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < A->cols; j++)
            C.data[i*A->cols+j] = A->data[i*A->cols+j] - B->data[i*A->cols+j];
    return C;
}

Matrix divide_matrix(const Matrix *A, double scalar) {
    Matrix C = create_matrix(A->rows, A->cols);
    for (int i = 0; i < A->rows * A->cols; i++)
        C.data[i] = A->data[i] / scalar;
    return C;
}

double matrix_sum(const Matrix *A) {
    double s = 0.0;
    for (int i = 0; i < A->rows * A->cols; i++) s += A->data[i];
    return s;
}

double* col_sum(const Matrix *A) {
    double *sum = (double*)calloc(A->cols, sizeof(double));
    if (!sum) { perror("calloc"); exit(EXIT_FAILURE); }
    for (int j = 0; j < A->cols; j++)
        for (int i = 0; i < A->rows; i++)
            sum[j] += A->data[i*A->cols+j];
    return sum;
}

double* row_sum(const Matrix *A) {
    double *sum = (double*)calloc(A->rows, sizeof(double));
    if (!sum) { perror("calloc"); exit(EXIT_FAILURE); }
    for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < A->cols; j++)
            sum[i] += A->data[i*A->cols+j];
    return sum;
}

Matrix append_column(Matrix A, const double* col) {
    Matrix B = create_matrix(A.rows, A.cols+1);
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < A.cols; j++) B.data[i*(A.cols+1)+j] = A.data[i*A.cols+j];
        B.data[i*(A.cols+1)+A.cols] = col[i];
    }
    free_matrix(&A);
    return B;
}

Matrix append_row(Matrix A, const double* row) {
    Matrix B = create_matrix(A.rows+1, A.cols);
    for (int i = 0; i < A.rows; i++)
        for (int j = 0; j < A.cols; j++)
            B.data[i*A.cols+j] = A.data[i*A.cols+j];
    for (int j = 0; j < A.cols; j++) B.data[A.rows*A.cols+j] = row[j];
    free_matrix(&A);
    return B;
}

Matrix transpose_matrix(const Matrix *A) {
    Matrix T = create_matrix(A->cols, A->rows);
    for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < A->cols; j++)
            T.data[j*A->rows+i] = A->data[i*A->cols+j];
    return T;
}

Matrix mat_multiply(const Matrix *A, const Matrix *B) {
    Matrix C = create_matrix(A->rows, B->cols);
    for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < B->cols; j++) {
            double s = 0.0;
            for (int k = 0; k < A->cols; k++)
                s += A->data[i*A->cols+k] * B->data[k*B->cols+j];
            C.data[i*B->cols+j] = s;
        }
    return C;
}

double norm_vector(const double* v, int n) { double s=0; for(int i=0;i<n;i++) s+=v[i]*v[i]; return sqrt(s); }
void normalize_vector(double* v, int n) { double nrm=norm_vector(v,n); if(nrm>1e-8) for(int i=0;i<n;i++) v[i]/=nrm; }
double* sqrt_vector(const double* v, int n) {
    double *r=(double*)malloc(n*sizeof(double)); if(!r){perror("malloc");exit(EXIT_FAILURE);} for(int i=0;i<n;i++) r[i]=sqrt(v[i]); return r;
}
double* elementwise_multiply_vectors(const double* a, const double* b, int n) {
    double *r=(double*)malloc(n*sizeof(double)); if(!r){perror("malloc");exit(EXIT_FAILURE);} for(int i=0;i<n;i++) r[i]=a[i]*b[i]; return r;
}
Matrix broadcast_divide(const Matrix *M, const double* v) {
    Matrix R=create_matrix(M->rows,M->cols);
    for(int j=0;j<M->cols;j++) for(int i=0;i<M->rows;i++) R.data[i*M->cols+j]=M->data[i*M->cols+j]/v[j];
    return R;
}
double* col_elementwise_product_sum(const Matrix *A, const Matrix *B) {
    double *r=(double*)calloc(A->cols,sizeof(double)); if(!r){perror("calloc");exit(EXIT_FAILURE);} for(int j=0;j<A->cols;j++) for(int i=0;i<A->rows;i++) r[j]+=A->data[i*A->cols+j]*B->data[i*A->cols+j]; return r;
}
Matrix broadcast_multiply(const Matrix *M, const double* v) {
    Matrix R=create_matrix(M->rows,M->cols);
    for(int j=0;j<M->cols;j++) for(int i=0;i<M->rows;i++) R.data[i*M->cols+j]=M->data[i*M->cols+j]*v[j];
    return R;
}
Matrix elementwise_multiply_matrix(const Matrix *A, const Matrix *B) {
    Matrix R=create_matrix(A->rows,A->cols);
    for(int i=0;i<A->rows;i++) for(int j=0;j<A->cols;j++) R.data[i*A->cols+j]=A->data[i*A->cols+j]*B->data[i*A->cols+j];
    return R;
}
Matrix elementwise_subtract(const Matrix *A, const Matrix *B) {
    Matrix R=create_matrix(A->rows,A->cols);
    for(int i=0;i<A->rows;i++) for(int j=0;j<A->cols;j++) R.data[i*A->cols+j]=A->data[i*A->cols+j]-B->data[i*A->cols+j];
    return R;
}
Matrix add_matrix(const Matrix *A, const Matrix *B) {
    Matrix R=create_matrix(A->rows,A->cols);
    for(int i=0;i<A->rows;i++) for(int j=0;j<A->cols;j++) R.data[i*A->cols+j]=A->data[i*A->cols+j]+B->data[i*A->cols+j];
    return R;
}
double frobenius_norm(const Matrix *M) { double s=0; for(int i=0;i<M->rows*M->cols;i++) s+=M->data[i]*M->data[i]; return sqrt(s); }
double matrix_elementwise_sum(const Matrix *M) { double s=0; for(int i=0;i<M->rows*M->cols;i++) s+=M->data[i]; return s; }
double clamp_val(double x,double lo,double hi){ return x<lo?lo:(x>hi?hi:x); }

double rand_normal(void) {
    static int have=0; static double spare;
    if(have){ have=0; return spare; }
    have=1;
    double u,v,s;
    do{ u=rand()/(double)RAND_MAX*2-1; v=rand()/(double)RAND_MAX*2-1; s=u*u+v*v; }while(s>=1||s==0);
    s=sqrt(-2*log(s)/s);
    spare=v*s;
    return u*s;
}

Matrix random_normal_matrix(int rows,int cols){ Matrix M=create_matrix(rows,cols); for(int i=0;i<rows*cols;i++) M.data[i]=rand_normal(); return M; }

double* col_square_sum(const Matrix *M) {
    double *r=(double*)calloc(M->cols,sizeof(double)); if(!r){perror("calloc");exit(EXIT_FAILURE);} for(int j=0;j<M->cols;j++) for(int i=0;i<M->rows;i++) r[j]+=M->data[i*M->cols+j]*M->data[i*M->cols+j]; return r;
}

CutNormQuadResult cut_norm_quad(const Matrix *V, const Matrix *A) {
    int n=A->rows, p=V->rows, n2=V->cols, half=n2/2;
    Matrix U=create_matrix(p,half), V2=create_matrix(p,half);
    for(int i=0;i<p;i++) for(int j=0;j<half;j++){ U.data[i*half+j]=V->data[i*n2+j]; V2.data[i*half+j]=V->data[i*n2+j+half]; }
    Matrix AT=transpose_matrix(A);
    Matrix part1=mat_multiply(&V2,&AT);
    Matrix part2=mat_multiply(&U,A);
    Matrix g=create_matrix(p,2*n);
    for(int i=0;i<p;i++) for(int j=0;j<n;j++){
        g.data[i*2*n+j]   = 2*part1.data[i*n+j];
        g.data[i*2*n+j+n] = 2*part2.data[i*n+j];
    }
    double f=0;
    for(int i=0;i<p;i++) for(int j=0;j<n;j++) f+=g.data[i*2*n+j]*U.data[i*half+j] + g.data[i*2*n+j+n]*V2.data[i*half+j];
    f/=2;
    free_matrix(&AT); free_matrix(&part1); free_matrix(&part2);
    free_matrix(&U); free_matrix(&V2);
    CutNormQuadResult res={f,g};
    return res;
}

OptimizeResult optimize(Matrix x, const Matrix *A,
                        double xtol,double ftol,double gtol,
                        double rho,double eta,double gamma,
                        double tau,int nt,int mxitr) {
    int n=x.rows, p=x.cols;
    double *crit=(double*)calloc(mxitr*3,sizeof(double));
    double *nrx=col_square_sum(&x);
    if(norm_vector(nrx,p)>1e-8){ double *sr=sqrt_vector(nrx,p); Matrix t=broadcast_divide(&x,sr); free_matrix(&x); x=t; free(sr); }
    free(nrx);
    CutNormQuadResult cnq=cut_norm_quad(&x,A);
    double f=cnq.f;
    Matrix g=cnq.g;
    double *xtg=col_elementwise_product_sum(&x,&g);
    double *gg=col_square_sum(&g);
    double *xx=col_square_sum(&x);
    double *xxgg=elementwise_multiply_vectors(xx,gg,p);
    Matrix dtX=broadcast_multiply(&x,xtg);
    Matrix sub=elementwise_subtract(&dtX,&g);
    free_matrix(&dtX); dtX=sub;
    double nrmG=frobenius_norm(&dtX);
    double Q=1, Cval=f, tau0=tau;
    for(int itr=0;itr<mxitr;itr++){
        Matrix xp=copy_matrix(&x);
        double fp=f;
        Matrix gp=copy_matrix(&g);
        Matrix dtXP=copy_matrix(&dtX);
        int nls=1;
        double deriv=rho*nrmG*nrmG;
        while(1){
            double tau2=tau/2;
            double *beta=(double*)malloc(p*sizeof(double));
            for(int j=0;j<p;j++) beta[j]=1+tau2*tau2*(-xtg[j]*xtg[j]+xxgg[j]);
            double *a1=(double*)malloc(p*sizeof(double));
            for(int j=0;j<p;j++){ double tv=1+tau2*xtg[j]; a1[j]=(tv*tv - tau2*tau2*xxgg[j])/beta[j]; }
            double *a2=(double*)malloc(p*sizeof(double));
            for(int j=0;j<p;j++) a2[j] = -tau*xx[j]/beta[j];
            Matrix bp1=broadcast_multiply(&xp,a1);
            Matrix bp2=broadcast_multiply(&gp,a2);
            Matrix xnew=add_matrix(&bp1,&bp2);
            free_matrix(&x); x=xnew;
            CutNormQuadResult cnq2=cut_norm_quad(&x,A);
            f=cnq2.f;
            free_matrix(&g);
            g=cnq2.g;
            free(beta); free(a1); free(a2);
            if(f<=Cval - tau*deriv || nls>=5) break;
            tau*=eta;
            nls++;
        }
        free(xtg);
        xtg=col_elementwise_product_sum(&x,&g);
        free(gg); gg=col_square_sum(&g);
        free(xx); xx=col_square_sum(&x);
        free(xxgg); xxgg=elementwise_multiply_vectors(xx,gg,p);
        free_matrix(&dtX);
        dtX=broadcast_multiply(&x,xtg);
        Matrix sub2=elementwise_subtract(&dtX,&g);
        free_matrix(&dtX); dtX=sub2;
        nrmG=frobenius_norm(&dtX);
        Matrix s=elementwise_subtract(&x,&xp);
        double XDiff=frobenius_norm(&s)/sqrt(n);
        double FDiff=fabs(fp-f)/(fabs(fp)+1);
        crit[itr*3+0]=nrmG;
        crit[itr*3+1]=XDiff;
        crit[itr*3+2]=FDiff;
        int start=itr<nt?0:itr-nt;
        double m0=0,m1=0,m2=0;
        for(int k=start;k<=itr;k++){ m0+=crit[k*3]; m1+=crit[k*3+1]; m2+=crit[k*3+2]; }
        int cnt=itr-start+1;
        m0/=cnt; m1/=cnt; m2/=cnt;
        if((XDiff<xtol && FDiff<ftol) || nrmG<gtol || (m1<10*xtol && m2<10*ftol)){
            free_matrix(&xp); free_matrix(&gp); free_matrix(&dtXP); free_matrix(&s);
            break;
        }
        Matrix y=elementwise_subtract(&dtX,&dtXP);
        Matrix se=elementwise_multiply_matrix(&s,&y);
        double sy=fabs(matrix_sum(&se));
        free_matrix(&se);
        tau=tau0;
        if(sy>0){
            Matrix ss=elementwise_multiply_matrix(&s,&s);
            double sum_ss=matrix_elementwise_sum(&ss);
            free_matrix(&ss);
            Matrix yy=elementwise_multiply_matrix(&y,&y);
            double sum_yy=matrix_elementwise_sum(&yy);
            free_matrix(&yy);
            tau = ((itr%2)==0? sum_ss/sy : sy/sum_yy);
            tau = clamp_val(tau,1e-20,1e20);
        }
        double Qp=Q;
        Q = gamma*Qp + 1;
        Cval = (gamma*Qp*Cval + f)/Q;
        free_matrix(&xp); free_matrix(&gp); free_matrix(&dtXP); free_matrix(&s);
    }
    free(xtg); free(gg); free(xx); free(xxgg); free(crit);
    free_matrix(&dtX);
    OptimizeResult R={x,g};
    return R;
}

double cut_norm(const Matrix *A_input) {
    Matrix A = copy_matrix(A_input);
    int n1 = A.rows;
    double *A_col = col_sum(&A);
    double *A_row = row_sum(&A);
    double A_tot = matrix_sum(&A);

    double *negA_row = (double*)calloc(n1, sizeof(*negA_row));
    if (!negA_row) { perror("calloc"); exit(EXIT_FAILURE); }
    for (int i = 0; i < n1; i++) negA_row[i] = -A_row[i];
    free(A_row);
    A = append_column(A, negA_row);
    free(negA_row);

    double* newRow = (double*)malloc((n1+1)*sizeof(*newRow));
    if (!newRow) { perror("malloc"); exit(EXIT_FAILURE); }
    for (int i = 0; i < n1; i++) newRow[i] = -A_col[i];
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
            U.data[i*half+j]   = x.data[i*n2+j];
            V.data[i*half+j]   = x.data[i*n2+j+half];
        }
    }
    Matrix U_trans = transpose_matrix(&U);
    Matrix prod    = mat_multiply(&U_trans, &V);
    double sdp = 0.0;
    int r = A.rows, c = A.cols;
    for (int i = 0; i < r && i < prod.rows; i++)
        for (int j = 0; j < c && j < prod.cols; j++)
            sdp += A.data[i*A.cols+j] * prod.data[i*prod.cols+j];

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
    Matrix A    = divide_matrix(&diff, (double)(G->rows * G->rows));
    free_matrix(&diff);
    double cd = cut_norm(&A);
    free_matrix(&A);
    return cd;
}

Matrix erdos_renyi_graph(int n, double p) {
    Matrix A = create_matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            if ((rand()/(double)RAND_MAX) < p) {
                A.data[i*n+j] = 1.0;
                A.data[j*n+i] = 1.0;
            }
        }
    }
    return A;
}

void normalize_columns(Matrix *X) {
    for (int j = 0; j < X->cols; j++) {
        double sum_sq = 0.0;
        for (int i = 0; i < X->rows; i++)
            sum_sq += X->data[i*X->cols+j] * X->data[i*X->cols+j];
        double nrm = sqrt(sum_sq);
        if (nrm > 1e-8)
            for (int i = 0; i < X->rows; i++)
                X->data[i*X->cols+j] /= nrm;
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
