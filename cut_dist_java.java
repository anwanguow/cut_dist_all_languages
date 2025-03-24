package cut_dist;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.lang.Math;

public class CutDistance {
    public static class Pair<F, S> {
        public F first;
        public S second;
        public Pair(F first, S second) {
            this.first = first;
            this.second = second;
        }
    }

    private static final Random rand = new Random();

    public static int numRows(List<List<Double>> M) {
        return M.size();
    }

    public static int numCols(List<List<Double>> M) {
        return M.isEmpty() ? 0 : M.get(0).size();
    }

    public static List<List<Double>> subtractMatrix(List<List<Double>> A, List<List<Double>> B) {
        int r = numRows(A), c = numCols(A);
        List<List<Double>> C = new ArrayList<>();
        for (int i = 0; i < r; i++) {
            ArrayList<Double> row = new ArrayList<>();
            for (int j = 0; j < c; j++) {
                row.add(A.get(i).get(j) - B.get(i).get(j));
            }
            C.add(row);
        }
        return C;
    }

    public static List<List<Double>> divideMatrix(List<List<Double>> A, double scalar) {
        int r = numRows(A), c = numCols(A);
        List<List<Double>> C = new ArrayList<>();
        for (int i = 0; i < r; i++) {
            ArrayList<Double> row = new ArrayList<>();
            for (int j = 0; j < c; j++) {
                row.add(A.get(i).get(j) / scalar);
            }
            C.add(row);
        }
        return C;
    }

    public static double matrixSum(List<List<Double>> A) {
        double s = 0.0;
        for (List<Double> row : A) {
            for (double val : row) {
                s += val;
            }
        }
        return s;
    }

    public static List<Double> colSum(List<List<Double>> A) {
        int r = numRows(A), c = numCols(A);
        List<Double> sum = new ArrayList<>();
        for (int j = 0; j < c; j++) {
            sum.add(0.0);
        }
        for (int j = 0; j < c; j++){
            for (int i = 0; i < r; i++){
                sum.set(j, sum.get(j) + A.get(i).get(j));
            }
        }
        return sum;
    }

    public static List<Double> rowSum(List<List<Double>> A) {
        int r = numRows(A), c = numCols(A);
        List<Double> sum = new ArrayList<>();
        for (int i = 0; i < r; i++){
            double s = 0.0;
            for (int j = 0; j < c; j++){
                s += A.get(i).get(j);
            }
            sum.add(s);
        }
        return sum;
    }

    public static void appendColumn(List<List<Double>> A, List<Double> col) {
        int r = numRows(A);
        for (int i = 0; i < r; i++){
            A.get(i).add(col.get(i));
        }
    }

    public static void appendRow(List<List<Double>> A, List<Double> row) {
        A.add(row);
    }

    public static List<List<Double>> transpose(List<List<Double>> A) {
        int r = numRows(A), c = numCols(A);
        List<List<Double>> T = new ArrayList<>();
        for (int j = 0; j < c; j++){
            ArrayList<Double> row = new ArrayList<>();
            for (int i = 0; i < r; i++){
                row.add(A.get(i).get(j));
            }
            T.add(row);
        }
        return T;
    }

    public static List<List<Double>> matMultiply(List<List<Double>> A, List<List<Double>> B) {
        int r = numRows(A), c = numCols(B), inner = numCols(A);
        List<List<Double>> C = new ArrayList<>();
        for (int i = 0; i < r; i++){
            ArrayList<Double> row = new ArrayList<>();
            for (int j = 0; j < c; j++){
                double sum = 0.0;
                for (int k = 0; k < inner; k++){
                    sum += A.get(i).get(k) * B.get(k).get(j);
                }
                row.add(sum);
            }
            C.add(row);
        }
        return C;
    }

    public static double normVector(List<Double> v) {
        double s = 0.0;
        for (double x : v)
            s += x * x;
        return Math.sqrt(s);
    }

    public static void normalizeVector(List<Double> v) {
        double n = normVector(v);
        if (n > 1e-8) {
            for (int i = 0; i < v.size(); i++) {
                v.set(i, v.get(i) / n);
            }
        }
    }

    public static void normalizeColumns(List<List<Double>> X) {
        int rows = numRows(X), cols = numCols(X);
        for (int j = 0; j < cols; j++){
            double sum_sq = 0.0;
            for (int i = 0; i < rows; i++){
                sum_sq += X.get(i).get(j) * X.get(i).get(j);
            }
            double nrm = Math.sqrt(sum_sq);
            if (nrm > 1e-8) {
                for (int i = 0; i < rows; i++){
                    X.get(i).set(j, X.get(i).get(j) / nrm);
                }
            }
        }
    }

    public static List<List<Double>> randomNormalMatrix(int rows, int cols) {
        List<List<Double>> M = new ArrayList<>();
        for (int i = 0; i < rows; i++){
            ArrayList<Double> row = new ArrayList<>();
            for (int j = 0; j < cols; j++){
                row.add(rand.nextGaussian());
            }
            M.add(row);
        }
        return M;
    }

    public static List<Double> colSquareSum(List<List<Double>> M) {
        int rows = numRows(M), cols = numCols(M);
        List<Double> res = new ArrayList<>();
        for (int j = 0; j < cols; j++){
            res.add(0.0);
        }
        for (int j = 0; j < cols; j++){
            for (int i = 0; i < rows; i++){
                res.set(j, res.get(j) + M.get(i).get(j) * M.get(i).get(j));
            }
        }
        return res;
    }

    public static List<Double> sqrtVector(List<Double> v) {
        List<Double> res = new ArrayList<>();
        for (double x : v) {
            res.add(Math.sqrt(x));
        }
        return res;
    }

    public static List<List<Double>> broadcastDivide(List<List<Double>> M, List<Double> v) {
        int rows = numRows(M), cols = numCols(M);
        List<List<Double>> res = deepCopyMatrix(M);
        for (int j = 0; j < cols; j++){
            for (int i = 0; i < rows; i++){
                res.get(i).set(j, res.get(i).get(j) / v.get(j));
            }
        }
        return res;
    }

    public static List<Double> colElementwiseProductSum(List<List<Double>> A, List<List<Double>> B) {
        int rows = numRows(A), cols = numCols(A);
        List<Double> res = new ArrayList<>();
        for (int j = 0; j < cols; j++){
            res.add(0.0);
        }
        for (int j = 0; j < cols; j++){
            for (int i = 0; i < rows; i++){
                res.set(j, res.get(j) + A.get(i).get(j) * B.get(i).get(j));
            }
        }
        return res;
    }

    public static List<Double> elementwiseMultiplyVector(List<Double> a, List<Double> b) {
        int n = a.size();
        List<Double> res = new ArrayList<>();
        for (int i = 0; i < n; i++){
            res.add(a.get(i) * b.get(i));
        }
        return res;
    }

    public static List<List<Double>> elementwiseMultiplyMatrix(List<List<Double>> A, List<List<Double>> B) {
        int rows = numRows(A), cols = numCols(A);
        List<List<Double>> res = new ArrayList<>();
        for (int i = 0; i < rows; i++){
            ArrayList<Double> row = new ArrayList<>();
            for (int j = 0; j < cols; j++){
                row.add(A.get(i).get(j) * B.get(i).get(j));
            }
            res.add(row);
        }
        return res;
    }

    public static List<List<Double>> broadcastMultiply(List<List<Double>> M, List<Double> v) {
        int rows = numRows(M), cols = numCols(M);
        List<List<Double>> res = deepCopyMatrix(M);
        for (int j = 0; j < cols; j++){
            for (int i = 0; i < rows; i++){
                res.get(i).set(j, res.get(i).get(j) * v.get(j));
            }
        }
        return res;
    }

    public static List<List<Double>> elementwiseSubtract(List<List<Double>> A, List<List<Double>> B) {
        int rows = numRows(A), cols = numCols(A);
        List<List<Double>> res = new ArrayList<>();
        for (int i = 0; i < rows; i++){
            ArrayList<Double> row = new ArrayList<>();
            for (int j = 0; j < cols; j++){
                row.add(A.get(i).get(j) - B.get(i).get(j));
            }
            res.add(row);
        }
        return res;
    }

    public static List<List<Double>> add(List<List<Double>> A, List<List<Double>> B) {
        int rows = numRows(A), cols = numCols(A);
        List<List<Double>> res = new ArrayList<>();
        for (int i = 0; i < rows; i++){
            ArrayList<Double> row = new ArrayList<>();
            for (int j = 0; j < cols; j++){
                row.add(A.get(i).get(j) + B.get(i).get(j));
            }
            res.add(row);
        }
        return res;
    }

    public static double frobeniusNorm(List<List<Double>> M) {
        double sum = 0.0;
        for (List<Double> row : M)
            for (double val : row)
                sum += val * val;
        return Math.sqrt(sum);
    }

    public static double matrixElementwiseSum(List<List<Double>> M) {
        double s = 0.0;
        for (List<Double> row : M)
            for (double val : row)
                s += val;
        return s;
    }

    public static double clamp_val(double x, double lower, double upper) {
        return Math.max(lower, Math.min(x, upper));
    }

    public static Pair<Double, List<List<Double>>> cut_norm_quad(List<List<Double>> V, List<List<Double>> A) {
        int p = numRows(V);
        int n2 = numCols(V);
        int half = n2 / 2;
        List<List<Double>> U = new ArrayList<>();
        List<List<Double>> V2 = new ArrayList<>();
        for (int i = 0; i < p; i++){
            ArrayList<Double> rowU = new ArrayList<>();
            ArrayList<Double> rowV2 = new ArrayList<>();
            for (int j = 0; j < half; j++){
                rowU.add(V.get(i).get(j));
                rowV2.add(V.get(i).get(j + half));
            }
            U.add(rowU);
            V2.add(rowV2);
        }
        List<List<Double>> A_T = transpose(A);
        List<List<Double>> part1 = matMultiply(V2, A_T);
        List<List<Double>> part2 = matMultiply(U, A);
        
        int n = numRows(A);
        List<List<Double>> g = new ArrayList<>();
        for (int i = 0; i < p; i++){
            ArrayList<Double> row = new ArrayList<>();
            for (int j = 0; j < 2 * n; j++){
                row.add(0.0);
            }
            g.add(row);
        }
        
        for (int i = 0; i < p; i++){
            for (int j = 0; j < n; j++){
                g.get(i).set(j, 2 * part1.get(i).get(j));
                g.get(i).set(j + n, 2 * part2.get(i).get(j));
            }
        }
        
        double f = 0.0;
        for (int i = 0; i < p; i++){
            for (int j = 0; j < n; j++){
                f += g.get(i).get(j) * U.get(i).get(j) + g.get(i).get(j + n) * V2.get(i).get(j);
            }
        }
        f /= 2.0;
        return new Pair<>(f, g);
    }

    public static Pair<List<List<Double>>, List<List<Double>>> optimize(List<List<Double>> x, List<List<Double>> A,
            double xtol, double ftol, double gtol, double rho,
            double eta, double gamma, double tau,
            int nt, int mxitr) {
        List<List<Double>> crit = new ArrayList<>();
        for (int i = 0; i < mxitr; i++){
            ArrayList<Double> row = new ArrayList<>();
            for (int j = 0; j < 3; j++){
                row.add(1.0);
            }
            crit.add(row);
        }
        
        int n = numRows(x);
        int p = numCols(x);
        List<Double> nrmx = colSquareSum(x);
        if (normVector(nrmx) > 1e-8) {
            List<Double> sqrt_nrmx = sqrtVector(nrmx);
            x = broadcastDivide(x, sqrt_nrmx);
        }
        double f;
        List<List<Double>> g;
        Pair<Double, List<List<Double>>> temp = cut_norm_quad(x, A);
        f = temp.first;
        g = temp.second;
        List<Double> xtg = colElementwiseProductSum(x, g);
        List<Double> gg = colSquareSum(g);
        List<Double> xx = colSquareSum(x);
        List<Double> xxgg = elementwiseMultiplyVector(xx, gg);
        List<List<Double>> dtX = broadcastMultiply(x, xtg);
        dtX = elementwiseSubtract(dtX, g);
        double nrmG = frobeniusNorm(dtX);
        double Q = 1.0;
        double Cval = f;
        double tau_orig = tau;
        tau = tau_orig;
        
        for (int itr = 0; itr < mxitr; itr++){
            List<List<Double>> xp = deepCopyMatrix(x);
            double fp = f;
            List<List<Double>> gp = deepCopyMatrix(g);
            List<List<Double>> dtXP = deepCopyMatrix(dtX);
            int nls = 1;
            double deriv = rho * nrmG * nrmG;
            while (true) {
                double tau2 = tau / 2.0;
                List<Double> beta = new ArrayList<>();
                for (int j = 0; j < p; j++){
                    beta.add(1.0 + tau2 * tau2 * (- xtg.get(j) * xtg.get(j) + xxgg.get(j)));
                }
                List<Double> a1 = new ArrayList<>();
                for (int j = 0; j < p; j++){
                    double tempVal = 1.0 + tau2 * xtg.get(j);
                    a1.add(((tempVal * tempVal) - tau2 * tau2 * xxgg.get(j)) / beta.get(j));
                }
                List<Double> a2 = new ArrayList<>();
                for (int j = 0; j < p; j++){
                    a2.add(-tau * xx.get(j) / beta.get(j));
                }
                List<List<Double>> x_new = add(broadcastMultiply(xp, a1), broadcastMultiply(gp, a2));
                x = x_new;
                temp = cut_norm_quad(x, A);
                f = temp.first;
                g = temp.second;
                if (f <= Cval - tau * deriv || nls >= 5) {
                    break;
                }
                tau = eta * tau;
                nls = nls + 1;
            }
            xtg = colElementwiseProductSum(x, g);
            gg = colSquareSum(g);
            xx = colSquareSum(x);
            xxgg = elementwiseMultiplyVector(xx, gg);
            dtX = broadcastMultiply(x, xtg);
            dtX = elementwiseSubtract(dtX, g);
            nrmG = frobeniusNorm(dtX);
            List<List<Double>> s = elementwiseSubtract(x, xp);
            double XDiff = frobeniusNorm(s) / Math.sqrt(n);
            double FDiff = Math.abs(fp - f) / (Math.abs(fp) + 1.0);
            crit.get(itr).set(0, nrmG);
            crit.get(itr).set(1, XDiff);
            crit.get(itr).set(2, FDiff);
            int start = Math.max(0, itr - Math.min(nt, itr));
            List<Double> mcrit = new ArrayList<>();
            for (int j = 0; j < 3; j++){
                mcrit.add(0.0);
            }
            int count = 0;
            for (int k = start; k <= itr; k++){
                for (int j = 0; j < 3; j++){
                    mcrit.set(j, mcrit.get(j) + crit.get(k).get(j));
                }
                count++;
            }
            for (int j = 0; j < 3; j++){
                mcrit.set(j, mcrit.get(j) / count);
            }
            if ((XDiff < xtol && FDiff < ftol) || nrmG < gtol || (mcrit.get(1) < 10 * xtol && mcrit.get(2) < 10 * ftol)) {
                break;
            }
            List<List<Double>> y = elementwiseSubtract(dtX, dtXP);
            List<List<Double>> sy_mat = elementwiseMultiplyMatrix(s, y);
            double sy = Math.abs(matrixSum(sy_mat));
            tau = tau_orig;
            if (sy > 0) {
                double sum_ss = matrixElementwiseSum(elementwiseMultiplyMatrix(s, s));
                double sum_yy = matrixElementwiseSum(elementwiseMultiplyMatrix(y, y));
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
        return new Pair<>(x, g);
    }

    public static double cut_norm(List<List<Double>> A_input) {
        List<List<Double>> A = deepCopyMatrix(A_input);
        int n1 = numRows(A);
        List<Double> A_col = colSum(A);
        List<Double> A_row = rowSum(A);
        double A_tot = matrixSum(A);
        List<Double> negA_row = new ArrayList<>();
        for (int i = 0; i < n1; i++){
            negA_row.add(-A_row.get(i));
        }
        appendColumn(A, negA_row);
        List<Double> newRow = new ArrayList<>();
        for (double val : A_col) {
            newRow.add(-val);
        }
        newRow.add(A_tot);
        appendRow(A, newRow);
        int p_val = Math.max(Math.min((int)Math.round(Math.sqrt(2 * n1) / 2.0), 100), 1);
        int n2 = 2 * n1 + 2;
        List<List<Double>> x0 = randomNormalMatrix(p_val, n2);
        normalizeColumns(x0);
        Pair<List<List<Double>>, List<List<Double>>> optRes = optimize(x0, A, 1e-8, 1e-10, 1e-8, 1e-4, 0.1, 0.85, 1e-3, 5, 600);
        List<List<Double>> x = optRes.first;
        int half = n2 / 2;
        List<List<Double>> U = new ArrayList<>();
        List<List<Double>> V = new ArrayList<>();
        for (int i = 0; i < p_val; i++){
            ArrayList<Double> rowU = new ArrayList<>();
            ArrayList<Double> rowV = new ArrayList<>();
            for (int j = 0; j < half; j++){
                rowU.add(x.get(i).get(j));
                rowV.add(x.get(i).get(j + half));
            }
            U.add(rowU);
            V.add(rowV);
        }
        List<List<Double>> U_trans = transpose(U);
        List<List<Double>> prod = matMultiply(U_trans, V);
        double sdp = 0.0;
        int r = numRows(A);
        int c = numCols(A);
        int prodRows = numRows(prod);
        int prodCols = numCols(prod);
        for (int i = 0; i < r && i < prodRows; i++){
            for (int j = 0; j < c && j < prodCols; j++){
                sdp += A.get(i).get(j) * prod.get(i).get(j);
            }
        }
        double sdp_result = Math.abs(sdp) / 4.0;
        return sdp_result;
    }

    public static double cut_distance(List<List<Double>> G, List<List<Double>> H) {
        int n = numRows(G);
        List<List<Double>> diff = subtractMatrix(G, H);
        List<List<Double>> A = divideMatrix(diff, (double)(n * n));
        return cut_norm(A);
    }

    public static List<List<Double>> erdosRenyiGraph(int n, double p) {
        List<List<Double>> A = new ArrayList<>();
        for (int i = 0; i < n; i++){
            ArrayList<Double> row = new ArrayList<>();
            for (int j = 0; j < n; j++){
                row.add(0.0);
            }
            A.add(row);
        }
        for (int i = 0; i < n; i++){
            for (int j = i + 1; j < n; j++){
                if (rand.nextDouble() < p) {
                    A.get(i).set(j, 1.0);
                    A.get(j).set(i, 1.0);
                }
            }
        }
        return A;
    }
    
    public static List<List<Double>> deepCopyMatrix(List<List<Double>> M) {
        List<List<Double>> copy = new ArrayList<>();
        for (List<Double> row : M) {
            ArrayList<Double> newRow = new ArrayList<>(row);
            copy.add(newRow);
        }
        return copy;
    }
    
    public static void main(String[] args) {
        int n = 30;
        double p1 = 0.2, p2 = 0.5;
        List<List<Double>> A = erdosRenyiGraph(n, p1);
        List<List<Double>> B = erdosRenyiGraph(n, p2);
        
        double s = cut_distance(A, B);
        System.out.println("The cut distance between A and B: " + s);
    }
}

