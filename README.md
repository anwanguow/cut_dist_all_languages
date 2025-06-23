The Implementation of Cut Distance in All Languages
==============

This repository contains the demonstration of the cut distance in all programming languages, including C, C++, Java, Python, Fortran, Matlab, and Mathematica (both imperative and functional programming).

The implementation of cut distance is based on computing the cut norm via a semidefinite programming (SDP) relaxation on the Stiefel manifold. As the number of nodes increases, the computed results tend to converge more closely to their true values, due to reduced randomness and better structural resolution—reflecting the asymptotic nature of cut distance in graph limit theory. The optimization algorithm relies on a random seed, which in the code is dynamically assigned based on the current time. To ensure consistent results and deterministic output across runs, the random seed should be fixed to a constant value by users. Note that the original method is proposed by Zaiwen Wen and Wotao Yin, see https://link.springer.com/article/10.1007/s10107-012-0584-1.

Algorithm Implementation:
-----------------

For a matrix $A \in \mathbb{R}^{N \times N}$, its $(\infty \to 1)$-norm is defined as:

$$
|A|_{\infty \to 1} = x_{X \in M_{N,p}},
$$

where $M_{N,p} = \{ X \in \mathbb{R}^{N \times p} : X^T X = I \}$ is called Stiefel manifold, and the following is the Frobenius inner product:

$$
\langle A, X\rangle_F = \sum_{i,j}A_{ij}X_{ij}, \quad X\in \mathbb{R}^{N\times p}.
$$

Note that when $p=1$ one recovers the usual $\ell_\infty\!\to\!\ell_1$ operator norm:

$$
\|A\|_{\infty\to1}
  = \max_{\|x\|_\infty\le1}\|Ax\|_1
  = \max_{\|x\|_\infty=1,\;\|y\|_\infty=1} x^T A y,
$$

where $x,y\in\mathbb{R}^N$ are the original symbolic $\{\pm1\}$ vectors. To obtain the orthogonal relaxation, we replace each scalar sign variable $x_i,y_j\in\{\pm1\}$ by a unit vector $u_i,v_j\in\mathbb{R}^d$ with $\|u_i\|_2=\|v_j\|_2=1$, and then stack $\{u_i\}$ into $U\in\mathbb{R}^{d\times N}$, and stack $\{u_i\}$ into $V\in\mathbb{R}^{d\times N}$.

This transforms the discrete optimization over $\{\pm1\}$ into a continuous problem on the Stiefel manifold $M_{N,p}$. 

The Grothendieck inequality offers a upper bound of the cut norm of a matrix, converting the calculation into a so-called ($\infty \rightarrow 1$)- norm. In detail, for a matrix $A \in \mathbb{R}^{N \times N}$:

$$
|A|_{\square} \leq |A|_{\infty \rightarrow 1} \leq K_G |A|_{\square},
$$

where $K_G$ is the Grothendieck constant.

The operation for finding a trace can be written as the form of a matrix multiplication by the inner product of two vectors, as follow:

$$
 \langle A, X\rangle_F
  = \sum_{i=1}^N\sum_{j=1}^p A_{ij}\,X_{ij}.
$$

Of course, when $p=1$ (so $X=xy^T$, $\|x\|_\infty=\|y\|_\infty=1$), we have:

$$
\langle A,\,X\rangle_F
  = \sum_{i,j}A_{ij}x_i y_j
  = x^T A y.
$$

In fact, the Grothendieck constant $K_G$ is a function of the $d$, that is, $K_G(d)$, where $d$ is the dimension of the real Hilbert space spanned by the vectors $\textbf{x}_i$ and $\textbf{y}_j$. Grothendieck proved that $K_G(d)$ is increasing to $d$ and bounded, so that the limit of $K_G(d)$ exists, as well as gave an effective estimate for $K_G(d)$.

For the estimation of $K_G$, we have:

$$
1.57 \approx \frac{\pi}{2} \leq K_G = \sup_{d} K_G(d) \leq \text{sinh} \frac{\pi}{2} \approx 2.3.
$$

In practical calculations, we take $K_G=1.782$.

Therefore, calculating the cut norm of a matrix is essentially calculating its $\infty \rightarrow 1$- norm. [N. Alon and A. Naor](https://web.math.princeton.edu/~naor/homepage%20files/cutnorm.pdf) proposed a method that convert this optimization into a problem of semidefinite programming (SDP). For the SDP solving, [Z. Wen and W. Yin](https://link.springer.com/article/10.1007/s10107-012-0584-1) designed a solver for linear searching on Stiefel manifold by using the Cayley transformation to do this kind of $|A|_{\infty \rightarrow 1}$ SDP relaxation.

Therefore, the relaxation of the cut norm of a matrix can be formulated as the following optimization:

1. The objective function is:
   $$
   \max_X \langle A,\,X\rangle_F.
   $$
2. The feasible domain is:
   $$
   M_{N,p} = \{X \in \mathbb{R}^{N \times p} \mid X^TX=I \}.
   $$

After the optimization is completed, we have:

$$
|A|_{\square} \approx \frac{|A|_{\infty \rightarrow 1}}{K_G}.
$$

However, in the realistic implementation, in order to express the objective function in the form of "taking the trace of a square matrix", we need to extend $A$ and $X$ as follows:

1. **Extending $A$ to $\mathfrak{A}$:**

   $$
   \mathfrak{A} = 
   \begin{pmatrix}
       A & -A \mathbf{1} \\
       -\mathbf{1}^T A & \mathbf{1}^T A \mathbf{1}
   \end{pmatrix}
   \in \mathbb{R}^{(N+1)\times (N+1)},
   $$
   where $\mathbf{1}$ denotes the column vector of all ones.

2. **Extending $X$ to $\mathfrak{X}$:**
   
   a. Let $\mathscr{X} \in \mathbb{R}^{p \times \mathfrak{n}}$ be a matrix to be initialized, where $\mathfrak{n} = 2(N+1)$ and $\mathscr{X}^T \mathscr{X} = I_{\mathfrak{n}}$. The initialization method for $\mathscr{X}$ will be introduced later. Then, $\mathscr{X}$ is partitioned into horizontal blocks as follows:
   
   $$
   \mathscr{X} = \begin{bmatrix} U \mid V \end{bmatrix},\quad \text{where}\; U, V \in \mathbb{R}^{p \times (N+1)}.
   $$
   
   b. In the subsequent optimization process, the Cayley update and tangent space projection are first applied to $\mathscr{X}$. When computing the objective function and gradient, we take:
   
   $$
   \mathfrak{X} = U^T V \in \mathbb{R}^{(N+1) \times (N+1)}.
   $$

After extending the matrix, we need to use the new inequality:

$$
 |\mathfrak{A}|_{\square} \leq |\mathfrak{A}|_{\infty\to 1} \leq 4|\mathfrak{A}|_{\square}.
$$

Besides, the optimization process is transformed into:

1. The objective function is:
   
   $$
   \max_{\mathscr{X}^T\mathscr{X} = I_{\mathfrak{n}}} \operatorname{Tr} (\mathfrak{A}^T \mathfrak{X}).
   $$
   
2. The feasible domain is:
   
   $$
   M_{p,\mathfrak{n}} = \left\{ \mathscr{X} \in \mathbb{R}^{p \times \mathfrak{n}} \mid \mathscr{X}^T\mathscr{X} = I_{\mathfrak{n}} \right\}.
   $$

The optimization algorithms contains two steps:

**1. Initialization:**  
The initialization is intended to provide a starting point for subsequent iterative updates. To start with, a random matrix $Z \in \mathbb{R}^{\mathfrak{n} \times p}$ is generated. Then, the QR decomposition is performed on $Z$, resulting in $Z = QR$, where $Q$ is an orthogonal matrix. The initial value of $\mathscr{X}$ is set from $Q$, i.e., $\mathscr{X}_0 = Q$. Note that this procedure corresponds precisely to relaxing each original scalar variable $x_i \in \{\pm1\}$ (and similarly $y_j$) into a continuous vector $u_i \in \mathbb{R}^d$ of unit norm, so that the combinatorial problem on $\{\pm1\}^N$ becomes optimization over the Stiefel manifold $M_{p,\mathfrak{n}}$ via the matrix $\mathscr{X}$.

For convenience, we denote $\operatorname{Tr}(\mathfrak{A}^T\mathscr{X})$ as $f(\mathscr{X})$. Considering $\nabla f(\mathscr{X}) = \mathfrak{A}$, then the gradient at the initial state can be calculated as follows (tangent space projection):

$$
G = \mathfrak{A} - \mathscr{X} \frac{\mathscr{X}^T \mathfrak{A} + \mathfrak{A}^T \mathscr{X}}{2}.
$$

Equation \eqref{equ:init_G_gradient} offers the direction for gradient update at the start of the algorithm, where the term after the minus sign ensures that the optimization process adheres to the constraints of the Stiefel manifold $M_{p,\mathfrak{n}}$.

In order to make sure the matrix $\mathscr{X}$ always satisfies $\mathscr{X}^T \mathscr{X} = I$, the Cayley transformation is used to implement the process of assignment. Let $G = \mathfrak{A} - \mathscr{X} \frac{\mathscr{X}^T \mathfrak{A} + \mathfrak{A}^T \mathscr{X}}{2}$, and $W = G \mathscr{X}^T - \mathscr{X} G^T$, so that $W^T = -W$ and $W$ is skew-symmetric. Then the update is as follows:

$$
C(\tau) := \left(I - \frac{\tau}{2} W\right)^{-1} \left(I + \frac{\tau}{2} W\right), \quad
Y(\tau) = C(\tau)\mathscr{X}.
$$

Through this process, a new trial matrix $Y(\tau)$ can be generated, which is then used to decide whether to replace $\mathscr{X}$ with $Y(\tau)$ according to the following Armijo-Wolfe condition:

$$
\left\{
\begin{aligned}
    &f(\mathscr{X}+\tau \Delta \mathscr{X}) \leq f(\mathscr{X}) + c_1 \tau \nabla f(\mathscr{X})^T \Delta \mathscr{X},\\
    &\nabla f(\mathscr{X}+\tau \Delta \mathscr{X})^T \Delta \mathscr{X} \geq  c_2 \nabla f(\mathscr{X})^T \Delta \mathscr{X},
\end{aligned}
\right.
$$

where $c_1 \in [0,0.5]$ and $c_2 \in [c_1, 1]$.

Specifically, the first inequality, known as the Armijo condition, ensures that each update achieves a sufficient decrease in the gradient. The second inequality, called the Wolfe condition, keeps the step length aligned with a reasonable gradient descent direction.

If this condition is satisfied, $Y(\tau)$ will be used to update $\mathscr{X}$. Otherwise, $\mathscr{X}$ will remain unchanged, and the value of $\tau$ will be adjusted as follows:

$$
\tau = \frac{s^T y}{y^T y} \quad \text{or} \quad \tau = \frac{s^T s}{s^T y},
$$

where $s = \mathscr{X}_k - \mathscr{X}_{k-1}$ and $y = \nabla f(\mathscr{X}_k) - \nabla f(\mathscr{X}_{k-1})$. Note that, the former of these two formulas tends to keep the step length stable, while the latter tends to speed up the update process.

Then, $Y(\tau)$ needs to be recalculated according to above assignment, and a suitable $\tau$ should be found to meet the Armijo-Wolfe condition. When the first assignment from $Y(\tau)$ to $\mathscr{X}$ occurs, the initialization is complete. At this point, the algorithm will enter an iterative loop.


**2. Iteration:**  
Each step of this iteration generates a new trial matrix $Y(\tau)$ based on the current $\mathscr{X}$ and uses the Armijo-Wolfe condition to evaluate whether to use $Y(\tau)$ to update $\mathscr{X}$:

1. Calculating the new direction for gradient update:
   
   $$
   W = \nabla f(\mathscr{X}) - \mathscr{X} \nabla f(\mathscr{X})^T;
   $$

2. Generating the new trial matrix:
   
   $$
   Y(\tau) = (I - \frac{\tau}{2}W)^{-1}(I + \frac{\tau}{2}W)\mathscr{X};
   $$

3. Assigning $Y(\tau)$ to $\mathscr{X}$ if the Armijo-Wolfe condition meets; otherwise, adjusting $\tau$ and recalculating $Y(\tau)$ until the Armijo-Wolfe condition is satisfied.

The iteration ends when the redefined maximum number of steps is reached, or the gradient update $|\nabla f(\mathscr{X})|$ is less than the predefined threshold (tolerance $\varepsilon$).

The time complexity of the constraint-preserving update scheme is $8Np^2 + O(p^3)$. By selecting $p$ as the column number of unextended $X$, by:

$$
p = \max\left(\min\left(\operatorname{round}\left(\sqrt{2N}/2\right),\, 100\right),\, 1\right),
$$

the computational complexity can be reduced to $4N^2 + O(N^{3/2})$.

Here, $p$ refers to the dimension of the identity matrix $I = X^T X$ used to maintain the orthogonality constraint during the original $|A|_{\infty \rightarrow 1}$ SDP relaxation before the extension, where the unextended $X \in \mathbb{R}^{N \times p}$ [Wen et al., 2013]. Because the method utilizes gradient methods on Stiefel manifolds combined with line search, achieving linear convergence to the optimum, the algorithm will converge in $O(\varepsilon^{-1})$ iterations, where $\varepsilon$ is the tolerance [Wen et al., 2013]. Therefore, the overall computational cost for cut distance is $O(N/\varepsilon)$, which simplifies to $O(N^2)$ with a constant tolerance $\varepsilon$ [[Z. Wen and W. Yin, 2013](https://link.springer.com/article/10.1007/s10107-012-0584-1)].

After completing the process of optimization, we have:

$$
|A|_{\square} \approx \frac{|\mathfrak{A}|_{\infty \rightarrow 1}}{4}.
$$

The relationship between these variables and the objective is:

$$
\mathrm{d}_{\square}(G, H) = |M^G - M^H|_{\square} = |A|_{\square} \approx \frac{|A|_{\infty \rightarrow 1}}{K_G} \approx \frac{|\mathfrak{A}|_{\infty \rightarrow 1}}{4}.
$$



Related Article:
-----------------
The related article "Graph theory based approach to identify phase transitions in condensed matter" is published in Physical Review B (PRB) and can be accessed at https://journals.aps.org/prb/abstract/10.1103/PhysRevB.111.054116.

In this article, cut distance is used for constructing two kinds of pseudo-order parameters, i.e., structural order parameter $\mathcal{D}_s$ and dynamical order parameter $\mathcal{D}_d$.

The repo of this article is https://github.com/anwanguow/graph_phase_transition.

Contact:
-----------------
An Wang: amturing@outlook.com


