# Mathematical Computation I: Matrix Computation Course 
## UChicago STAT 30900, Autumn 2023 

\newcommand{\N}{\mathbb{N}}
\newcommand{\R}{\mathbb{R}}
\newcommand{\Q}{\mathbb{Q}}
\newcommand{\Z}{\mathbb{Z}}
\newcommand{\C}{\mathbb{C}}

\DeclareMathOperator{\Jac}{Jac}
\DeclareMathOperator{\Ker}{Ker}
\DeclareMathOperator{\Spec}{Spec}
\DeclareMathOperator{\Tr}{Tr}
\DeclareMathOperator{\diag}{diag}
\DeclareMathOperator{\argmin}{argmin}

\let\temp\phi
\let\phi\varphi
\let\varphi\temp

\newcommand{\pa}[1]{\left(#1\right)}
\newcommand{\bra}[1]{\left[#1\right]}
\newcommand{\cbra}[1]{\left\{#1\right\}}
\newcommand{\norm}[1]{\|#1\|}
\newcommand{\lrnorm}[1]{\left\|#1\right\|}

\newcommand{\mat}[1]{\begin{matrix}#1\end{matrix}}
\newcommand{\pmat}[1]{\pa{\mat{#1}}}
\newcommand{\bmat}[1]{\bra{\mat{#1}}}

\newcommand{\pfrac}[2]{\pa{\frac{#1}{#2}}}
\newcommand{\bfrac}[2]{\bra{\frac{#1}{#2}}}
\newcommand{\psfrac}[2]{\pa{\sfrac{#1}{#2}}}
\newcommand{\bsfrac}[2]{\bra{\sfrac{#1}{#2}}}

### Linear Algebra Review 

By convention, vectors are column vectors. 

Set $V$ a vector space (almost always real or complex).

_Def_: $\| \cdot \| : V \to \mathbb R$ is a **norm** if it satisfies the following properties:

- $\|v\| \geq 0$ for any $v \in V$,
- $\|v\| = 0$ iff $v = 0$,
- $\|\alpha v\| = |\alpha| \|v\|$ for $\alpha \in \mathbb R$ and $v \in V$,
- $\|v + w\| \leq \|v\| + \|w\|$ for any $v, w \in V$.

Now, since this is a computational class, we only care about specific norms, almost all of which we can quickly qrite down. 

#### Vector Norms 

Set $V = \mathbb R^n$ or $\mathbb C^n$ equivalently.

_Def_: The **Minkowski** or $p$*-norm* is given by
$$
  \| x \|_p = \left(\sum_{i=1} x_i^p\right)^{1/p}
$$
and we call the $2$-norm the **Euclidean norm** and the $1$-norm the **Manhattan norm**.

_Def_: The $\infty$*-norm* is the limit of $p$-norms as $p \to \infty$, and is given by
$$
  \| x \|_\infty = \max_{i = 1, \dots, n} |x_i| = \lim_{p \to \infty} \| x \|_p.
$$

_Def_: For a weight vector $\underline w = \begin{bmatrix}w_1, \dots, w_n\end{bmatrix}^T \in \mathbb R^n$, with each $w_i > 0$, we have that the **weighted** $p$*-norm* is
$$
  \| v \|_{\underline w, p} = \left(\sum_{i=1} w_i x_i^p\right)^{1/p}.
$$

_Def_: In general, for any positive definite matrix $A$ (that is, $x^TAx > 0$ for all $x \neq 0$), we may consider the **Mahalanobis norm**
$$
  \| v \|_{A} = \left(x^T A x\right)^{1/2}.
$$

As convention, the "default" norm when a subscript is omitted is the Euclidean norm.

#### Matrix Norms 

Set $V = \mathbb R^{m \times n}$ or $\mathbb C^{m \times n}$ equivalently.

_Def_: The **Hölder** $p$*-norms* are given by
$$
  \| X \|_{H, p} = \left(\sum_{i=1}^m\sum_{j=1}^m |x_{ij}|^p \right)^{1/p},
$$
and the Hölder $2$-norm is called the Frobenius norm, which is also defined on infinite dimensional vector spaces as
$$
  \| X \|_F = \left(\Tr(XX^*)\right)^{1/2}
$$
where $^*$ is the conjugate transpose.

_Def_: As before, we can take $p \to \infty$ to get the **Hölder** $\infty$*-norm* given by
$$
  \| X\|_{H, \infty} = \max_{\substack{i = 1, \dots n \\ j = 1, \dots, n}} |x_{ij}|.
$$

_Def_: We can also define norms on matrices by viewing them as linear maps $A: \mathbb R^n \to \mathbb R^m$; in particular, if we have some norm $\| \cdot \|_a$ on $\mathbb R^n$ and some norm $\| \cdot \|_b$ on $\mathbb R^m$, we may define the **operator norm** (or **induced norm**)
$$
  \| A \|_{a, b} = \max_{x \neq 0} \frac{\| A x \|_b}{\| x \|_a}.
$$
In particular, if the norms on the domain and codomain are just $p$-norms, we write
$$
  \| A \|_{p} = \max_{x \neq 0}\frac{\| A x\|_p}{\| x \|_p}
$$
and call it the $p$*-norm* of $A$. In particular, we call the $2$-norm the spectral norm. Further, the $1$-norm and $\infty$-norm are just
$$
  \| A \|_1 = \max_{j = 1, \dots, n} \left(\sum_{i=1}^m |a_{ij}| \right),
$$
which is the max column sum, and
$$
  \| A \|_\infty = \max_{i = 1, \dots, m} \left(\sum_{j=1}^n |a_{ij}| \right),
$$
which is the max row sum; both facts are easy to check.

In general, for $p \notin \{1, 2, \infty\}$, computing $\| A \|_p$ is NP-hard, and if we consider $\| A \|_{p,q}$ then $\|A\|_{\infty, 1}$ is hard and $\|A\|_{1, \infty}$ is easy.

#### Properties of Norms 

We may also want to consider some other desirable properties on our norms.

- For example, we might want **submultiplicativity**:
$$
  \| A B\| \leq \|A \| \| B \|.
$$
The Frobenius norm is submultiplicative.
- Take also **consistency**:
$$
  \| Ax \|_b \leq \| A \|_{a,b} \|x\|_a.
$$
This is true for $p$-norms, but not in general.

Some properties always hold.

**Prop**: Every norm is Lipschitz.

_Proof_: Let our norm be $\| \cdot \| : V \to \R$. The triangle inequality immediately implies
$$
  | \| u \| - \| v \| | \leq \| u - v \|.
$$

**Theorem (Equivalence of Norms)**: [Link](https://kconrad.math.uconn.edu/blurbs/gradnumthy/equivnorms.pdf) Set $V$ a finite-dimensional vector space. Then every pair of norms $\| \cdot \|_a, \| \cdot \|_b$ are equivalent to each other, e.g. there are constants $c_1, c_2$ such that for any $v \in V$,
$$
  c_1\| v \|_b \leq \| v \|_a \leq c_2 \| v \|_b.
$$

_Proof_: Induct on the dimension of $V$ and see that every norm is equivalent to the infinity norm.

_Def_: We say that a sequence $\{ x_k \}_{k=1}^\infty$ of vectors **converges** to $x$ if
$$
  \lim_{k \to \infty} \| x_k - x \| = 0.
$$

Then, the above clearly shows that convergence in one norm implies convergence in every norm.

#### Inner, Outer, Matrix Products 

_Def_: Set $V$ a $K$-vector space (where $K = \R, \C$). An **inner product** is a binary operation $\langle \cdot, \cdot \rangle: V \times V \to \R$ which satisfies that

- $\left\langle v, v \right\rangle \geq 0$ for all $v \in V$,
- $\left\langle v, v,\right\rangle = 0$ if and only if $v = 0$,
- $\left\langle u, v \right\rangle = \overline{ \left\langle v, u\right\rangle }$ for all $u, v \in V$,
- $\left\langle \alpha_1 u_1 + \alpha_2 u_2, v\right\rangle = \alpha \left\langle u_1, v \right\rangle + \alpha_2 \left\langle u_2, v\right\rangle$ for all $u_1, u_2, v \in V, \alpha_1, \alpha_2 \in K$.

**Prop**: For an inner product $\left\langle \cdot, \cdot \right\rangle$,
$$
  \| v \| = \sqrt{ \left\langle v, v \right\rangle }
$$
is a norm. Furthermore, an arbitrary norm $\| \|$ is induced by an inner product if and only if it satisfies the parallelogram law
$$
  \norm{u + v}^2 + \norm{u - v}^2 = 2 \norm u ^ 2 + 2 \norm v ^2.
$$

**Theorem (Cauchy-Schwarz)**: [Link](https://en.wikipedia.org/wiki/Cauchy%E2%80%93Schwarz_inequality) Let $\norm \cdot$ be induced by $\left\langle \cdot, \cdot \right\rangle$. Then
$$
  \sqrt{\left\langle u, v \right\rangle} \leq \norm u \norm v.
$$


_Def_: The **standard Euclidean inner product** on $\C^n$ is
$$
  \left\langle x, y\right\rangle = x^*y.
$$

_Def_: The **Frobenius** norm on $\C^{m \times n}$ is 
$$
  \left\langle X, Y\right\rangle = \sum_{i=1}^m \sum_{j=1}^n x_{ij} y_{ij} = \Tr(X^*Y).
$$

**Theorem (Hölder Inequality)**: [Link](https://en.wikipedia.org/wiki/H%C3%B6lder%27s_inequality) For $x, y \in \C^n$ and $p^{-1} + q^{-1} = 1$, we have
$$
  |x^*y| \leq \norm x_p \norm y_q.
$$

**Theorem (Bessel Inequality)**: [Link](https://en.wikipedia.org/wiki/Bessel%27s_inequality) For $x \in \C^n$ and an orthonormal basis $e_1, \dots, e_n$, we have 
$$
  \sum_{k=1}^n \left| \left\langle x, e_k \right\rangle^2\right|  \leq \norm x_2.
$$

_Def_: The **outer product** is a binary operator $\C^m \times \C^n \to \C^{m \times n}$ taking 
$$
  (x, y) \mapsto xy^*.
$$

**Prop**: $A \in \R^{m \times n}$ is an outer product iff and only if it has rank 1.

_Def_: The matrix product is a binary operator $\C^{m \times n} \times \C^{n \times p} \to \C^{M \times p}$. Set $A = \bmat{\alpha_1 & \alpha_2 & \cdots & \alpha_m}^T$ and $B = \bmat{\beta_1 & \beta_2 & \cdots & \beta_p}$. Then,
$$
  AB = \bmat{
    \alpha_1^T \beta_1 & \cdots & \alpha_1^T \beta_n \\
    \alpha_2^T \beta_1 & \cdots & \alpha_2^T \beta_n \\
    \vdots & \ddots & \vdots \\
    \alpha_m^T \beta_1 & \cdots & \alpha_m^T \beta_n 
  } 
  = \bmat{A\beta_1 & A\beta_2 & \cdots & A\beta_n} 
  = \bmat{\alpha_1^TB \\ \alpha_2^TB \\ \vdots \\ \alpha_m^TB}.
$$
Alternatively, it is uniquely characterized as the matrix representing the composition of $A$ and $B$ as operators.

**Prop**: Let $D = \diag(d_1, \dots, d_n)$, the diagonal matrix with entries $d_1, \dots, d_n$; then
$$
  AD = \bmat{d_1\alpha_1, \dots, d_n \alpha_n}.
$$
Simiar for the other direction of multiplication.

### Eigenvalues and Eigenvectors 

_Def_: For a complex matrix $A$, an **eigenvalue** $\lambda \in \C$ and **eigenvector** $v \neq 0$ satisfy
$$
  Av = \lambda v.
$$ 

_Def_: Furthermore, an **eigenspace** is the span of all eigenvectors correspnding to a single eigenvalue, the **spectrum** of $A$, $\Spec(A)$ is the set of all eigenvalues of $A$, and the **spectral radius** is $\rho(A) = \max_{\lambda \in \Spec(A)} |\lambda| = |\lambda_{\text{max}}|$. Sometimes we will call this the top eigenvector/eigenvalue.

As a convention, usually we implcitly order eigenvalues, e.g. $|\lambda_1| \geq |\lambda_2| \geq \cdots \geq \lambda_n$.

_Def_: More generally, if $v$ is a eigenvector of $A^T$, we say that $v$ is a **left eigenvector** of $A$ (and thus usual eigenvectors are **right eigenvectors**).

There are a few warnings about these things:

- real matrices may have complex eigenvectors;
- we normally normalize eigenvectors to have unit length;
- left and right eigenvectors are usually not the same.

_Def_: A square matrix $A \in \C^{n \times n}$ is diagonalizable if it is simiarly to a diagonal matrix.

**Prop**: $A$ is diagonalizable if and only if it has $n$ linearly independent eigenvectors.

_Proof_: The eigenvectors form a basis by intertibility; change basis from the eigenvectors, scale by eigenvalues, and change basis back. In particular, let $A \in \C^{n \times n}$ have eigenvalues $\lambda_1, \dots, \lambda_n$, with corresponding eigenvectors $v_1, \dots, v_n$; set $\Lambda = \diag(\lambda_1, \dots, \lambda_n)$ and $V = \bmat{v_1, \dots, v_n}$. Then
$$
  A = V \Lambda V^{-1}.
$$

_Def_: The above $A = X \Lambda X^{-1}$ is called the **eigenvalue decomposition**, or **EVD**.

**Prop**: The diagonal entries in a triangular matrix are just the eigenvalues.

_Def_: A matrix $A \in \C^{n \times n}$ is **normal** if it commutes with its adjoint $A^*$ (the conjugate transpose), and **unitary** if $AA^* = A^*A = I$.

**Prop**: $A$ is unitary if and only if it's columns (or rows) are orthonormal.

**Theorem**: $A \in \C^{n \times n}$ is normal if and only if it is unitarily diagonalizable, e.g. if $A = V\Lambda V^*$ with $V$ unitary.

_Def_: $A \in \C^{n \times n}$ is **Hermitian** if $A^* = A$.

**Theorem (Spectral Theorem)**: [Link](https://en.wikipedia.org/wiki/Spectral_theorem ) $A \in \C^{n \times n}$ is Hermetian if and only if $A$ is unitarily diagonalizable with all real eigenvalues.

**Corollary**: $A \in \R^{n \times n}$ is symmetric if and only if it is orthgonally diagonalizable with all eigenvalues real.

#### Jordan Canonical Form 

_Def_: [ Link](https://en.wikipedia.org/wiki/Jordan_normal_form) Any matrix can be written in **Jordan canonical form**, e.g.
$$
  A = XJX^{-1}
$$
where $J$ is nonzero only on the diagonal and the superdiagonal, which has values only in $\{0, 1\}$, such as
$$
  J = \begin{bmatrix}
    \lambda_1 & 1         &           &   & & & \\
              & \lambda_1 & 1         &   & & & \\
              &           & \lambda_1 & 0 & & & \\
              &           & & \lambda_2 & 0 & & \\
              &           &         &  & \lambda_3 & 1 & \\
              &           &         &  & & \lambda_3  & \\
  \end{bmatrix}
$$
for example.

The way it's written above is not by coincidence: you can permute everything so that $J$ is composed of Jordan blocks, which have the same entry all the way down the diagonal and have a superdiagonal of only ones, e.g.
$$
  J = \begin{bmatrix}
    J_1 & & & \\
    & J_2 & & \\
    & & \ddots & \\
    & & & J_k \\
  \end{bmatrix}
$$
where each $J_i = \lambda I + N$, where $N$ has all zero entries except on the superdiagonal, on which it is always one.

Unfortunately, the JCF is pretty useless in application.

**Theorem (Golub-Wilkinson)**: The Jordan canonical form cannot be computed in finite precision.

#### Spectral Radius 

Let's return to the spectral radius,
$$
  \rho(A) = \max_{\lambda \in \Spec(A)} |\lambda|.
$$

Any nilpotent matrix shows that $\rho$ is not a norm, but the following is true.

**Prop**: If $\| \cdot \|: \C^{m \times m} \to \R$ is any consistent matrix norm, then $\rho(A) \leq \| A \|$.

_Proof_: Look at the norm of the image of a top eigenvector.

**Theorem**: Given any $A \in \C^{n \times n}$, any positive $\epsilon$, there is a consistent norm (in fact, an operator norm) $\| \cdot \|: \C^{n \times m} \to \R$ such that 
$$
  \| A \| \leq \rho(A) + \epsilon.
$$

**Prop**: Given any matrix norm $\| \cdot \|: \C^{n \times m} \to R$, we have
$$
  \rho(A) = \lim_{k \to \infty} \| A^k \|^{1/k}.
$$

**Lemma**: $\lim_{k \to \infty} A^k = 0$  if and only if $\rho(A) < 1$.

_Proof_: $(\implies)$ Set $\lambda$ to be a top eigenvalue of $A$, and $x$ a corresponding eigenvector. Then, $A^k x = \lambda^k x$. Send $k \to \infty$ and conclude.

$(\impliedby)$ By the above theorems, there is some operator norm $\| \cdot \|$ such that $\|A\| \leq \rho(A) + \epsilon < 1$. Send $k \to \infty$ and use the fact that operator norms imply $\| A^k \| \leq \|A\|^k$ and win.

#### Finding Eigenvalues 

**Theorem (Gershgorin Circle Theorem)**: [Wikipedia](https://en.wikipedia.org/wiki/Gershgorin_circle_theorem) Let $A \in \C^{n \times n}$, with entries $a_{ij}$ and for $1 \leq i \leq n$ set $r_i = \sum_{j\neq i} |a_{ij}|$. Then, every eigenvalue of $A$ lies within the union of the Gershgorin discs
$$
  G_i = \left\{ z \mid z \in \C, |z - a_{ii}| \leq \sum_{j \neq i} |a_{ij}| \right\}
$$
and the number of eigenvalues in each connected component is equal to the number of Gershgorin disks that constitute that component.

_Proof_: If $A \in \C^{n \times n}$ is strictly diagonally dominant, e.g.
$$
  |a_{ii}| > \sum_{j \neq i}|a_{ij}|
$$
for all $1 \leq i \leq n$, then $A$ is invertible. To see this, take any $x \in \ker(A)$ so that $\sum_{j=1}^n a_{ij}x_i = 0$. But look at the index $k$ that witnesses $\| x \|_\infty$:
$$
  -a_{kk}x_k = \sum_{j \neq k} a_{kj}x_j \implies |a_{kk}||x_k| \leq \sum_{j \neq k} |a_{kj}||x_k|
$$
so $x = 0$ and $\ker(A)$ is trivial.

We proceed to prove that any $z \notin \bigcup_{i=1}^n G_i$ cannot be an eigenvalue by showing that $A - zI$ is invertible. But this is clear, since $A - zI$ in that case is strictly diagonally dominant.

_Def_: [Wikipedia](https://en.wikipedia.org/wiki/Schur_decomposition) The **Schur decomposition** of a matrix $A$ is $A = QUQ^*$ such that $Q$ is unitary and $U$, the **Schur form** of $A$ is upper triangular.

Note that you can read off the eigenvalues of $A$ from the diagonal of the Schur form, which is numerically stable to compute.

As an aside, how one numerically finds the roots of a polynomial $p(x) = \sum_{k=0}^m c_kx^k$ is by forming the companion matrix
$$
  A = \begin{bmatrix}
    0 & \cdots & 0 & -c_0 \\
    1 & \ddots &0 & -c_1 \\
    \vdots & \ddots & \vdots & \vdots \\
    0 & \cdots & 1 & -c_{d-1} 
  \end{bmatrix}.
$$

### Singular Value Decomposition 

_Def_: [Link](https://en.wikipedia.org/wiki/Singular_value_decomposition) The **SVD decomposition** of a real (resp. complex) matrix $A \in \C^{m \times n}$) is
$$
  A = U \Sigma V^*
$$
where $U, V$ are orthogonal (resp. unitary) and $\Sigma$ is diagonal (on the shorter diagonal) with nonnegative entries; that is, if $m > n$, then
$$
  \Sigma = \bmat{\diag(\sigma_1, \dots, \sigma_n) \\ \mathbf 0} \in \R^{m \times n}_{\geq 0}
$$
and if $m < n$,
$$
  \Sigma = \bmat{\diag(\sigma_1, \dots, \sigma_n) & \mathbf 0} \in \R^{m \times n}_{\geq 0}
$$

_Def_: We may arrange $\sigma_{1} \geq \sigma_2 \geq \dots \geq \sigma_{\min(m, n)} \geq 0$; we then call top few singular values the principal singular values. The columns of $U$ are left singular vectors, and the columns of $V$  are right singular vectors.

**Prop**: The largest $r$ for which $r > 0$  is the rank of $A$.

_Def_: We sometimes use the **reduced SVD**, which is the same as the normal SVD but we drop the extra rows of $\Sigma$ to force it to be a square diagonal matrix in $\C^{\min\{m, n\} \times \min\{m, n\}}$. Then,
$$
  A = \sum_{k=1}^r \sigma_k u_r v_r^*
$$
where $u_r, v_r$ are the columns of $U, V$ in the reduced SVD and $r = \operatorname{rank}(A)$, or the number of nonzero singular values.

_Def_: In a similar vein, the **condensed SVD** is attained by removing all of the nonzero rows/columns of $\Sigma$, so $\Sigma \in \C^{r \times r}$.

**Theorem**: Given any $A \in \C^{m \times n}$, there are $U \in \C^{m \times m}$, $V \in \C^{n \times n}$, and 
$$
  \Sigma = \diag(\sigma_1, \dots, \sigma_r, 0, \dots, 0) \in \R_{\geq 0}^{m \times n}
$$
such that $U, V$ are unitary and $A = U \Sigma V^*$.

_Proof_: Form the matrix
$$
  W = \bmat{0 & A \\ A^* & 0}
$$
which is Hermitian. The spectral theorem gives an EVD, and if you look carefully you get the compact SVD:
$$
  W = \bmat{ U_r & U_r \\ V_r & -V_r } \bmat{
  \sigma_1 & & & & &  &\\
   & \sigma_2 & & & & &\\
   & & \ddots & & & &\\
   & &  & \sigma_r & & &\\
   & &  & & -\sigma_1 & & \\
   & &  & & & \ddots & \\
   & &  & & & & -\sigma_r\\ }
  \bmat{ U_r & U_r \\ V_r & -V_r }^*.
$$

Singular values have a nice property: we have for any singular value $\sigma$ and left singular vector $u$ and right singular value $v$, $Av = \sigma u$ and $A^*v = \sigma v$; thus $u, v$ are eigenvalues of $AA^*$ and $A^*A$ respectively.

#### Applications 

We can read off many important quantities/spaces from the SVD. Let $A = U \Sigma V^*$ be a SVD, with only $r$ nonzero singular values.

- The rank is the number of nonzero singular values.
- The absolute value of the determinant is the product of the signular values.
- The two norm is the maximum singular value.
- Set $\sigma = (\sigma_1, \dots, \sigma_n)$; then the Frobenius norm of $A$ is $\| \sigma \|_2$, and we call the general case of this the **Schatten norm**, e.g. $\| A \|_{S, p} = \| \sigma \|_p$; the case of $p = 1$ is also called the **nuclear norm**.
- The Ky Fan $(p, k)$ norm for an integer $1 \leq k \leq \infty$ is 
  $$
    \| A \|_{S, p, k} = \left( \sum_{i=1}^p \sigma_i^p \right)^{1/p}.
  $$
- The kernel of $A = U \Sigma V^*$ is spanned by $v_{r+1}, \dots, v_n$.
- The image of $A$ is $u_1, \dots, u_r$.
- The kernel of $A^T$ (this is the cokernel) is spanned by $v_{1, \dots, v_r}$.
- The image of $A^T$ is spanned by $u_{r+1}, \dots, u_m$.

You can also solve fundamental problems. For example, if you consider the linear system $Ax = b$, it might not have exactly one solution, so we can translate this into the least squares system $\min_{x \in \R^n} \| Ax - b \|_2^2$; this still might not have a unique solution when $A$ is singular, so we consider the **minimum norm least squares** problem: $\min \{ \| x \|_2 \mid x \in \argmin \| A_x - b \|_2 \}$.

SVD solves all of these simultaneously; set $A = U \Sigma V^T$; this means that we want $\Sigma V^T x = U^T b$; equivalently, we just need to solve $\Sigma y = c$, where $b = Uc, x = Vy$. But this is trivial since $\Sigma$ is diagonal.

 We could define the pseudoinverse of $A$ to be $A^+$ to be the matrix $A^+$ that sends any $b$ to the minimum norm least squares solution of $Ax = b$.

**Prop**: $A^+$ is just $V \Sigma^{-1} U^*$, where by $\Sigma^{-1}$ is just gotten by flipping all the nonzero diagonal elements.

_Proof_: Clear.

_Def_: In general, the **pseudoinverse** or **Moore-Penrose** inverse of $A \in K^{m \times n}$ is a matrix $X \in K^{n \times m}$ that satisfies 

- $(AX)^{*} = AX$,
- $(XA)^* = XA$,
- $AXA = A$,
- $XAX = X$.

**Prop**: The above gives a unique matrix and is in fact the same as $A^+$ from earlier.

**Prop**: There following statements are true.

- $A^+A$ and $AA^+$ are not necessarily $I$.
- If $A$ has full column rank $A^+ = (A^*A)^{-1}A$.
- If $A$ has full row rank, $A^+ = A^*(AA^*)^{-1}$.
