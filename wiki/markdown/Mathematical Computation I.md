---
title: 'Mathematical Computation I: Matrix Computation Course'
subtitle: 'UChicago STAT 30900, Autumn 2023'
...

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
\DeclareMathOperator{\rank}{rank}
\DeclareMathOperator{\tr}{tr}
\DeclareMathOperator{\fl}{fl}
\DeclareMathOperator{\rcond}{rcond}

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

## Linear Algebra Review 

-------

By convention, vectors are column vectors. 

Set $V$ a vector space (almost always real or complex).

_Def_: $\| \cdot \| : V \to \mathbb R$ is a **norm** if it satisfies the following properties:

- $\|v\| \geq 0$ for any $v \in V$,
- $\|v\| = 0$ iff $v = 0$,
- $\|\alpha v\| = |\alpha| \|v\|$ for $\alpha \in \mathbb R$ and $v \in V$,
- $\|v + w\| \leq \|v\| + \|w\|$ for any $v, w \in V$.

Now, since this is a computational class, we only care about specific norms, almost all of which we can quickly write down. 

### Vector Norms 

Set $V = \mathbb R^n$ or $\mathbb C^n$ equivalently.

_Def_: The **Minkowski** or $p$**-norm** is given by
$$
  \| x \|_p = \left(\sum_{i=1} x_i^p\right)^{1/p}
$$
and we call the $2$-norm the **Euclidean norm** and the $1$-norm the **Manhattan norm**.

_Def_: The $\infty$**-norm** is the limit of $p$-norms as $p \to \infty$, and is given by
$$
  \| x \|_\infty = \max_{i = 1, \dots, n} |x_i| = \lim_{p \to \infty} \| x \|_p.
$$

_Def_: For a weight vector $\underline w = \begin{bmatrix}w_1, \dots, w_n\end{bmatrix}^T \in \mathbb R^n$, with each $w_i > 0$, we have that the **weighted** $p$**-norm** is
$$
  \| v \|_{\underline w, p} = \left(\sum_{i=1} w_i x_i^p\right)^{1/p}.
$$

_Def_: In general, for any positive definite matrix $A$ (that is, $x^TAx > 0$ for all $x \neq 0$), we may consider the **Mahalanobis norm**
$$
  \| v \|_{A} = \left(x^T A x\right)^{1/2}.
$$

As convention, the "default" norm when a subscript is omitted is the Euclidean norm.

### Matrix Norms 

Set $V = \mathbb R^{m \times n}$ or $\mathbb C^{m \times n}$ equivalently.

_Def_: The **Hölder** $p$**-norms** are given by
$$
  \| X \|_{H, p} = \left(\sum_{i=1}^m\sum_{j=1}^m |x_{ij}|^p \right)^{1/p},
$$
and the Hölder $2$-norm is called the Frobenius norm, which is also defined on infinite dimensional vector spaces as
$$
  \| X \|_F = \left(\Tr(XX^*)\right)^{1/2}
$$
where $^*$ is the conjugate transpose.

_Def_: As before, we can take $p \to \infty$ to get the **Hölder** $\infty$**-norm** given by
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
and call it the $p$**-norm** of $A$. In particular, we call the $2$-norm the spectral norm. Further, the $1$-norm and $\infty$-norm are just
$$
  \| A \|_1 = \max_{j = 1, \dots, n} \left(\sum_{i=1}^m |a_{ij}| \right),
$$
which is the max column sum, and
$$
  \| A \|_\infty = \max_{i = 1, \dots, m} \left(\sum_{j=1}^n |a_{ij}| \right),
$$
which is the max row sum; both facts are easy to check.

In general, for $p \notin \{1, 2, \infty\}$, computing $\| A \|_p$ is NP-hard, and if we consider $\| A \|_{p,q}$ then $\|A\|_{\infty, 1}$ is hard and $\|A\|_{1, \infty}$ is easy.

### Properties of Norms 

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

### Inner, Outer, Matrix Products 

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

## Eigenvalues and Eigenvectors 

-------

_Def_: For a complex matrix $A$, an **eigenvalue** $\lambda \in \C$ and **eigenvector** $v \neq 0$ satisfy
$$
  Av = \lambda v.
$$ 

_Def_: Furthermore, an **eigenspace** is the span of all eigenvectors correspnding to a single eigenvalue, the **spectrum** of $A$, $\Spec(A)$ is the set of all eigenvalues of $A$, and the **spectral radius** is $\rho(A) = \max_{\lambda \in \Spec(A)} |\lambda| = |\lambda_{\text{max}}|$. Sometimes we will call this the top eigenvector/eigenvalue.

As a convention, usually we implicitly order eigenvalues, e.g. $|\lambda_1| \geq |\lambda_2| \geq \cdots \geq \lambda_n$.

_Def_: More generally, if $v$ is a eigenvector of $A^T$, we say that $v$ is a **left eigenvector** of $A$ (and thus usual eigenvectors are **right eigenvectors**).

There are a few warnings about these things:

- real matrices may have complex eigenvectors;
- we normally normalize eigenvectors to have unit length;
- left and right eigenvectors are usually not the same.

_Def_: A square matrix $A \in \C^{n \times n}$ is diagonalizable if it is similarly to a diagonal matrix.

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

**Theorem (Spectral Theorem)**: [Link](https://en.wikipedia.org/wiki/Spectral_theorem ) $A \in \C^{n \times n}$ is Hermitian if and only if $A$ is unitarily diagonalizable with all real eigenvalues.

**Corollary**: $A \in \R^{n \times n}$ is symmetric if and only if it is orthogonally diagonalizable with all eigenvalues real.

### Jordan Canonical Form 

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

### Spectral Radius 

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

### Finding Eigenvalues 

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

## Singular Value Decomposition 

-------

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

Singular values have a nice property: we have for any singular value $\sigma$ and left singular vector $u$ and right singular vector $v$, $Av = \sigma u$ and $A^*v = \sigma v$; thus $u, v$ are eigenvalues of $AA^*$ and $A^*A$ respectively.

## Applications of the SVD

----

We can read off many important quantities/spaces from the SVD. Let $A = U \Sigma V^*$ be a SVD, with only $r$ nonzero singular values.

- The rank is the number of nonzero singular values.
- The absolute value of the determinant is the product of the singular values.
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

**Prop**: The following statements are true.

- $A^+A$ and $AA^+$ are not necessarily $I$.
- If $A$ has full column rank $A^+ = (A^*A)^{-1}A$.
- If $A$ has full row rank, $A^+ = A^*(AA^*)^{-1}$.

### Ridge Regressions

Suppose that we are looking for
$$
  \min \|Ax - b\|^2_2 \text{ such that } \|x\|^2_2 \leq \alpha^2
$$
which is a ridge regression. Remember that the pseudoinverse $A^+$ gives that $A^+b$ is the minimal norm solution; so if $\|A^+b\| \leq \alpha$, this works. On the other hand, if we have $\|A^+b\| > \alpha$, we first use the fact that the above optimization reduces to one on the boundary $\|x\|_2 = \alpha$. Then take the SVD decomposition $A = U\Sigma V^T$, and look at the Lagrange multiplier condition
$$
  A^TAx - A^Tb - \lambda x = 0 \implies x = (A^TA - \lambda I)^{-1}A^Tb
$$
and thus, setting $c = U^b$,
$$
  \begin{align*}
    \alpha^2 & = x^Tx \\
             & = b^TA(A^TA- \lambda I)^{-2}A^Tb \\
             & = c^T \Sigma(\Sigma^T\Sigma + \lambda I)^{-2}\Sigma^Tc \\
             & = \sum_{i=1}^r \frac{c_i^2\sigma_i^2}{(\sigma_i^2 + \lambda)^2} = f(\lambda).
  \end{align*}
$$ 
Then, solve $f(\lambda) = \alpha^2$ to get a solution (or really, $f(\lambda)^{-1} = \alpha^{-2}$ since most univariate root finders are bad when things go to infinity). In fact, there are many solutions, but the one we are looking for is actually the biggest solution.

### Matrix Approximation Problems

Suppose that $A \in \R^{n \times n}$; then, we want to solve
$$
  \min_{X^T = X} \|A - X\|_F
$$
or
$$
  \min_{X^TX = 1} \|A - X\|_F
$$
or if $A \in \R^{m \times n}$,
$$
  \min_{\rank(X) \leq r} \|A - X\|_F.
$$
The last problem is very important since rounding error sends every matrix to an invertible matrix, so the above approximates a solution to this problem.

Set
$$
  A = \frac{A + A^T}{2} + \frac{A - A^T}{2} = X + Y;
$$
then, $X = \frac{A + A^T}{2}$ is the solution to the first problem. Further, we can compute
$$
  \tr(X^TY) = 0 
$$
and so
$$
  \|A\|^2 = \|X + Y\|^2 = \tr(X^TX) + 2\tr(X^TY) + \tr(Y^TY) = \|X\|^2 + \|Y\|^2.
$$
In general, we see $\R^{n \times n} = S^2(\R^n) \oplus \Lambda^2(\R^n)$, e.g. the space of all matrices is the direct sum of symmetric and skew symmetric matrices.

In the second problem, take the SVD $A = U\Sigma V^T$ and $Z = U^T X V$ so that
$$
\begin{align*}
  \min_{X^TX = I}\|U \Sigma V^T - X\|_F^2 & = \min_{X^TX = I}\|\Sigma - U^TXV\|^2_F \\
                                          & = \min_{Z^TZ = I}\| \Sigma - Z \|^2_F \\
                                          & = \max_{Z^TZ = I} \| \Sigma \|_F^2 - 2 \tr(\Sigma^T Z) + \|Z\|^2_F
                                          & = \max_{Z^TZ = I} \| \Sigma \|_F^2 - 2 \tr(\Sigma^T Z) + n
\end{align*}
$$
so we just need to maximize $\tr(\Sigma^T Z) = \sum_{i=1}^n \sigma_i z_ii \leq \sum_{i=1}^n \sigma_i$; but this is attained when $Z = I$, so we need $X = UV^2$.

The solution to the last problem is called the Eckhart-Young theorem.

**Theorem (Eckhart-Young)**: Let $A \in \R^{m \times n}$ and let $A = U \Sigma V^T$ be the SVD. Then the solution to $\min_{\rank(X) \leq r} \|A - X\|_2$ is given by dropping all but the first $r$ entries in $\Sigma$.

_Proof_: Suppose not. Then let the solution be called $B$. Then
$$
  \|A - X\|_2 = \|U \diag(0, \dots, 0, \sigma_{r+1}, \dots, \sigma_{\min(m, n)}) V^T\|_2 = \sigma_{r + 1}.
$$
However, the nullity of $B$ is at least $n - r$; but if we take $W$ to be the span of the first $r+1$ right singular vectors, then any $w \in W$ is
$$
  w = V_{r+1}\alpha
$$
where $V_{r+1}$ is the first $r+1$ columns of $X$ and $\alpha$ is some vector in $\R^{r+1}$. Then,
$$
  \|Aw\|^2 = \| U\Sigma V^T V_{r+1}\alpha \|_2^2 = \sum_{i=1}^{r+1}\sigma_i^2 \|\alpha_i\|^2 \geq \sigma_{r+1}^2 \|w\|^2.
$$
But if we pick some $w \in \ker(B)$, then $Aw = (A - B)w$, so
$$
  \|Aw\|_2 \leq \|A - B\|\|w\|_2 < \sigma_{r+1}\|w\|_2
$$
and since the ranks sum above $r$, we have a contradiction.

The Eckhart-Young-Mirsky theorem extends this to all unitarily invariant norms.

### Total Least Squares / Errors in Variables Regression

Consider a linear system $Ax = b$ with $A$ having full column rank. Then, if this system is not consistent, then we can "solve it" by taking $Ax= b + r$ for some minimal $r$, which is just least squares; alternatively, we can take $(A + E)x = b$, which is the **data least squares**; if we take both, we get the **total least squares** $(A + E)x = b + r$.

In this last problem, we are minimizing $\|E\|_F^2 + \|r\|_2^2$; take $C = \bmat{A & b}$ and $F = \bmat{E & b}$. Then we can restate the constraint to $(C + F)\bmat{x \\ -1} = 0$. Then either $\rank(C) = n$, in which case $b$ is in the span of $A$ and we may take $E = 0, r = 0$, or $\rank(C) = n + 1$. In the latter case, the kernel of $C + F$ is nontrivial, and so $\rank(C + F) \leq n$.

Taking the SVD of both $C$ and $C + F$, we must have that
$$
  F = U \diag(0, \dots, 0, \sigma_{n+1}, 0, \dots, 0) V^T.
$$
Solve for $z$ in $(C + F)z = 0$, and take a solution so that the last entry is $-1$.

## Floating Point Arithmetic

_Def_: We have a few different types of errors; let $x$ be the real solution and $\hat x$ the computed solution.

- The **forward error** is $\|\hat x - x\|$, 
- the **relative error** is $\frac{\|\hat x - x\|}{\|x\|}$,
- the **pointwise error** is $\left\| \frac{x_\cdot - \hat x_{\cdot}}{x_i} \right\|$.

Error is inherent in floating point arithmetic. If there are $n$-bits, there can only be $2^n$ real numbers that are representable; the way this works is by taking each number to be of the form
$$
    \mat{ \pm & e_1 & e_2 & \dots & e_l & a_1 & a_2 & \dots & a_k} = \pm a_1.a_2a_3\dots a_k \cdot 2^{e_1e_2\dots e_l}
$$
where the $e_i$ are called the exponent and the $a_i$ are called the mantissa.

The exponent is stored as two's complement, and is computed as
$$
    e = -e_1\cdot 2^{l - 1} + \sum_{k=1}^{l} e_{k+1}\cdot 2^{l - k}.
$$

If $a_1 = 1$, it is called normal, and if $a_1 = 0$, it is called subnormal (and since they are in some ways pathological, we ignore them for now).

Double precision corresponds to $n = 64$, $l = 11$, $k = 52$. This includes numbers from $2^{-1024}$ to $2^{1024}$.

Let $F$ be the set of floating point numbers. We have some rounding scheme $\fl: \R \to F$ and some machine $\epsilon$, $\epsilon_{m}$ such that
$$
    \epsilon_m = \inf \{ x \in \R \mid x > 0, \fl(1 + x) \neq 1 \}.
$$
In double precision, $\epsilon_m = 2^{-52}$.

**Theorem (Kahan)**: For any $x \in [-2^{M}, -2^{m}] \cup [2^m, 2^M]$, there is $x' \in F$ such that $|x - x'| \leq \epsilon_m |x|$. Furthermore, for any $x, y \in F$,

- $\fl(x \pm y) = (x \pm y)(1 + \epsilon_1)$ where $|\epsilon_1| \leq \epsilon_m$,
- $\fl(xy) = xy(1 + \epsilon_2)$ where $|\epsilon_2| \leq \epsilon_m$,
- $\fl(x/y) = x/y \cdot (1 + \epsilon_3)$ where $|\epsilon_3| \leq \epsilon_m$.

Unfortunately, $F$ is terrible otherwise: it is not commutative, not associative, and not distributive. We have all types of errors here: for example, round off error is stuff like $\fl(1.1112 \times 10^5) = 1.111 \times 10^5$.

There is also overflow and underflow: when your computations exceed the defined limits of the floating point standard; also cancellation error: double precision claims that
$$
    844487^5 + 1288439^5 - 1318202^5 = 0
$$
which is about 200 billion off the real answer, but this fits within the non-significant bits of the floating point representation.

### Avoiding Floating Point Error

To compute $\|x\|_2$, you compute $\|x\|_\infty$, compute $y = \frac{x}{\|x\|_\infty}$, and then finally $\|x\|_2 = \|x\|_\infty\|y\|$. The reason to do this is to avoid underflow/overflow.

As another example, consider the sample variance;
$$
    s^2 = \frac{1}{n-1} \left( \sum_{i=1}^n x_i^2 - n^{-1} \left( \sum_{i=1}^n x_i \right)^2 \right)
$$
is a terrible formula, since it suffers from cancellation error; but
$$
    s^2 = \frac{1}{n-1} \sum_{i=1}^n (x_i - \bar x)^2
$$
is fine. Similarly, $x^2 - y^2$ is bad, but $(x + y)(x - y)$ is good. As another example, you would rather compute
$$
    x_1 = \frac{-b + \operatorname{sign}(b) \sqrt{b^2 - 4ac}}{2a}
$$
and $x_2 = \frac{c}{ax_1}$.

## Conditioning

Conditioning will be a property of the problem.

Suppose you have some problem, such as finding a solution to a linear system $Ax = b$ for some invertible $A$. One thing that we can do is perturb $A$ by a little bit, e.g.
$$
    (A + \Delta A)(x + \Delta x) = b + \Delta b;
$$
then we want to bound $\|\Delta x\| = \|\hat x - x\|$. A rough answer in this case is possible; in particular suppose $\Delta A = 0$, so that
$$
    \|\Delta x\| = \|A^{-1}\Delta b\|.
$$
Further, $\|b\| \leq \|A\|\|x\|$, so $\frac{1}{\|x\|} \leq \frac{\|A\|}{\|b\|}$ and
$$
    \frac{\|\Delta x\|}{\|x\|} \leq \|A\|\|A^{-1}\|\frac{\|\Delta b\|}{\|b\|}.
$$
The quantity $\kappa(A) = \|A\|\|A^{-1}\|$ is so important that it is called the condition number of $A$. In the case of the two norm, it is called the spectral condition number.

Now if the error is instead in $A$, we have that
$$
    \frac{\|\Delta x\|}{\|x\|} = \frac{\kappa(A)\frac{\|\Delta A\|}{\|A\|}}{1 - \kappa(A)\frac{\|\Delta A\|}{\|A\|}}.
$$
Most generally, we have
$$
    \frac{\|\Delta x\|}{\|x\|} \leq \frac{\kappa(A)\left(\frac{\|\Delta A\|}{\|A\|} + \frac{\|\Delta b\|}{\|b\|}\right)}{1 - \kappa(A)\frac{\|\Delta A\|}{\|A\|}}.
$$
Thus, if $\frac{\|\Delta A\|}{\|A\|}, \frac{\|\Delta b\|}{\|b\|} \leq \epsilon$ then
$$
    \frac{\|\Delta x\|}{\|x\|} \leq \frac{2\epsilon}{1 - \rho}\kappa(A)
$$
where $\rho = \kappa(A)\frac{\|\Delta A\|}{\|A\|}$.


_Def_: Given any norm $\|\cdot \|: \C^{m \times n} \to \R$, the **condition number of a matrix** $A \in \C^{m \times n}$ is defined by
$$
    \kappa_{\|\cdot\|}(A) =  \|A\|\|A^+\|.
$$

In the case of the spectral norm, this is just the ratio of the largest singular value to the smallest singular value.

We have that in general, $1 \leq \kappa_2(A) < \infty$, and $\kappa_2(A) = 1$ (perfectly conditioned) if and only if it is unitary.

Since you can get invertible matrices that we terribly conditioned, we never care about the determinant or invertibility: we only ever compute $\rcond(A)$ to check for near-singularity, never the determinant.

_Def_: A problem is **well posed** if a solution exists and is unique, and **ill-posed** otherwise. An **instance** of a problem is a selection of parameters that describes the problem further.

_Def_: The **condition number of a problem** is the normalized reciprocal of the distance  to the nearest ill-posed instance.

For example, in the case of solving $Ax = b$, this will be $\frac{\|A\|}{d(A, M)}$, where $M = \{X \in \C^{n \times n} \mid \det(X) = 0 \}$ is the ill-posed manifold.

(Quiz problem: show that $d(A, M) = \|A^{-1}\|^{-1}$).

Luckily, if $A = U\Sigma V^T$, we can see that $\|\Delta A\| = \|\Delta \Sigma\|$ so the SVD is perfectly conditioned and so we will be able to find it to machine precision.

In some sense, every problem gives rise to a map $f: X \to Y$ from input data to an output solution; when the problem is well posed this is a function. When $X, Y$ have norms, then
$$
    \kappa_f(x) = \lim_{\delta \to 0} \sup_{\operatorname{RelErr}(x) < \delta}\frac{\operatorname{RelErr}(f(x))}{\operatorname{RelErr}(x)}
$$
is the (relative) condition number. Recall that the relative error is 
$$
    \operatorname{RelErr}(x) = \frac{\|\Delta x\|}{\|x\|}.
$$
Then, we immediately see that
$$
    \operatorname{RelErr}(f(x)) \leq \kappa_f(x) \operatorname{RelErr}(x) + o(\operatorname{RelErr}(x))
$$
or, with names,
$$
    \text{forward error} \lessapprox \text{condition number } \cdot \text{ backward error}
$$

## Stability

Stability, on the other hand will be a property of an algorithm. As above, each problem is some $f:X \to Y$ sending inputs to solutions; on the other hand, we only have some algorithm $\hat f:X \to Y$.

Suppose that we have some backward error $x + \Delta x$ and some forward error $y + \Delta y$, such that $\hat f(x + \Delta x) = f(x) = y + \Delta y$.

For example, if $A \in \operatorname{GL}(n)$, then the solution to $Ax = b$ is $f(A, b) = A^{-1}b$; for any algorithm $\hat f(A, b) = \hat x$, if we know $A$ exactly, then the forward error is
$$
    \|f(A, b) - \hat f(A, b)\| = \|A^{-1} b - \hat x\|
$$
and the backward error is 
$$
    \|A\hat{x} - b\|.
$$

If the problem is SVD, then we have $f(A) = (U, \Sigma, V)$; the forward errors are
$$
    \|U - \hat U\|, \|V - \hat V\|, \|\Sigma - \hat \Sigma\|
$$
but the backward error is just
$$
    \|A - \hat U \hat \Sigma \hat V^T\|.
$$

_Def_: We say an algorithm $\hat f$ is **backwards stable**, if for any $x \in X$, the computed $\hat y = \hat f(x)$ satisfies that
$$
    \hat y = f(x + \Delta x), \ \  \|\Delta x\| \leq \delta \|x\|
$$
for $\delta$ small.

_Def_: We say that $\hat f$ is **numerically stable** if for any $x \in X$, the computed $\hat y = \hat f(x)$ satisfies
$$
    \hat y + \Delta y = f(x + \Delta x), \ \ \|\Delta x\| \leq \delta \|x\|, \|\Delta y\| \leq \delta \|y\|
$$
for $\delta, \epsilon$ small.


## QR Decomposition

_Def_: Let $A \in \R^{m \times n}$ with $n \leq m$; then we can find a decomposition $A = QR$, where $Q \in O(m)$ is orthogonal and $R \in \R^{m \times n}$ is upper triangular. In particular, we have $R = \bmat{R_1 \\ 0}$ where $R_1$ is upper triangular in $\R^{n \times n}$. If $A$ is of full column rank, then $R_1$ is invertible. This is called the **full $QR$ decomposition**.

_Def_: In the same setup as above, partition $A = \bmat{Q_1 & Q_2}\bmat{R_1 \\ 0} = Q_1R_1$. This is the **condensed $QR$ decomposition**.

_Def_: In the same setup as above, if we take $\rank(A) = r$, then there is a permutation matrix $\Pi$ such that
$$
    A\Pi = Q \bmat{R_1 & S \\ 0 & 0}
$$
and so
$$
    A = Q \bmat{R_2^T & 0 \\ 0 & 0}Z^T \Pi^T
$$
where $R_1 \in \C^{r \times r}$ and we get $R_2$ by doing the full $QR$ decomposition on $\bmat{R_1^T \\ S^T} = Z \bmat{R_2 \\ 0}$. This is the **rank-retaining $QR$ decomposition**.

_Def_: For any $A \in \R^{m \times n}$, we may find $A = Q\bmat{L & 0 \\ 0 & 0}U^T$. This is the **complete orthogonal decomposition**.

### Solving Linear Equations

To solve $Ax = b$, we can take $A = QR$ so $Rx = Q^Tb$ and solve using back substitution. Alternatively, let $A\Pi = LU$; then we solve $Ly = b$ with back substitution, $Uz = y$ with forward substitution, and set $x = \Pi z$. 

Linearly constrained least squares problems are problems of the form 
$$
    \min\| Ax - b \|, \  \ \text{s.t. } C^Tx = d.
$$
where $A \in \R^{m \times n}, b \in \R^n, C \in \R^{n \times p}$, and $d \in \R^p$. 

- We could form the Lagrangian to get
$$
    \begin{cases}
        A^Ax - A^Tb + C\lambda = 0\\
        C^Tx - d = 0
    \end{cases}   
$$
which is a KKT constraint. Sometimes we care about $\lambda$, and this is fine (though $A^TA$ can often be ill-conditioned).
- Instead, we may use the QR decomposition: note that $A^TAx = A^Tb - C\lambda$ and thus $x = \hat x - (A^TA)^{-1}C\lambda$, where $\hat x = \argmin \|Ax - b\|_2$. Then, since $C^Tx = d$, we must have that
$$
    C^T(A^TA)^{-1}C\lambda = C^T\hat x - d.
$$
To avoid $A^TA$, set $A = Q\bmat{R \\ 0}$, set $W = R^{-T}C$ via backsubstitution, take the QR $W = Q_1R_1$, set $\eta = C^T \hat x - d$ and finally solve $R_1^TR_1 \lambda = \eta$ via backsubstitution.

- If we do not want to compute $\lambda$, if we take $p \leq n$, let $C = Q_2 \bmat{R_2 \\ 0}$, so that $\bmat{R_2^T & 0}Q_2^Tx = d$. Thus, if we take $Q_2^Tx = \bmat{u \\ v}$, we backsolve for $u$ such that $R_2^Tu = d$, we know that
$$
    \|b - Ax\|^2 = \|b - AQ_2Q_2^Tx\|^2 = \left\|b - \bmat{A_1 & A_2} \bmat{u \\ v}\right\|^2 = \|b - A_1u - A_2v\|_2^2.
$$

### Computing QR, LU, CO Decompositions

Either use Householder reflections or Givens rotations to compute QR. Compute the complete orthogonal decomposition via two QR decompositions. Compute the LU decomposition via Gauss/elimination matrices.

### Multiple RHS

Suppose that we have some sequence of equations $Ax_i = b_i$. In this case, we form $B = \mat{b_1 & b_2 & \dots & b_n}$ and solve $AX = B$. Everything from before carries through to $X = A^+B$, where $A^+$ is the pseudoinverse.

### Low Rank Updates

If you have a low-rank error, e.g. you solved $Ax = b$ but in fact needed $\hat A = A + uv^T$ instead, then you can solve the new system efficiently by using the Sherman-Morrison formula, or in general the Woodbury formula for higher rank corrections.
