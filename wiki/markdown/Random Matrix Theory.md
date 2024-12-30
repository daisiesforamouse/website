---
title: 'Topics in Random Matrix Theory'
subtitle: 'UChicago STAT 38520, Autumn 2024'
...

\newcommand{\N}{\mathbb{N}}
\newcommand{\R}{\mathbb{R}}
\newcommand{\Q}{\mathbb{Q}}
\newcommand{\Z}{\mathbb{Z}}
\newcommand{\C}{\mathbb{C}}

\DeclareMathOperator{\Jac}{Jac}
\DeclareMathOperator{\Ker}{Ker}
\DeclareMathOperator{\mesh}{mesh}
\DeclareMathOperator{\Var}{Var}
\DeclareMathOperator{\hm}{hm}
\DeclareMathOperator{\Tr}{Tr}

\let\temp\phi
\let\phi\varphi
\let\varphi\temp

\newcommand{\pa}[1]{\left(#1\right)}
\newcommand{\bra}[1]{\left[#1\right]}
\newcommand{\cbra}[1]{\left\{#1\right\}}

\newcommand{\mat}[1]{\begin{matrix}#1\end{matrix}}
\newcommand{\pmat}[1]{\pa{\mat{#1}}}
\newcommand{\bmat}[1]{\bra{\mat{#1}}}

\newcommand{\pfrac}[2]{\pa{\frac{#1}{#2}}}
\newcommand{\bfrac}[2]{\bra{\frac{#1}{#2}}}
\newcommand{\psfrac}[2]{\pa{\sfrac{#1}{#2}}}
\newcommand{\bsfrac}[2]{\bra{\sfrac{#1}{#2}}}

## Wigner's Semicircle Law

-------

Given a symmetric matrix $M$, we will always denote its eigenvalues in increasing order, i.e. as
$$
\lambda_1(M) \leq \lambda_2(M) \leq \cdots \leq \lambda_n(M)
$$
and the associated orthonormal eigenvectors
$$
M\psi_{\lambda_k(M)} = \lambda_k(M) \psi_{\lambda_k(M)}.
$$

_Def_: We say a matrix is **Rademacher** if is symmetric and has Rademacher entries on all off-diagonal entries and vanshes on the diagonal, e.g. of the form
$$
\bmat{
    0 & R_{12} & \cdots & R_{1n} \\
    R_{12} & 0 & \cdots & R_{2n} \\
    \vdots & \vdots & \ddots & \vdots \\
    R_{1n} & R_{2n} & \cdots & 0 \\
}
$$
where $R_{ij} = \pm 1$ with equal probability. If the diagonals are also Rademacher, we call it **fully** Rademacher.

The diagonals are zero to simplify computations - we will see later that due to universality, it does not really matter.

_Def_: The **empirical spectral distribution** (ESD) of a symmetric matrix $M$ is given as
$$
L_M(A) = \frac{1}{n} \sum_{k=1}^n 1\{\lambda_k(M) \in A\}.
$$
For over, in a single component $1 \leq i \leq n$, we say that the **ESD at $i$** is
$$
L_M^i(A) = \sum_{k=1}^n |\psi_{\lambda_k(M)}(i)|^2 1\{\lambda_k(M) \in A\}.
$$

**Prop**: The above are both probability measures.

_Proof_: Look at it.

_Def_: The **semicircle measure** is denoted as $\mu_{sc}$ and is given as
$$
\mu_{sc}(A) = \int_A \frac{\sqrt{4-x^2}}{2\pi} 1\{|x| \leq 2\}dx.
$$

**Theorem (Wigner's Semicircle Law)**: For any interval $I \subset \R$ and $\widetilde R_n = \frac{R_n}{\sqrt{n}}$ where $R_n$ are Rademacher,
$$
\lim_{n \to \infty} E[L_{\widetilde R_n}(I)] = \mu_{sc}(I)
$$
and for any sequence $i_1, i_2, \dots, i_n$ satisfying $1 \leq i_n \leq n$, 
$$
\lim_{n \to \infty} E[L_{\widetilde R_n}^{i_n}(I)] = \mu_{sc}(I).
$$

**Lemma (The Method of Moments)**: Let $\mu_n$ and $\mu$ be probability measures on $\mathbb R$ such that
$$
\int_{\R} x^p d\mu_n = \int_{\R} x^p d\mu,
$$
and let $\mu$ be determined by its moments. Then $\mu_n \to \mu$ weakly. Moreover, $\mu$ having bounded support is sufficient for it to be determined by its moments.

**Lemma**: For any symmetric $M \in \R^{n \times n}$ and integer $p \in \mathbb N$, we have that
$$
\int_{\R} x^p dL_M = \frac{1}{n} \Tr(M^p) = \frac{1}{n} \sum_{i_1, \dots, i_p = 1}^n M_{i_1, i_2}M_{i_2, i_3}\cdots M_{i_p, i_1}.
$$
Similarly,
$$
\int_{\R} x^p dL^i_M = \frac{1}{n} \Tr(M^p) = \frac{1}{n} \sum_{i_2, \dots, i_p = 1}^n M_{i, i_2}M_{i_2, i_3}\cdots M_{i_p, i}.
$$

_Proof (Semicircle Law)_: From the above lemmas, we simply need to show that
$$
\lim_{n \to \infty} n^{- \left(\frac{p}{2} + 1\right)} \sum_{i_1, \dots, i_p = 1}^n E \left[R_{i_1, i_2}R_{i_2, i_3}\cdots R_{i_p, i_1}\right] \to \int x^p d\mu_{sc}
$$
and similarly for the individual components.

In fact, it will be easier to define the complete graph of order $n$, denoted as $K_n$, as the graph of $n$ vertices and every possible edge. Then, let $W_p^{i,j}(K_n)$ be the set of all length $p$ walks starting from $i$ and ending at $j$. Then, the above reduces to
$$
\lim_{n \to \infty} n^{- \left(\frac{p}{2} + 1\right)} \sum_{i=1}^n \sum_{\nu \in W_p^{i,i}(K_n)}^n E \left[R_n(\nu_1, \nu_2)\cdots R_n(\nu_{p}, \nu_{p+1})\right] \to \int x^p d\mu_{sc}
$$

By independence, we see that the summand is $1$ if every edge is crossed an even amount of times, and $0$ otherwise, so this question reduces to counting the number of walks that cross each edge an odd amout of times. To count this, introduce $SW^{i,j}_p(K_n)$ as the set of standard walks, which visit nodes in increasing vertex order (e.g. you cannot visit vertex $1$, then $2$, and so forth); clearly each walk is equivalent to a standard walk by permuting the verticies, so we just need to compute
$$
\lim_{n \to \infty} n^{-\frac{p}{2}} \sum_{\nu \in SW_{p}(K_n)} P_{n-1, w(\nu) - 1} = \lim_{n \to \infty} n^{-\frac{p}{2}} \sum_{\ell = 2}^{p+1} P_{n-1, \ell - 1} \cdot |\{\nu \in SW_p(K_n), w(\nu) = \ell\}|
$$
where $w(\nu)$, the weight of $\nu$, is the total amount of vertices visited and $P_{a,b}$ is the number of ways of picking $b$ elements from a set of size $a$ with replacement.


Since we require walks to visit each edge twice, we can reduce this to
$$
\lim_{n \to \infty} n^{-\frac{p}{2}} \sum_{\ell = 2}^{\frac{p}{2}+1} P_{n-1, \ell - 1} \cdot |\{\nu \in SW_p(K_{\frac{p}{2} + 1}), w(\nu) = \ell\}| = \left|\left\{\nu \in SW_p(K_{\frac{p}{2}+1}), w(\nu) = \frac{p}{2} + 1 \right\}\right| = C_{\frac{p}{2}}
$$
which is the $\frac{p}{2}$-th Catalan number; one can check by straightfoward computation that this matches the moments of the semi-circle distribution.

We can say more as well. The convergence can be strengthened to almost sure convergence (which can be shown by looking at the variance of the above sums. Additionally, we have the following.

**Theorem (Bai-Silverstein)**: With the same notation as above,
$$
\sup_{x \in \R} |E[\mathcal L_{\widetilde R_n}((-\infty, x])] - \mu_{sc}((-\infty, x])| \lesssim n^{-\frac{1}{2}}.
$$

### Covariance Matrices

Consider a random vector $X \in \R^d$ with covariance matrix $\Sigma$. Then, given a random sample of i.i.d. copies $X_1, \dots, X_n$, we can naively estimate the covariance as
$$
\Sigma_n = \frac{1}{n} \sum_{k=1}^n X_k X_k^\top.
$$
We are interested in how well the spectrum of $\Sigma_n$ approximates the spectrum of $\Sigma$.

For fixed $d$, it is immediate from the law of large numbers that the spectra converge. However, this is not the case when $d$ is comprable to $n$.

_Def_: Let $\gamma \in (0, \infty)$, and let 
$$
a_\gamma = (1 - \sqrt{\gamma})^2 \ \text{ and } \ b_\gamma = (1 + \sqrt{\gamma})^2.
$$
The **Marchenko-Pastur measure** with parameter $\gamma$ is the measure which satisfies, for any $A \subset \R$,
$$
\mu^\gamma_{mp} (A) = 
\begin{cases}
    \int_A \frac{\sqrt{(x - a_\gamma)(b_\gamma - x)}}{2\pi \gamma x}1\{a_\gamma \leq x \leq b_\gamma\}dx & 0 < \gamma \leq 1 \\
    \left(1 - \frac{1}{\gamma} \right)1\{0 \in A\} + \int_A \frac{\sqrt{(x - a_\gamma)(b_\gamma - x)}}{2\pi \gamma x}1\{a_\gamma \leq x \leq b_\gamma\}dx & \gamma > 1
\end{cases}
$$


**Theorem**: Let the entries of $R \in \R^d$ be Rademacher; moreover let $R_1, \dots, R_n$ be i.i.d. copies of $R$ and let
$$
R_{d,n} = \frac{1}{n} \sum_{k=1}^n R_k R_k^\top
$$
be the covariance estimator. Then, if $\lim_{n \to \infty}\frac{d}{n} = \gamma$,
$$
\lim_{n \to \infty} E[L_{R_{d,n}}] \to \mu_{mp}^\gamma.
$$

_Proof_: The proof proceeds via the method of moments, but we now have walks on bipartite graphs of order $(d, n)$ instead.

As before, this convergence can be made almost sure.

### A Microscopic Result

Recall that the spectral norm $\|M\|_2$ of a symmetric matrix is just its top eigenvalue.

**Theorem**: For $\tilde R_n$ as in the semicircle law, we have
$$
\lim_{n \to \infty} \|\widetilde R_n\|_2 = 2
$$
almost surely.

_Proof_: The lower bound is clear: if $\liminf_{n \to \infty} \|\widetilde R_n\|_2 < 2$, then the strong semicircular law would be violated. The upper bound is much harder, and proceeds by the bound
$$
\|\widetilde R_n\|_2 = \max_{k} |\lambda_k(\widetilde R_n)| \leq \left(\sum_{k}\lambda_k(\widetilde R_n)^{2p}\right)^{1/2p} = \Tr(\widetilde R_n^{2p})^{1/2p}
$$
where you pick $p \asymp n^\alpha$ for some $0 < \alpha < 1/6$.

## Universality

------------------------

To prove a more general result, we have the following strategy.

1. Identify a canidate for a universal limit by checking a simple example. 
2. Identify a set of functions which characterizes convergence in distribution to that universal limit, which is relatively easy to control, and which places minimal assumptions on our sequence of random variables.
3. Combine the first two by showing that the behaviour of a general sequence is close to the simple sequence.

_Def_: A matrix is of **Wigner-type** if it is symmetric with independent entries.

**Theorem (Wigner's Universality)**: Let $X_n \in \R^{n \times n}$ be a sequence of Wigner-type matrices. Then, if for any $n \in \mathbb N, 1 \leq i,j, \leq n$, 

1. $X_{i,j}$ is symmetric with unit variance,
2. and for any $p \geq 3$, $\sup_{n} \sup_{i,j} E[|X^p_{i,j}|] < \infty$,

then $E[\mathcal L_{X_n / \sqrt{n}}] \to \mu_{sc}$ as $n \to \infty$.

The previous result on Rademacher matrices has provided a canidate for universality, under the usual heuristic that functions of independent random variables do not really depend on the distributions of those variables when the amount is large. To give an explicit bound on this heuristic (oft called the Lindeberg universality principle), we have the following bound.

**Theorem**: Let $X, Y \in \R^n$ be random vectors with independent entries, and let $f: \R^n \to \C$ be continuous. If there is some $p\in \mathbb N$ such that the following hold:

1. for every $1 \leq i \leq n$ and $1 \leq \ell \leq p - 1$,
$$
E[X_i^\ell] = E[Y_i^\ell];
$$
2. there is a constant $M < \infty$ such that
$$
\sup_{1 \leq i \leq n} (E[|X_i|^p] + E[|Y_i|^p]) \leq M;
$$ 
3. $f$ is $p$-times differentiable in every coordinate;

then we have the estimate
$$
  \left|E[f(X) - E[f(Y)]] \right| \leq \frac{M}{p!} \sum_{i=1}^n \left\|\frac{\partial^p f}{\partial x_i^p}\right\|_\infty.
$$

_Proof_: Induct on $n$, take a Taylor expansion, and it falls out.

### The Stieltjes Transform

_Def_: For a finite measure $\mu$ on $\R$, the Stieltjes transform of $\mu$ is
$$
S_\mu(z) = \int_\R \frac{1}{x-z}d \mu
$$
for $z \in \C \setminus \R$.

**Prop**: We have that:

1. for any $\mu$ and finite $z$, $|S_\mu(z)| < \infty$;
2. for any two $\mu, \nu$ which satisfy that $S_\mu(z) = S_\nu(z)$ for all $z \in \C \setminus \R$, we have $\mu = \nu$;
3. for any two reals $a < b$, 
$$
\lim_{\epsilon \to 0} \frac{1}{\pi} \int_a^b \text{Im}(S_\mu(t + i\epsilon))dt = \mu((a,b)) + \frac{1}{2}\left(\mu(\{a\}) + \mu(\{b\})\right).
$$

This third fact is called the **Stieltjes inversion formula**.

_Proof_: Propositions 1 and 3 are immediate from direct computation, and proposition 2 follows from 3 since 3 implies that $\nu, \eta$ assign identical mass to all intervals and atoms.

_Def_: We say that a measure $\nu$ on $\R$ is a **sub-probability** measure if $\nu(\R) \leq 1$; moreover, we say that a sequence of sub-probability measures $\nu_n$ **converges vaguely** to $\nu$ if one of the two equivalent conditions hold.

1.
$$
\lim_{n \to \infty} \nu_n((x, y]) = \nu((x, y])
$$
for every $x < y$ amd $\nu(\{x\}) = \nu(\{y\}) = 0$.
2. 
$$
\lim_{n \to \infty} \int_\R f(x) d\nu_n(x) = \int_\R f(x) d\nu(x)
$$
for every continuous function with $\lim_{x \to \pm \infty} f(x) = 0$.

**Theorem (Helly Selection)**: For every infinite sequence of probability measures $\mu_n$, there exists a sub-probability measure $\nu$ and a subsequence $n_i$ along which $\mu_{n_i}$ converges vaguely to $\nu$ as $i \to \infty$.

_Proof_: Look at the corresponding CDFs: you can find a subsequence along which they converge on a countable dense subset of $\R$.

**Theorem**: Let $\mu_n$ be a sequence of probability measures and $\mu$ be a probability measure; then $\mu_n \to \mu$ weakly if and only if $\lim_{n \to \infty} S_{\mu_n}(z) = S_\mu(z)$ for all $z \in \C \setminus \R$.

_Proof_ Forward direction is clear. In the other dirction, use Helly selection and uniqueness of the Stieltjes transform to show convergence of subsequences of subsequences.

**Prop**: Let $\mu$ be a compactly supported measure, with support contained in $[-R , R]$. Then $S_\mu$ satisfies the expansion
$$
S_\mu(z) = - \sum_{n=0}^\infty \frac{m_n}{z^{n+1}}
$$
where $m_n$ is the $n$-th moment of $\mu$.

_Proof_: Take a Taylor expansion after using the fact that
$$
\frac{1}{x - z} = -\frac{1}{z} \cdot \frac{1}{1 - \frac{x}{z}}.
$$

**Prop**: $S_{\mu_{sc}}(z) = \frac{-z + \sqrt{z^2 - 4}}{2}$.

_Proof_: Check that $S_{\mu_{sc}}(z)$ satisfies that $S_{\mu_{sc}}(z)^2 + zS_{\mu_{sc}}(z) + 1 = 0$.

### Proof of Universality

**Theorem (Wigner's Universality)**: Let $X_n \in \R^{n \times n}$ be a sequence of Wigner-type matrices. Then, if for any $n \in \mathbb N, 1 \leq i,j, \leq n$, 

1. $X_{i,j}$ is symmetric with unit variance,
2. and for any $p \geq 3$, $\sup_{n} \sup_{i,j} E[|X^p_{i,j}|] < \infty$,

then $E[\mathcal L_{X_n / \sqrt{n}}] \to \mu_{sc}$ as $n \to \infty$.

_Proof_:
We only need to show that
$$
\lim_{n \to \infty} |S_{E[\mathcal L_{X_n / \sqrt{n}}]}(z) - S_{E[\mathcal L_{R_n / \sqrt{n}}]}(z)| = 0.
$$

We can check that
$$
S_{E[\mathcal L_{X_n / \sqrt{n}}]} = E[f(X_n(i, j), 1 \leq i, j \leq n)]
$$
where
$$
f(X_n(i, j), 1 \leq i, j \leq n) = \frac{1}{n} \Tr \left[\left( \frac{X_n}{\sqrt{n}} - z I_n \right)^{-1} \right]
$$

Applying the Lindeberg universality principle, combined with careful derivative bounds, gives us what we want. 


## The Gaussian Orthogonal Ensemble

_Def_: Let $g_{ij}$ for $1 \leq i \leq j \leq n$ be independent centered Gaussians with variance $1 + \delta_{ij}$. Then, any $A_n = [g_{ij}]$ with $g_{ij} = g_{ji}$ is called a **Gaussian orthogonal ensemble**.

Sometimes we will write $A_n / \sqrt{n} = \tilde A_n$.

_Def_: For any matrix $A \in \R^{n \times n}$, the resolvent of $A$ is
$$
R_A(z) = (A - zI)^{-1}.
$$

**Prop**: The Stieltjes transform of $\mu_n = E[\mathcal L_{A_n / \sqrt{n}}]$ is given by
$$
S_{\mu_n}(z) = -\frac{1}{z}\left( 1 - \frac{1}{n} E[\Tr(\tilde A_n R_{\tilde A_n}(z)] \right)
$$

**Theorem**: Denote the eignevalues of a GOE $A_n$ by $\lambda_1, \dots, \lambda_n$. Then their joint eigenvalue density is given by
$$
\frac{1}{z_n} \cdot 1\{\lambda_1 \leq \cdots \leq \lambda_n\} \cdot \prod_{i < j} |\lambda_i - \lambda_j|e^{-\frac{1}{4}\sum_{i=1}^n \lambda_i^2}
$$
where $z_n = (2\pi)^{-\frac{n}{2}} \cdot 2^{-\frac{n(n+1)}{4}} \cdot \prod_{i=1}^n \frac{\Gamma(1/2)}{\Gamma(i / 2)}$ is a normalizing constant.
