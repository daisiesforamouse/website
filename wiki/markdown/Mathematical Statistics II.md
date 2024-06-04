---
title: 'Mathematical Statistics II'
subtitle: 'UChicago STAT 30200, Spring 2024'
...

\newcommand{\N}{\mathbb{N}}
\newcommand{\R}{\mathbb{R}}
\newcommand{\Q}{\mathbb{Q}}
\newcommand{\Z}{\mathbb{Z}}
\newcommand{\C}{\mathbb{C}}
\newcommand{\M}{\mathbb{M}}

\DeclareMathOperator{\Var}{Var}
\DeclareMathOperator{\Cov}{Cov}
\DeclareMathOperator{\RSS}{RSS}
\DeclareMathOperator{\Jac}{Jac}
\DeclareMathOperator{\Ker}{Ker}
\DeclareMathOperator{\Im}{Im}
\DeclareMathOperator*{\argmin}{argmin}
\DeclareMathOperator*{\argmax}{argmin}
\DeclareMathOperator{\proj}{proj}
\DeclareMathOperator{\rank}{rank}
\DeclareMathOperator{\subG}{subG}
\DeclareMathOperator{\subE}{subE}
\DeclareMathOperator{\Lap}{Lap}
\DeclareMathOperator{\supp}{supp}

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

# Statistics In High Dimensions 

-------

## Concentration Inequalities

In general, a **concentration inequality** is just something that looks like
$$ P[|X - E[X]| > t] \leq \phi(t) $$
for some function $\phi$ that hopefully decays quickly.

**Theorem (Markov)**:  If $X \geq 0$, then
$$ P(X \geq t) \leq \frac{E[X]}{t} $$

_Proof_: Condition on $X > t$.

**Theorem (Chebyshev)**: For square integrable $X$,
$$ P((X - E[X])^2 \geq t) \leq \frac{\Var(X)}{t^2}. $$

_Proof_: Apply Markov.

**Theorem (Chernoff)**: For centered $X$ which the following is defined,
$$
  P(X \geq \lambda) \leq E[e^{tX}]e^{-t\lambda}.
$$

_Proof_: Apply Markov.

### Sub-Gaussian Random Variables

_Def_: A random variable $X$ is said to be $\sigma^2$-**sub-Gaussian** if
$$
  E[e^{\lambda(X - E[X])}] \leq e^{\frac{\lambda^2 \sigma^2}{2}}
$$
for all $\lambda$. We put $X \sim \subG(\sigma^2)$.

_Example_: Gaussians, clearly; also any bounded random variable.

**Prop:** For $X \sim \subG(\sigma^2)$,
$$
P(|X - \mu| \geq t) \leq 2e^{-\frac{t^2}{2\sigma^2}}.
$$

**Lemma**: For $X_1, \dots, X_n$ independent and $\sigma_i^2$-sub-Gaussian, $\sum_{i=1}^n X_i$ is $\sum_{i=1}^n \sigma_i^2$-sub-Gaussian.

**Theorem**: Given any centered random variable $X$, the following are equivalent.

1. $X$ is $\sigma^2$-sub-Gaussian.
2. There is a constant $c$ and a centered Gaussian $Z$ such that 
$$
P(|X| > s) \leq cP(|Z| > s).
$$
3. There is a constant $\theta$ such that
$$
  E[X^{2k}] \leq \frac{(2k)!\theta^{2k}}{2^k k!}.
$$
4. There is a constant $\sigma$ such that
$$
  E[e^{-\frac{sX^2}{2\sigma^2}}] \leq \frac{1}{\sqrt{1 - s}}, \ \ \forall s \in (0, 1).
$$

_Proof_:

- $(1) \implies (3)$. We will have for $Z \sim N(0, 2\sigma^2)$ that
$$
\frac{P(X > t)}{P(Z > t)} \leq 4\sqrt{\pi} e.
$$
Now apply Chernoff to $X$, and use Mill's ratio for $Z$, which yields 
$$
P(Z > t) \geq \left( \frac{\sqrt{2} \sigma}{t} - \left( \frac{\sqrt{2}\sigma}{t} \right)^3 \right) \frac{e^{\frac{-t^2}{4\sigma^2}}}{\sqrt{2\pi}}.
$$
Now simplify with $t = \sqrt{2}\sigma$ so that
$$
\frac{P(X > t)}{P(Z > t)} \leq 4\sqrt{\pi}e
$$
holds for all $t \in [0, \sqrt{2}\sigma]$. Similar for $t > 2\sigma$.

- $(2) \implies (3)$. Here we simply compute
$$
E[X^{2k}] = \int_0^\infty P(X^{2k} > s)ds  \leq \int_0^\infty cP(|Z| > s^{1/2k}) ds = cE[|Z^{2k}|] = c\frac{(2k)!\tau^{2k}}{2^k k!}
$$
as desired.

- $(3) \implies (1)$. Taylor expand the MGF if $X$ is symmetric. Otherwise, you can apply Cauchy-Schwarz to get very crude bounds on the odd moments from the even moments.


- $(1) \implies (4)$. Just compute the integrals.

- $(4) \implies (1)$. We have a bound
$$
  e^u \leq u + e^{gu^2/16}
$$
that we apply to $e^{\lambda X}$. Apply this and do a bunch of miscellanous bounds to get what you want (briefly, split into cases depending on $\lambda$, and use $\frac{1}{\sqrt{1 - s}} \leq e^s$ for one half; for the other half, simply use $2ab \leq c a^2 + c^{-1}b^2$).

**Theorem (Hoeffding)**: Let $X_i$ be independent $\sigma^2_i$-sub-Gaussian random variables. Then 
$$
\bar X = \frac{1}{n}\sum_{i=1}^n (X_i - E[X_i])
$$ is sub-Gaussian with parameter $\frac{1}{n^2} \sum_{i=1}^n \sigma_i^2$.

_Def_: A $X \in \R^d$ is sub-Gaussian with parameter $\sigma^2$ if
$$
u^\top (X - E[X])
$$
is $\sigma^2$ sub-Gaussian for all $u \in S^{d-1}$. Similarly, for $X \in \R^{d \times T}$, $X$ is sub-Gaussian with parameter $\sigma^2$ if $u^\top (X - E[X])\sigma$ is $\sigma^2$ sub-Gaussain. In these cases, we use the notation $X \sim \subG_{d \times T}(\sigma^2)$.

### Sub-Exponential Random Variables

Let $X \sim \Lap(1)$. Then $P(|X| > t) \leq e^{-t}$, and
$$
E[e^{\lambda X}] = \frac{1}{1 - \lambda^2}
$$
is well-defined for $\lambda \in [0, 1)$; in fact for $|\lambda| < 1/2$,
$$
E[e^{\lambda X}] \leq e^{2\lambda^2}.
$$

In fact, something more general is true. 

**Theorem**: Suppose that for some random variable $X$ and positive constant $b$,
$$
P(|X| > t) \leq 2e^{-2t/b}.
$$
Then,

1. $E[|X|^k] \leq b^k k!$;
2. $E[|X|^k]^{1/k} \leq 2bk$;
3. $E[e^{\lambda X}] \leq e^{2\lambda^2b^2}$ for $|\lambda| \leq 1/2b$.

_Proof_:

- $(1)$. Use the fact that
$$
E[|X|^k] = \int P(|X|^k > t)dt
$$
- $(2)$. $k! \leq k^k$.
- $(3)$. Taylor expand the exponential.

_Def_: $X$ is said to be **sub-Exponential** with parameters $(\tau^2, b)$ if
$$
    E[e^{\lambda (X - E[X])}] \leq e^{\tau^2\lambda^2 / 2} 
$$
for all $|\lambda| \leq 1/b$. We also say $X \sim \subE(\tau^2, b)$.

**Prop**: If $X \sim \subE(\tau^2, b)$, then
$$
P(X - E[X] > t) \leq \exp\left(-\frac{1}{2} \min\left(\frac{t^2}{2b}, \frac{t}{b}\right) \right) \leq \begin{cases}e^{-t^2/2\tau^2} & t \in [0, \tau^2/b] \\ e^{-t/2\tau^2} & t > \tau^2/b\end{cases}.
$$

_Proof_: Just apply Chernoff.

**Lemma (Bernstein's Condition)**: Let $b > 0$; if
$$
  E[(X - E[X])^k] \leq \frac{1}{2} k! b^{k-2} \sigma^2,
$$
then $X$ is sub-Exponential.

_Proof_: Taylor expand the MGF.

In fact, if $X$ satisfies Bernstein's condition, we have
$$
E[e^{\lambda (X - E[X])}] \leq \exp \left(\frac{\lambda^2 \sigma^2}{2(1 - b|\lambda|)} \right)
$$
for all $|\lambda| < 1/b$ and
$$
P(|X - E[X]| > t)\leq 2\exp \left( -\frac{t^2}{2(\sigma^2 + bt)} \right)
$$
for all $t > 0$.

**Theorem:** If $X$ is a centered random variable, the following are equivalent.

1. $X \sim \subE(u^2, \alpha)$.
2. There is some $c_0 > 0$ such that $E[e^{\lambda X}] < \infty$ for all $|\lambda| \leq c_0$.
3. There are $c_1, c_2 > 0$ such that
$$
P(|X| \geq t) \leq c_1e^{-c_2t}, \ \ \forall t > 0.
$$
4. 
$$
\gamma = \sup_{k} \left[ \frac{E[|X|^k]}{k!} \right]^{1/k} < \infty.
$$

**Lemma (One-sided Berstein's Inequality)**: If $X \leq b$ almost surely,
$$
E[e^{\lambda(X - E[X])}] \leq \exp \left( \frac{\lambda^2E[X^2]}{2(1 - b\lambda / 3)} \right).
$$

**Corollary**: 
$$
P \left(\frac{1}{n} \sum_{i=1}^n (X_i -E[X_i]) > \delta \right) \leq \exp \left(\frac{-n\delta}{2(n^{-1}\sum E[X_i^2] + b\delta/3)} \right).
$$

### Maximal Inequalities

We first start by think about maximums over finite sets.

**Theorem**: Let $X_1, \dots, X_n$ be a set of $n$ random variables such that $X_i \sim \subG(\sigma^2)$. Then
$$
  E[\max X_i] \leq \sigma \sqrt{2 \log(n)}, \ \ E[\max |X_i|] \leq \sigma \sqrt{2 \log(2n)}.
$$

_Proof_: For $s > 0$, consider $s^{-1}E[\max sX_i]$ and use a Chernoff bound alongside a judicious choice of $s$.

Another realtively simple case is maxima over convex polytopes.

_Def_: A convex polytope is a compact set with a finite set of vertices $V$ such that 
$$
P = \left\{ X \mid X = \sum_{i=1}^{|V|} \lambda_i v_i \right\}.
$$

**Lemma**: The minima/maxima of a linear form over a polytope is achieved at a vertex.

_Proof_: Remember 8th grade.

**Theorem**: Let $P$ be a polytope with $n$ vertices $v_i$, and $X \in \R^d$ a random variable such that $v_i^\top X \sim \subG(\sigma^2)$. Then the earlier bound holds again.

_Proof_: Apply the last lemma.

We can now talk about more general bodies.

_Def_: Fix a set $T \subset \R^d$, and $\epsilon > 0$. A set $N(\epsilon, T, d)$ is said to be an $\epsilon$-**net** with respect to a distance $d$ if there are $\theta_i \in N(\epsilon)$ such that $d(\theta_i, z) < \epsilon$ or any $z \in T$. The $\epsilon$-**covering number** is the minimal cardinality of such a set.

_Def_: A $\delta$**-packing** of a set $T$ with a distance $d$ is a set $P(\delta, T, d)$ is a set such that
$$
  d(\theta_i, \theta_j) > \delta \ \ \forall \theta_i, \theta_j \in P(\delta, T, d).
$$
The $\delta$-packing number is the maximal cardinality of such a set.

**Lemma**: The $L^2$ unit ball $B_2$ has a $\epsilon$ net of cardinality
$$
N \leq \left(\frac{3}{\epsilon} \right)^d.
$$

_Proof_: Compute a loose bound on the packing number and note that a packing is also a covering.

**Theorem**: Let $X \in \R^d$ be $\subG(\sigma^2)$. Then,
$$
E[\max_{\theta \in B_2} \theta^\top X] = E[\max_{\theta \in B_2} |\theta^T X|] \leq 4\sigma \sqrt{d}
$$
and with probability at least $1 - \delta$,
$$
\max_{\theta \in B_2} \theta^\top X  \leq 4\sigma\sqrt{d} + 2\sigma \sqrt{2 \log(1/\delta)}.
$$

_Proof_ Let $N$ be a $\epsilon = 1/2$ net for $B_2$. Then $|N| \leq 6^d$ and for $\theta \in B_2$, $\theta = z + x$ for some $z \in N$, $x$ of size $1/2$. Then bound
$$
\max \theta^\top X \leq \max z^\top X + \max x^\top X = \max z^\top X + \frac{1}{2} \max \theta^\top X.
$$
Now use the finite set bound.


## Linear Regression

The setup is as follows: set $Y = X \beta + \epsilon$ for some **fixed** design matrix $X$ and $\epsilon \sim \subG(\sigma^2)$ with independent components. We abbreviate this model by $M$.

### Unconstrained OLS

In this case, we define
$$
    \hat \theta^{LS} \in \argmin_{\theta} \|Y - X \theta\|^2
$$
and
$$
    \hat \mu = X\hat \theta^{LS}.
$$
**Prop**: It must be that
$$
X^\top \hat \mu^{LS} = X^\top y
$$
and so we may take
$$
\hat \theta^{LS} = (X^\top X)^\dagger X^\top Y.
$$

_Proof_: I mean, just look.

**Theorem**: Assume that $M$ holds. The least squares estimate $\hat \theta^{LS}$ satisfies
$$
E[MSE(\hat \theta^{LS})] = E \left[ \frac{1}{n} \|X (\theta^* - \hat \theta^{LS}) \|^2 \right] \leq \sigma^2 \cdot \frac{r}{n}
$$
where $r$ is the rank of $X^\top X$. Similarly,
$$
MSE(\hat \theta^{LS}) \leq \frac{\sigma^2(r + \log(1 / \delta))}{n}
$$
with probability at least $1 - \delta$.

_Proof_: We have
$$
\|Y - X\hat \theta^{LS}\|^2 \leq \|Y - X\theta^* \|^2 \leq \|\epsilon\|^2 
$$
so
$$
\|X(\theta^* - \hat \theta^{LS}) \|^2 \leq 2\epsilon^\top X(\hat \theta^{LS} - \theta^*).
$$
Divide by $\|X(\theta^* - \hat \theta^{LS})\|$, so that we have the bound 
$$
\|X(\theta^* - \hat \theta^{LS}) \|^2 \leq 4\sup_{u \in B_2} (\epsilon^\top u)^2.
$$
But we need to be a little more careful: reduce to the column space of $X$ (by mapping $\epsilon$ into the column space) and apply the bounds in the previous sections.

### Sparsity Constraints

_Def_: Define the **support** of a vector $\theta$ as the number of nonzero entries.

In this case, let $\theta^* \in K$ and 
$$
\hat \theta^{LS}_K = \argmin_{\theta \in K} \|Y - X\theta\|^2.
$$
As before, we have
$$
\|X(\hat \theta - \theta^*)\|_2^2 \leq 2\epsilon^\top X(\hat \theta - \theta^*).
$$

**Theorem**: Let $K = B_1$; if $\theta^* \in K$ and $\max_{j \leq d}\|X_j\| \leq \sqrt{n}$, then
$$
E[MSE(\hat \theta^{LS}_K)] \leq \sigma \sqrt{\frac{\log d}{n}}
$$

### The Gaussian Sequence Model

In the model $Y = X\theta + \epsilon$, we take an orthogonality condition on $X$ for simplicity and thus reduce to $Z = \theta + \epsilon$ for $\epsilon \sim N(0, \sigma^2/n)$.

In the most basic case, we can just use a threshholding estimator, e.g. 
$$
\hat \theta_j^2 = \begin{cases}
    Z_j & |Z_j| > \tau \\
    0 & \text{otherwise}
\end{cases}.
$$

**Theorem**: Let $\tau^2 = 4\sigma^2 \log(n / k)$. Then
$$
E[\|\hat \theta^2 - \theta^*\|^2] \leq C\frac{k\sigma^2}{n}(1 + \log(n / k)).
$$


We can prove similar results for soft thresholding estimators as well.

To drop the orthogonality constraint on $X$ is more difficult. Consider first the noiseless setting; again $Y = X\theta$, where $X \in \R^{n \times d}$ for $d \gg n$.
In this case, we can consider this as an optimization problem of minimizing $\|\theta\|_1$ such that $Y = X\theta$.

_Def_: For a set $S \subset \{1, \dots, d\}$, **the critical cone**  is the set
$$
\C(S) = \{ \Delta \in \R^d \mid \|\Delta_{S^c}\| \leq \|\Delta_S\| \}
$$

_Def_: A matrix $X$ is said to satisfy the **restricted null space** property with respect to $S$, if
$$
\ker(X) \cap \C(S) = \{ 0 \}.
$$

**Theorem**: The following are equivalent.

- $X$ satisfies the restricted null space property with respect to $S$.
- For any $\theta^*$ with nonzero entries on only $S$, then $\theta^*$ is the unique solution to the earlier optimization problem.

_Proof_: Draw a picture.

We can give a few examples of matrices which obey the restricted null property (RNP). 

_Def_: The **pairwise incoherence** is 
$$
\delta_{PW}(X) = \max_{j, k} \left| \frac{\left\langle X_j, X_j \right\rangle}{n} - \delta_{jk}\right|.
$$

**Theorem**: If $\delta_{PW}(X) < \frac{1}{3s}$, the RNP holds for all subsets of cardinality less than $s$.

_Def_: $X$ satisfiesthe **restricted isometry property** with parameter $\delta_S(X)$ if
$$
    \left\| \frac{(X^\top X)_{S, S}}{n} - I_S \right\|_2 \leq \delta_S(X).
$$

**Theorem**:
$$
\delta_{PW}(X) \leq \delta_s(X) \leq s \delta_{PW}(X).
$$

**Theorem**: If the RIP constant of order $2s$ is bounded as $\delta_{2s}(X) \leq \frac{1}{3}$ the uniform RNP holds for any $|S| \leq s$.

_Proof_: Take a vector in $\mathbb C(S) \cap \ker(X)$ and split it into blocks of length $s$.

### Relaxed Basis Pursuit

Consider the problem of 
$$
\argmin \|X\theta - Y\|^2 \text{ subject to } \|\theta\|_1 \leq b
$$
where $Y = X\theta + \epsilon$.

If $\|\theta^*\|_1 = b$, then $\Delta = \theta^* - \hat \theta \in \C(S)$. In fact, if the noise is $\subG(\sigma^2)$, we can get that
$$
\|\Delta\|^2 \leq \frac{32}{\kappa_{min}}\frac{k\sigma^2}{n} \left(\log(2d) + t^2\right)
$$
with probability at least $1 - e^{-t^2}$. Here
$$
\frac{\|X\Delta\|^2_2}{n} \geq \kappa_{min} \|\Delta\|^2_2.
$$

Of course, this is equivalent to the LASSO. The proof looks roughly the same: you 1) localize the error and then 2) show that the error is controllable on that localization with high probability.

**Theorem**: Suppose $\hat \theta$ is the solution of the LASSO problem with $\lambda_n \geq 2\frac{\|X^\top \epsilon\|_\infty}{m}$
$$
\argmin \frac{1}{n}\|Y - X\theta\| + \lambda\|\theta\|_1.
$$
Then,

1. any optimal solution $\hat \theta$ is such that 
$$
\frac{\|X(\hat \theta - \theta^*)\|}{n} \leq 12\|\theta^*\|_1 \lambda_n
$$
2. and if $\theta^*$ is supported on a subset $S$ of cardinality $s$, and $X$ satisfies the $\mu$-restricted convexity property on $\mathbb C_3(S)$, then
$$
\frac{\|X(\hat \theta - \theta^*)\|_2^2}{n}  \leq \frac{9}{\mu} s\lambda_n^2.
$$

If you have random designs, then you need to use a theorem that implies the usage of the strong convexity codition; for example, with Gaussian matrices, this works out with a suitable high probability result.

## Regularized M-Estimators

$M$ is for minimization. 

Consider a sequence of random variables $Z_i \sim P_\theta$, as well as

1. some loss function $L$, so that we may take
$$
   \hat \theta = \argmin_{\theta \in \Theta} E[L_\theta(Z)]
$$
2.  and some regularizer $\Phi$, so that we may enforce a certain type of structure
$$ \hat \theta = \argmin_{\theta \in \Theta} E[L_\theta(Z)] + \lambda \Phi(\theta). $$

The regularizer admits many forms depending on the type of problem you are considering. If you have blocking stucture in $\theta$, for example, you can take norms of blocks in $\theta$, or if you want some sort of smoothness represented by a graph $G$, you can consider $\theta^\top L_G \theta$ where $L_G$ is the graph Laplacian.

We need a things to get well-controlled behavior. 

_Def_: Take some subset $\mathbb M$; we need that for $( \alpha, \beta ) \in ( \mathbb M, \overline{\overline{M}^\perp} )$,
$$ \Phi(\alpha + \beta) = \Phi(\alpha) + \Phi(\beta). $$
We then say that $\Phi$ is **decomposable**.

_Def_: We consider the **dual norm**, which is given as
$$ \Phi^*(v) = \sup_{\Phi(u) \leq 1} \left\langle u, v \right\rangle. $$

We now need a few different key pieces.

First, the corresponding "good event" that we need (analagous to $\lambda > \|X^\top \epsilon\|_\infty$) is
$$ \mathbb G(\lambda_n) = \left\{ \Phi^*(\nabla L_n(\theta)) \leq \frac{\lambda_n}{2} \right\}. $$

**Prop**: If $L_n$ is a convex function and $\Phi$ is a decomposable norm over $\mathbb M, \overline{\M}^\perp$, then on $\mathbb G(\lambda_n)$, we have
$$ \Delta \in C_{\theta^*}(\mathbb M, \overline{\M}^\perp) = \left\{ \Phi(\delta_{\overline{\M}^\perp}) \leq 3\Phi(\Delta_{\mathbb M}) + 4\Phi(\theta^*_{\overline{\M}^\perp}) \right\}.
$$

_Proof_: The basic inequality becomes
$$ L_n(\theta^* + \Delta) + \lambda \Phi(\theta^* + \Delta) \leq L(\theta^*) + \lambda \Phi(\theta^*). $$
Then, by decomposition,
$$ \Phi(\theta^* + \Delta) - \Phi(\theta^*) \geq \Phi(\Delta_{\overline{\M}^\perp}) - \Phi(\Delta_{\mathbb M}) - 2\Phi(\theta^*_{\overline{\M}^\perp}) $$
and on the good event, by convexity,
$$ L_n(\theta^* + \Delta) - L_n(\theta^*) \geq - \frac{\lambda_n}{2} \left( \Phi(\Delta_{\overline{}^\perp}) + \Phi(\Delta_{\M}) \right) $$
and we can get
$$ \lambda_n \Phi(\Delta_{\overline{}^\perp}) \leq 3\Phi(\Delta_{\M}) + 4\Phi(\theta^*_{\overline{}^\perp}). $$

We also need a well-behaved cost function
$$ \epsilon_n(\Delta) = L_n(\theta^* + \Delta) - L(\theta^*) - \left\langle \nabla L_n(\theta^*),\Delta \right\rangle $$
which satisfies a strong convexity condition.

For a given norm and regularizer denoted as above, the cost function must satisfy with some positive radius
$$ \epsilon(\Delta) \geq \frac{\kappa}{2} \|\Delta\|^2 - \tau^2_n \Phi^2(\Delta)$$

We also need a subspace Lipschitz constant.

_Def_: For any subspace $S$, the subspace Lipchitz constant is
$$ \Psi(S) = \sup_{u \neq 0}\frac{\Phi(u)}{\|u\|} $$

Then, we can finally prove a result that we like.

**Theorem**: Conditional on the good event,
1. any optimal solution satisfies the bound
$$ \Phi(\theta^* - \hat \theta) \leq C \Phi(\M)(\|\theta^* - \hat \theta\| +  \Phi(\theta^*_{\overline{\M}^\perp}))$$ and

2. if $\tau^2_n \leq \frac{\kappa}{C_0}$ and $\epsilon_n(\M, \overline{\M}^\perp) \leq R$ then
$$ \|\theta^* - \hat \theta\|_2^2 \leq C_1 \frac{\lambda_n^2}{\kappa^2} \Phi^2(\overline{\M}) + \frac{C_2}{\kappa}(\lambda_n \Phi(\theta^*_{\M^\perp}) + 16 \tau_n^2 \Phi^2(\theta^*_{\M}))$$

## Matrix Estimation

The key problem here is generally something of the form
$$ \|Y - \theta\|_F + \lambda \|\theta\|_{n}. $$

First, a few key lemmas.

**Lemma (Weyl)**: Let $A, B$ with singular values $\sigma_i(A), \sigma_i(B)$. Then,
$$ \max |\sigma_k(A) - \sigma_k(B)| \leq \|A - B\|_{op} $$.

**Lemma (Hoffman-Wielandt)**:
$$ \sum_{k=1}^n |\sigma_k(A) - \sigma_k(B)|^2 \leq \|A - B\|_F^2 $$

 **Lemma (Holder's Inequality)**:
$$ \left\langle A, B \right\rangle \leq \|A\|_p \|B\|_q$$
where the relevant norms are Schatten norms.

**Theorem (Eckart–Young–Mirsky)**: Let $A = \sum_{i=1}^R \lambda_i u_i v^\top_i$. Then, the best rank $r$ approximation is $\sum_{i=1}^r \lambda_i u_i v^\top_i$.

We start with the following.

**Theorem**: Let $Y \in \R^d$ be  a random vector with $E[YY^\top] = I$. Then, if $X = \Sigma_X^{1/2} Y$ and $\hat \Sigma = \frac{1}{n}X^\top X$m then
$$ \|\hat \Sigma  - \Sigma \|_{op} \leq \|\Sigma \|  \left( \sqrt{\frac{d + \log(1/\delta)}{n}}  \lor \frac{d + \log(1 / \delta)}{n}\right) $$
with probability at least $1 - \delta$.

## Covariance Estimation
 
The above bound is not particularly interesting if $d > n$, since it involves $\frac{d}{n}$. We can do better in the case of sparsity.

Suppose that $\Sigma$ is sparse; then we may consider something like
$$ \hat \Sigma = T_\lambda(S) = T_\lambda \left( \frac{X^\top X}{n} \right) $$
wherer $T_\lambda$ is some thresholding operator, e.g. $T_\lambda(u) = 1_{|u| \geq \lambda}$.

Then $\|T_\lambda(S)\|_{op} \leq s$, where $s$ is the maximal number of nonzero entries in a row.


**Theorem**: If $X_i$ is a sequence of $\subG(\sigma^2)$ i.i.d. random variables with covariance matrix $\Sigma$, then if $n > \log(d)$, for all $\delta > 0$, if 
$$ \lambda_n = 8\sigma^2 \sqrt{\frac{\log d}{n}}  + \sigma^2 \delta$$
then the probability of $\|T_\lambda(\hat \Sigma) - \Sigma\|_{op} \geq 2 \|\Sigma\|_{op}\lambda_n$ is small. 

Consider now a simple model. A spiked covariance matrix is 
$$ \Sigma = \theta v v^\top + I$$
which represents the data generating mechanism
$$ X_i = \sqrt{\theta}uv + \epsilon $$
where $u, \epsilon$ are normal.

_Def_: The cosine distance between two unit vector $u, v$ is $\angle(u, v) = \arccos(|u^\top v|)$.

**Theorem (Davis-Kahan)**: Let $A, B$ be 2 PSD matrices, and $(\lambda_i, u_i)$ be eigenvalue/vector pairs of $A$ and $(\mu_i, v_i)$ eigenvalue/vector pairs for $B$. Then,
$$\sin(\angle(u_i, v_i)) \leq \frac{2}{\max(\lambda_1 - \lambda_2, \mu_1 - \mu_2)} \|A - B\|_{op}$$
and moreover
$$\min_{\epsilon \in \{\pm 1\}} |\epsilon u_i - v_i|_2^2 \leq 2 \sin^2(\angle(u ,v)).$$

In the case of the spiked covariance matrix, the eigengap is exactly $\theta$. So with high probability
$$\sin(\angle(\hat v_1, v_1)) \leq \frac{2\|S - \hat \Sigma\|_{op}}{\theta} \leq \frac{2(\theta+1)}{\theta}\sqrt{\frac{d+\log(1/\delta)}{n}}.$$

**Theorem**: Under the spiked covariance model, if $\|v_1\|_0 \leq s$, then the $k$-sparse largest eigenvector of $\hat \Sigma$ satisfies
$$ \min_{\epsilon \in \{-1, 1\}}\|\epsilon \hat v_1 - v_1\|_2 \leq \frac{\theta + 1}{\theta} \sqrt{\frac{s\log(ed/s) + \log(1/\delta)}{n}} \lor 1.$$

Of course, this is not necessarily tractable (it's not a convex optimization problem); instead we can bound the 1-norm $\|v\|_1 \leq \lambda$ instead.

Alternatively, when we extend to more spikes, (and in some sense more principal components), we get some issues as eigenvectors are not uniquely determined: in that case you have to consider different losses, say.

# Classification

-------

Consider $(X, y)$ for $y \in \{0, 1\}$. The goal is to find an estimator $h: \mathcal X \to \{0, 1\}$ which minimizes the risk
$$ R(h) = P[y \neq h(X)]. $$

Note that the optimal one is clearly just the Bayes classifier $h^*$; however, we often don't have access to the joint or conditional law of $(X, y)$.

**Theorem**: For any classifier $h$,
$$R(h) - R(h^*) = E[|2\eta(x) - 1| \cdot 1_{h(X) \neq h^*(X)}]$$
where
$$\eta(X) = E[y \mid X].$$

_Def_: The excess risk is
$$\epsilon = E[R(\hat h) - R^*]$$
and the empirical risk is
$$R_n(h) = \frac{1}{n} \sum_{i=1}^n 1_{y_i \neq h(X_i)}.$$

_Def_: Let $\mathcal H$ be a set of classifiesr. The empirical risk of the classifier is
$$R_n(\hat h_n) = \argmin_{h \in \mathcal H} R_n(h).$$
Also, under the split
$$R(\hat h_n)  - R^* = [ R(\hat h_n) - \min_{h \in \mathcal H}R(h) ] + [ \min_{h \in \mathcal H}R(h) - R^* ]$$
the first term is the stochastic error and the latter the approximation error.

The goal is then to show that
$$E[R(\hat h_n)] \leq \min_{h \in \mathcal H} R(h) + \Delta_n(\mathcal H)$$
for some $\Delta_n \to 0$.

Alternatively, sometimes we wish for thing with high probability, i.e.
$$R(\hat h_n) \leq \min_{h \in \mathcal H}R(h) + \Delta_n(\mathcal H, \delta).$$

**Lemma**: Let $\mathcal H$ be a set of classifiers; the stochastic error of $\hat h_n$ satisfies that
$$R(\hat h_n) - R(h_{\mathcal H}) \leq 2\sup_{h \in \mathcal H} |R_n(h) - R(h)|$$

_Proof_: Consider $\epsilon > 0$ and $h_\epsilon$ such that $R(h_\epsilon) \leq \inf_{h \in \mathcal H} R(h) + \epsilon$; then we get that 
$$R(\hat h_n) - R(h_{\mathcal H}) =R(\hat h_n) - R_n(\hat h_n) + R_n(\hat h_n) - R(h_{\mathcal H} \leq R(\hat h_n) - R_n(\hat h_n) + R_n(h_\epsilon) - R(h_\epsilon) + \epsilon.$$
The result follows.

**Theorem**: Suppose $|\mathcal H| < \infty$. Then, with probability $1 - \delta$,
$$R(\hat h_n) \leq \min_{j} R(h_j) + \sqrt{\frac{2}{n} \log(2m\delta^{-1})}.$$

To do more than finite families of classifiers, it will be simpler to consider a classifier as simply a class of sets, e.g. $h(x) = 1_{x \in A}$; to this extent, write $\bar \mu(A) = P(y \neq 1_{X \in A})$ and $\bar \mu_n(A) = \frac{1}{n} \sum_{i=1}^n 1_{Y_i \neq 1_{X \in A}}$ and $\mu(A) = P(X \in A)$ and $\mu_n(A) = \frac{1}{n} \sum 1_{X_i \in A}$. 

Now, we can more easily control
$$\xi = \sup_{A \in \mathcal A}|\mu(A) - \mu_n(A)|$$
via a symmetrization argument. That is, if $X_1', \dots, X_n'$ are independent of observations $X_1, \dots, X_n$ with the same distribution,
$$E[\xi] = E[\sup |\mu_n(A) - \mu(A)|] = E[\sup |E[\mu_n(A) - \mu'_n(A) \mid X_1, \dots, X_n]|] \leq E[\sup |\mu_n(A) - \mu'_n(A)|].$$

_Def_: Let $B$ be a finite set of vectors; the Rademacher complexity of $B$ is given by
$$R_n(B) = \frac{1}{n} E_\sigma \left[ \sup_{b \in B} \sum_{i=1}^n \sigma_i b_i\right]$$
where $\sigma_i$ are i.i.d. Rademacher.

_Def_: Let $\mathcal A$ be a class of measurable sets $A$, the binary fingerprint of $X_1^n$ on $\mathcal A$ is the set of vectors defined as
$$\mathcal A(X_1^n) = \{b = (b_1, \dots, b_n) \in \{0, 1\}^n, b_i = 1_{X_i \in A}\}.$$
Then clearly 
$$E[\sup |\mu_n(A) - \mu(A)|] \leq 2E[R_n(\mathcal A(X_1^n))]$$

**Lemma**: Let $B = (b^{(1)}, \dots, b^{(n)})$ be as above; then
$$R_n(B) \leq \max_{j=1, \dots, n} \|b^{(j)}\|_2 \frac{2\log(n)}{n}.$$
Moreover, the maximum over the simplex spanned by $B$ is the same as the maximum over $B$.

_Def_: We call the **shattering coefficient** the function 
$$S_{\mathcal A}(n) = \max_{X_1^n} |\mathcal A(X_1^n)|.$$

**Theorem**: For any class of sets $\mathcal A$,
$$E[\sup_{A \in \mathcal A} |\mu_n(A) - \mu(A)|] \leq 2 \sqrt{\frac{2 \log(2S_{\mathcal A}(n))}{n}}.$$

_Def_: The **VC dimension** of $\mathcal A$ is
$$VC(\mathcal A) = \sup \{ v \in \mathcal N \mid S(\mathcal A) = 2^v\}.$$

**Lemma**: Let $\mathcal A$ be a class of sets of finite dimension $v$; then for any $n \geq v$. Then for $n \geq v$,
$$S_{\mathcal A}(n) = \sum_{j=1}^v \binom{n}{j} \leq (n+1)^v$$.

**Theorem (McDiarmid)**: Suppose $X_1, \dots, X_n$ are indepdendent RVs; put $Z = g(X_1, \dots, X_n)$. Then, if $g$ satisfies
$$\sup|g(X_1, \dots, X_i', \dots, X_n) - g(X_1, \dots, X_n)| \leq c_i$$
for some constants $c_i$, then 
$$P(Z - E[Z] > t) \leq e^{-2t^2 / \sum c_i^2}.$$

Applying this to the earlier result, we get that

Overall,
$$P(\mu(A) - \mu(A) > t) \leq e^{-2nt^2}.$$

Putting everything together gives us the bound
$$|\hat R_n(\hat h) - R(\hat h)| \leq 2 \sqrt{\frac{2V \log(n + 1) + \log 2}{n}} + \sqrt{\frac{\log(1/\delta)}{n}}$$
where $V$ is the VC-dimension.

But in practice, we need some sort of optimization procedure - and so we have to consider convex versions of the above. So we consider functions $f \in \mathcal F$ where $f: \mathcal X \to \R$ and a convex loss $\phi(Y_if(x_i))$.

**Lemma (Zhang)**: If $\phi$ is a convex function such that $\phi(0) = 1$ and $\phi(x) \geq \phi(-x)$ for $x \geq 0$, and 
$$H_\eta(\alpha) = \eta \phi(-\alpha) + (1-\eta)\phi(\alpha)$$
then for $\tau: [0, 1] \to \R$ satisfying $\tau(\eta) = \inf_\alpha H_\eta(\alpha)$, and $|0.5 - \eta| \leq c(1 - \tau(\eta))^\gamma$ for $\gamma \in (0,1)$, then
$$R(\text{sgn}(f)) - R^* \leq 2c(R_\phi(f) - R_\phi(f^*))$$
for all $f$.

**Theorem (Ledeux-Talagrand)**: If $\Psi: [-1, 1] \to \R$ is a contraction and satisfies $\Psi(0) = 0$, then for
$$\mathcal G = \{ \Psi \circ f, f \in \mathcal F \}$$
then
$$R_n(\mathcal G(z_i^n)) \leq 2R(\mathcal F(z_i^n))$$

Moreover, if $\psi$ is Lipschitz, then in fact we can use the above to see that
\begin{align*}
E[\sup_f |R_{n, \phi}(f) - R_{\phi}(f)|] &\leq E \left[ \sup_f \frac{1}{n} \sum_{i=1}^n \sigma_i(\phi(-Y_i f(x_i)) - \phi(-Y_i'f(x_i'))) \right] 
\end{align*}

We finally arrive at
$$E[R(\text{sgn}(\hat f)) - R^*] \leq 2c\left(8L \sqrt{\frac{2\log(2n)}{n}}\right)^\gamma + 2c\left(\inf_{\mathcal F} R_\psi(f) - R^*_\psi\right)^\gamma$$
where $L$ is the Lipschitz constant on $\psi$, $\gamma$ is as in the Zhang lemma, and the root-log factor depends on the Rademacher complexity of the family of functions under consideration.

# Multiple Hypothesis Testing

-----------------------

For ease, remember that a type I error is a false positive and a type II error is a false negative.

There are many possible approaches here. We give a few.

## The Global Null

Here we test the global null
$$H_0 = \bigcap_{j=1}^n H_{0, j}$$
where we just do some $p$-value adjustment.

1. The most obvious is Bonferonni, which proceeds simply by a union bound
$$P_{H_0}(\text{Type I Error}) = P_{H_0}\left(\min_{i=1,\dots,n} p_i \leq \frac{\alpha}{n}\right) \leq \sum_{i=1}^n P_{H_0}\left(p_i \leq \frac{\alpha}{n}\right) = \alpha.$$
Note that this bound isn't actually that conservative under independence:
$$P_{H_0}(\text{Type I Error}) = 1 - \left( 1 - \frac{\alpha}{n} \right)^n = 1 - e^{-\alpha + o(1)} \approx \alpha$$
for reasonably small $\alpha$. Note that this only cares about large deviations from then null.
2. We can also use the Fischer combination test, which rejects for large values of
$$T = -2 \sum_{i=1}^n \log(p_i).$$
This test requires independence between $p$-values, in which case $T \sim \chi^2_{2n}$.

The next thing to think about is the optimality of Boneferonni under sparse alternatives. Consider a Gaussian sequence model
$$Y_i \sim N(\mu_i, 1)$$
for $i \in [n]$ with the null $H_0: \mu_i = 0$ for all $i$ and alternative $\mu_i \neq 0$ for some $i$. Then, we consider the Bonferonni test
$$\max_{i \in [n]} y_i > |z(\alpha / n)|.$$

To examine the power, we have the result
$$|z(\alpha / n)| = \sqrt{2\log(n)}(1 + o(1)) \simeq \sqrt{2 \log(n)} \left( 1 + \frac{\log(\log(n))}{4n} \right)$$
and an approximation
$$|z(\alpha / n)| = \sqrt{B(1 - \log(B)/B)}$$
where $B = 2\log(n / \alpha) - \log(2\pi)$.

WLOG, suppose that $\mu_1 \neq 0$.

1. If $\mu_1 = (1 + \epsilon)\sqrt{2 \log(n)}$, then
$$P(\max Y_i > |z(\alpha / n)|) \geq P(Y_i > |z(\alpha / n)|) = P(z + (1 + \epsilon)\sqrt{2\log(n)} > |z(\alpha / n)|) \to 1$$
as $n \to \infty$.
2. On the other hand, by a similar argument we see that the power goes to $1 - \alpha$ with $\mu_1 = (1 - \epsilon)\sqrt{2\log(n)}$.

We can do the same analysis for the Fischer test. In this case, the test statistic $T_n = \sum_{i=1}^n Y_i^2$ is non-central $\chi^2_{n}(\|\mu\|^2)$. The CLT yields that asymptotically
$$\frac{T_n - (n + \|\mu\|^2)}{\sqrt{2n + 4\|\mu\|^2}} \sim N(0, 1).$$

Sime's modification of the Bonferonni is given by rejecting when there is some $i$ such that
$$p_{(i)} \leq \frac{\alpha \cdot i}{n}$$
rather than
$$p_{(1)} \leq \frac{\alpha}{n}$$
where $p_{(i)}$ are the ordered $p$-values.

Now under $H_0$ and independence of the $p$-values, then we have that
$$T_n = \min \frac{p_{(i)} n}{i} \sim U([0, 1])$$
which can be easily shown via induction.

We can test based on the empirical CDF as well.

_Def_: The **empirical CDF** is given by
$$\hat F_n(t) = \frac{|\{i \mid p_i \leq t\}|}{n} \sim_{H_0} \frac{1}{n}\text{Bin}(n, t).$$

We can consider the Kolmogorov Smirnov statistics
$$KS = \sup_{t \in [0, 1]} |\hat F_n(t) - t| \text{ and } KS^+ = \sup_{t \in [0, 1]} (\hat F_n(t) - t)$$
which, under independence, has distribution under the null satisfying
$$P(KS^+ \geq u) \leq e^{-2nu^2}.$$
Alternatively, just bootstrap the distribution of the KS-statistics under the null or something.

Alternatively, we consider the Anderson-Darling test, which computes
$$A^2 = n\int w(t)(\hat F_n(t) - t)^2 dt$$
for some weight function $w(t)$; common choices are $w(t) = 1$ (called the Cramer Von-Mises statistic) and $w(t) = t^{-1}(1-t)^{-1}$, which correponds to the expectation of $(\hat F_n(t) - t)^2$ under the null. Note that for this particular $w$, you get
$$A^2 = -n + \sum_{k=1}^n \frac{2k-1}{n} [\log(p_{(k)}) + \log(1 - p_{(k)})]$$

Lastly, we have the Tukey higher criticism test, which is given by the statistics
$$HC^* = \max_{0 \leq \alpha \leq \alpha_0} \frac{\hat F_n(\alpha) - \alpha}{\sqrt{\frac{\alpha(1-\alpha))}{n}}}$$

To understand the properties of these tests, often people will use a setup that looks something like
\begin{align*}
H_{0}&: x \sim N(0, 1) \\
H_{1}&: x \sim N(\mu, 1) \\
\end{align*}
for some proportion $\epsilon$ of non null results; in fact the higher criticism test (for unknown $\epsilon, \mu$ almost as good as the likelihood ratio test (for known $\epsilon, \mu$).

## The Family Wise Error Rate

For convenience, adopt the notation
\begin{align*}
  U &= |\{\text{true negatives}\}| \\
  V &= |\{\text{false positive}\}| \\
  T &= |\{\text{false negatives}\}| \\
  S &= |\{\text{true negatives}\}| \\
  R &= V + S \\
  n_0 &= U + V \\
  n &=  U + V + T + S 
\end{align*}

_Def_: The $k$**-family wise error rate** is the probability of $k$ type I errors.

As before, the simpliest possible thing is Bonferonni, e.g. reject if $p_i \leq \frac{\alpha}{n}$. Of course, this is very conservative. Under independence of the null $p$-values, you can get a bound of  $p_i \leq \frac{\alpha}{n(1 - \alpha / 2)}$ which is very close to Bonferonni.

_Def_ A testing procedure controls the FWER **weakly** if it controls it under the global null.

Consider the following LSD procedure from Fisher: test first the global null and then test each hypothesis at level $\alpha$.

Consider now Holm's procedure, which sorts the $p$-values $p_{(1)}, \dots, p_{(n)}$; then, you do a step-down where you reject if $p_{(i)} \leq \frac{\alpha}{n - i + 1}$ for all $i$ up to some bound and accept otherwise. This in fact controls the FWER strongly.

The closure principle is the following: suppose you have $\{H_1, \dots, H_n\}$; then the closure of these sets is $\{H_I \mid I \in P([n])\}$, where for any subset $I$, $H_I = \bigcap_{i \in I}H_i$. Now if you have some procedures $\phi_i$ that satisfy $P(\phi_I = 1 \mid H_I) \leq \alpha$, then we should reject $H_I$ if and only if for all $J$ such that $I \subset J$, $H_J$ is rejected at level $\alpha$ by $\phi_J$.

In fact, Holm is simply the closure of Bonferonni.
