---
title: 'Mathematical Statistics I'
subtitle: 'UChicago STAT 30100, Winter 2024'
...

\newcommand{\N}{\mathbb{N}}
\newcommand{\R}{\mathbb{R}}
\newcommand{\Q}{\mathbb{Q}}
\newcommand{\Z}{\mathbb{Z}}
\newcommand{\C}{\mathbb{C}}

\DeclareMathOperator{\Var}{Var}
\DeclareMathOperator{\Cov}{Cov}
\DeclareMathOperator{\RSS}{RSS}
\DeclareMathOperator{\Jac}{Jac}
\DeclareMathOperator{\Ker}{Ker}
\DeclareMathOperator{\Im}{Im}
\DeclareMathOperator{\argmin}{argmin}
\DeclareMathOperator{\argmax}{argmin}
\DeclareMathOperator{\proj}{proj}
\DeclareMathOperator{\rank}{rank}

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

## Sufficient Statistics

-------

We model the outcome of a statistical experiment by some $P_\theta$, where $\theta \in \Theta$ is a parameter space; we take some $X_1, \dots, X_n \sim P_\theta$, and we ask: can we summarize $X_1, \dots, X_n$ by some statistic $T = T(X_1, \dots, X_n)$ without losing information?

_Def_: [Link](https://en.wikipedia.org/wiki/Sufficient_statistic) A **sufficient statistic** is a statistic $T = T(X_1, \dots, X_n)$ whose conditional distribution $X \mid T$ does not depend on $\theta$.

The intuition is that a sufficient statistic recovers all information about $\theta$. In fact, if we sample $Y_1 \mid T, \dots, Y_n \mid T$, then $Y_1, \dots, Y_n$ has the same joint distribution as $X_1, \dots, X_n$.

_Def_: Let $X_1, \dots, X_n \sim N(\theta, 1)$; then $T = n^{-1} \sum_{i=1}^n X_i = \bar X$ is a sufficient statistic. In fact,
$$
    \bmat{X_1 \\ \vdots \\ X_n} \mid T \sim N \left( \bmat{ \bar X \\ \vdots \\ \bar X}, I - n^{-1} U \right)
$$
where $U$ is a matrix with ones on upper triangular entries.

_Def_: Let $X_1, \dots, X_n \sim P_\theta$ for some $\theta \in \Theta$; the **order statistic** is just $T = (X_{(1)}, \dots, X_{(n)})$ where $X_{(1)} \leq \dots \leq X_{(n)}$ is the ordered list of observations. This is a sufficient statistic when $X_1, \dots, X_n$ are exchangeable.

**Theorem (Factorization)**: Suppose $P_\theta$ is either discrete or continuous; then $T = T(X)$ is sufficient if and only if $p(X \mid \theta) = g_\theta(T(x)) h(x)$ where $p$ is the density function corresponding to $P_\theta$.

_Proof_: The continuous case is similar to the discrete one. The backwards direction is just computation. The forwards direction is just noting that 
$$ P(X = x) = P(X = x, T(X) = T(x)) = P(X = x \mid T(X) = T(x)) P(T(X) = T(x)). $$

_Corollary_: $T = T(X)$ is sufficient iff $\theta \to T \to X$ is a Markov chain, e.g. $\theta \perp X \mid T$.

_Def_: An **exponential family** is given by the distribution
$$
P(X \mid \theta) = \exp \left( \sum_{j=1}^d \eta_j(\theta)T_j(x) - B(\theta) \right)h(x).
$$
Moreover,

- $\eta_j$ is the natural parameter,
- $T_j$ is the sufficient statistic,
- $d$ is the dimension,
- $h$ is the base measure,
- and
$$ 
B(\theta) = \log \left( \int \exp\left(\sum_{j=1}^d \eta_j(\theta) T_j(x)\right) h(x) d\mu \right) 
$$
is the log partition function (though here it is merely a normalizing constant).

_Def_: $P_\eta, \eta \in H$ is an **exponential family of canonical form** if 
$$
p(x \mid \eta) = \exp \left( \sum_{j=1}^d \eta_j T_j(x) - A(\eta) \right)h(x)
$$

_Def_: $P_\eta, \eta \in H$ is a **minimal** exponential family (of canonical form) if the dimension cannot be further reduced.

_Def_: A minimal canonical exponential family $P_\eta, \eta \in H$ is **curved** if the natural parameters are nonlinearly related.

_Def_: A statistic $S$ is minimal sufficient if it is sufficient and for any sufficient $T$, $S$ is a function of $T$.

### Subfamily Methods


**Lemma**: Suppose $\Theta_0 \subset \Theta$. If $T$ is sufficient for $\theta \in \Theta$ and it is minimal sufficient for $\theta \in \Theta_0$, then $T$ is minimal sufficient for all $\theta \in \Theta$.

_Proof_: That $T$ is sufficient comes by assumption. Minimality holds since on $\Theta_0$, $S(x) = f(T(x))$ for some function $f$, and this automatically extends to $\Theta$.

**Theorem**: Let $P_\theta, \theta \in \{ \theta_0, \dots, \theta_d \}$ have the same support. Then 
$$
T = \left(  \frac{p(x \mid \theta_1)}{p(x \mid \theta_0)}, \dots, \frac{p(x \mid \theta_d)}{p(x \mid \theta_0)} \right)
$$
is minimal sufficient for $\theta$.

_Proof_: We have that
$$
    p(x \mid \theta_j) = T_j(x) p(x \mid \theta_0)
$$
so we are done by factorization.

**Theorem**: Let $(P_\eta, \eta \in H)$ be a minimal exponential family, i.e.
$$
    p(x \mid \eta) = \exp \left( \langle \eta, T(x) \rangle - A(\eta) \right)h(x).
$$
Then $T(x) \in \R^d$ is minimal sufficient.

_Proof_: We wish to find $\eta_0, \dots, \eta_d \in H$ such that the matrix $A = \bmat{(\eta_i - \eta_{i-1})^T} \in \R^{d \times d}$ is full rank. This possible in both the full rank and curved cases (draw a picture!). But we know that the likelihood ratios are minimal sufficient for the subfamily and since
$$
    \frac{p(x \mid \eta_j)}{p(x \mid \eta_0)} = \exp \left(\langle \eta_j - \eta_0, T(x) \rangle - (A(\eta_j) - A(\eta_0))\right)
$$
is just a function of $\langle \eta_j - \eta_0, T \rangle$. Then, $AT$ is minimal sufficient, and since $A$ is invertible, $T$ is minimal sufficient.

### Completeness Methods


Let $X_1, X_2$ be i.i.d. $N(\mu, 1)$. Then $(X_1 - X_2, X_1 + X_2)$ is a one-to-one transform, but $X_1 - X_2 \sim N(0, 2)$ is unrelated to $\mu$, so $X_1 + X_2$ is the only one that matters.

_Def_: $A = A(x)$ is **ancillary** if its distribution does not depend on $\theta \in \Theta$. $A$ is **first-order ancillary** if $E_\theta[A(x)]$ does not depend on $\theta \in \Theta$.

_Def_: $T = T(x)$ is a complete statistic if the following holds:
$$
    E_\theta[f(T(X))] = 0, \forall \theta \in \Theta \implies f(T(X)) = 0, a.s. \ \forall \theta \in \Theta.
$$
Morally, this means that there is no nontrivial transform that makes $T$ first-order ancillary.

**Theorem (Bahadur)**: If $T$ is sufficient and complete, then $T$ is minimal sufficient.

_Proof_: Suppose some minimal sufficient statistic $U = U(x)$ exists. Since $U$ is minimal, there is $h$ such that $U = h(T)$. Define $g(u) = E_\theta[T \mid U = u]$, which by sufficiency does not depend on $\theta$. Then, 
$$
E_\theta[g(h(T))] = E_\theta[g(U)] = E_\theta[T]
$$
so $E_\theta[g(h(T)) - T] = 0$ for all $\theta$. Since $T$ is complete, then $g(h(T)) = g(U) = T$ almost surely.

**Theorem (Basu)**: If $T$ is sufficient and complete, and $A$ is ancillary, then $T \perp A$.

_Proof_: We want to show that $P_\theta(A \in B \mid T) = P_\theta(A \in B)$. Let 
$$
g(t) = P_\theta(A \in B \mid T = t)
$$
which does not depend on $\theta$ by sufficiency, and similarly $P_\theta(A \in B \mid T = t)$ is also just a constant, now abbreviated to $c$. Then, $E[g(t) - c] = 0$, and by completeness $g(t) = c$ so we win.

**Theorem**: If $P_\eta(x) = \exp(\eta^\top T(x) - A(\eta))$ is a full-rank exponential family, then $T$ is complete.

_Proof_: Compare moment generating functions.

## Decision Theory

-------

_Def_: If we have some family $(P_\theta, \theta \in \Theta)$, and some observations $X \sim P_\theta$, then we may call any $\hat \theta = \hat \theta(X)$ a **decision/procedure/estimator**. We measure how good the estimator is using a **loss function** $L(\theta, \hat \theta)$; we also call the expectation $E_\theta[L(\theta, \hat \theta)]$ the **risk**.

_Def_: An estimator $\hat \theta$ is **minimax** if for any $\tilde \theta$,
$$
    \sup_{\theta \in \Theta} R(\hat \theta, \theta) \leq \sup_{\theta \in \Theta} R(\tilde \theta, \theta).
$$

_Def_: An estimator $\hat \theta$ is **Bayes** with respect to a prior distribution $\pi$ on $\Theta$ if for any $\tilde \theta$,
$$
    \int R(\hat \theta, \theta)\pi(\theta)d\theta \leq \int R(\tilde \theta, \theta)\pi(\theta)d\theta.
$$

**Theorem**: For any $\hat \theta$ and any sufficient statistic T, if $L(\hat \theta, \theta)$ is convex in $\hat \theta$, then $\tilde \theta = E_\theta[\hat \theta \mid T]$ must have no greater risk.

_Proof_: Apply Jensen to $L(\tilde \theta, \theta)$, and take a second expectation.

For the Bayes estimator, we can also frame it in terms of the marginal distribution of $x$, e.g.
$$
    \int R(\hat \theta, \theta) \pi(\theta) d\theta = \iint L(\hat \theta(x), \theta) p(x \mid \theta)\pi(\theta) dxd\theta = \iint L(\hat \theta(x), \theta), \pi(\theta \mid x) m(x)d\theta dx
$$
where $\pi(\theta \mid x), m(x)$ are determined by Bayes.

**Lemma**: Define 
$$
\hat \theta_\pi(x) = \argmin_a \int L(a, \theta) \pi(\theta \mid x)d\theta.
$$
Then $\hat \theta_\pi$ is Bayes.

_Example_: If $L(\hat \theta, \theta) = (\hat \theta - \theta)^2$, then the Bayes optimizer is just the posterior mean $E[\theta \mid X]$.

**Theorem**: Suppose that there is some prior $\pi$ such that the Bayes estimator satisfies
$$
    \sup_{\theta \in \Theta} R(\hat \theta, \theta) = \inf_{\tilde \theta} \int R(\hat \theta, \theta) \pi(\theta) d\theta.
$$
Then, $\hat \theta$ is minimax.

_Proof_: For all $\hat \theta$, 
$$
\sup_{\theta \in \Theta} R(\tilde \theta, \theta) \geq \int R(\tilde \theta, \theta) \pi(\theta) d\theta \geq \inf_{\hat \theta} \int R(\tilde \theta, \theta)\pi (\theta)d\theta = \sup_{\theta \in \Theta} R(\hat \theta, \theta).
$$

**Corollary**: Suppose $\hat \theta$ is Bayes with risk independent of $\theta$. Then it is also minimax.

**Theorem**: Suppose there exists some sequence of prior distributions $\pi_n$ such that 
$$
    \sup_{\theta \in \Theta} R(\hat \theta, \theta) = \lim_{n \to \infty} \inf_{\tilde \theta} \int R(\tilde \theta, \theta) \pi_n(\theta) d\theta.
$$

_Proof_: Exactly the same as above.

### Admissible Estimators

_Def_: We say $\hat \theta$ is **inadmissible** if there exists some other estimator $\tilde \theta$ such that
$$
    R(\tilde \theta, \theta) \leq R(\hat \theta, \theta), \ \ \forall \theta \in \Theta
$$
and there is some $\theta_0 \in \Theta$ such that
$$
    R(\tilde \theta, \theta_0) \leq R(\hat \theta, \theta_0).
$$
Otherwise, $\hat \theta$ is called **admissible**.

**Theorem**: If $\hat \theta$ is Bayes, then it is also admissible.

_Proof_: If $\hat \theta$ was inadmissible, then there is some $\tilde \theta$ satisfying the above conditions. Then, put $\tilde \theta$ in the condition for Bayes optimality and win (as long as $R(\tilde \theta, \theta_0) \leq R(\hat \theta, \theta_0)$ holds on a set of positive measure w.r.t. the prior distribution).

**Theorem (Complete Class, Brown)**: If $\hat \theta$ is admissible, then $\hat \theta_{\pi_n} \to \hat \theta$.

### Nonparametric Density Estimation

Suppose we have $X_1, \dots, X_n \sim f$, where $f \in S_\alpha(R) \subset L^2[0, 1]$ and $S_\alpha(R)$ is the Sobolev ball
$$
    S_\alpha(R) = \left\{ f \in L^2[0, 1] \mid f \geq 0, \int f = 1, f = \sum_{i=1}^\infty \theta_j \phi_j(x), \sum_{i=1}^\infty j^{2\alpha}\theta_j^2 \leq R^2\right\}.
$$
We will consider the loss function 
$$
L(\hat f, f) = \|\hat f-f\|^2 = \int_0^1 ( \hat f - f)^2
$$
and proceed with Fourier analysis, e.g. expressing
$$
    f(x) = a_0 + \sum_{j=1}^\infty \left(a_i \sin(2\pi j x) + b_j \cos(2\pi jx)\right) = \sum_{j=1}^\infty \theta_j \phi_j(x).
$$

The natural choice for estimating $\theta_j$ is the empirical coefficient
$$
    \hat \theta_j = \frac{1}{n} \sum_{i=1}^n \phi_j(X_i) \to E[\phi_j] = \int_0^1 f(x)\phi_j(x)dx = \theta_j.
$$
which has approximate distribution $N\left(\theta_j, \frac{\Var(\phi_j(X_i))}{n}\right)$.

### Gaussian Sequence Model

Suppose that we observe the following sequence $X_j = \theta + \frac{1}{\sqrt{n}} Z_j$, where $j \in \N$; set a loss function 
$$
L(\hat \theta, \theta) = \|\hat \theta - \theta\|^2 = \sum_{j=1}^\infty (\hat \theta_j - \theta_j)^2.
$$

This is the same as before (for details, look for asymptotic equivalence); in fact this is asymptotically equivalent to white noise models and nonparametric regression.

In this case, we may take $\hat \theta_j = X_j$ if $j \leq k$ and $0$ otherwise; this gives an optimal loss rate of $n^{-\frac{2\alpha}{2\alpha + 1}}$. We can do better however, since this relies on setting $k = n^{\frac{1}{2\alpha+1}}$, which relies on the smoothness parameter.

We can overcome this by an adaptive nonparametric estimation; we perform a blockwise James-Stein estimator for blocks of size $\frac{2}{3} 3^j$, up to $\log_3(n)$ blocks.

In fact this is rate-optimal: in general
$$
    \inf_{\hat \theta} \sup_{\theta \in \Theta_{\alpha}(R)} E_{\theta} \|\hat \theta - \theta\|^2 = (1 + o(1)) C_{\alpha, R}n^{-\frac{2\alpha}{2\alpha + 1}}
$$
for some Pinsker constant $C(\alpha, R)$. In fact the blockwise James-Stein achieves this constant as well.

We will prove a weaker version of this.

_Def_: In general, for any hypothesis test with
$$
\begin{align*}
    H_0: X \sim P \\
    H_1: X \sim Q \\
\end{align*}
$$
a **test** is a measurable $\phi: X \to \{ 0, 1 \}$; a type-1 error happens with probability $E_P[\phi] = P[\phi]$ and a type-2 with $E_Q[1 - \phi] = Q[1 - \phi]$ (we will generally use this shorthand for the expectation). The testing error is their sum, and the **optimal testing error** is
$$
    \inf_\phi P[\phi] + Q[1 - \phi].
$$
Moreover, the **total variation** is
$$
    TV(P, Q) = \sup_{B} |P[B] - Q[B]|
$$
where the supremum is over all events.

**Theorem**: 
$$
    TV(P, Q) = \frac{1}{2} \int |p  - q| = 1 - \int p \land q
$$
where the integral is w.r.t. a dominating measure, and $p, q$ are their densities w.r.t. that dominating measure.

_Proof_: Just look with $B = \{p(x) > q(x)\}$.

_Def_: We call $\int p \land q$ the affinity.

**Theorem (Neyman-Pearson)**: The optimal testing error is
$$
    \inf_{\phi} (P[\phi] - Q[1 - \phi]) = 1 - TV(P, Q) = \int p \land q.
$$
Note that this is a maximum likelihood estimator.

_Proof_: Same as above, but use the earlier theorem.

**Theorem (LeCam)**: For all $\theta_0, \theta_1 \in \Theta$,
$$
    \inf_{\theta} \sup_{\theta \in \Theta} E_{\theta}[(\hat \theta - \theta)^2] \geq \frac{1}{4}(\theta_0 - \theta_1)^2 \int p_{\theta_0} \land p_{\theta_1}.
$$

_Proof_: Bound by the average.
$$
    \inf_{\hat \theta} \sup_{\theta \in \Theta} E_{\theta}[(\hat \theta - \theta)^2] \geq \frac{1}{2} \inf_{\hat \theta} \left( E_{\theta_0}[(\hat \theta - \theta_0)^2] + E_{\theta_1}[(\hat \theta - \theta_1)^2] \right).
$$
It is not difficult to finish from here.

<!-- To know that the blockwise JS estimator is rate-optimal, apply LeCam to $\theta_0 = 0, \theta-1 = \frac{1}{\sqrt{n}}$. -->

## Estimation Under Constraints

-------

_Def_: We say $\delta(X)$ is **UMVUE (uniformly minimum-variance unbiased estimator)** for $g(\theta)$ if it is

- unbiased;
- has lower variance than any other unbiased estimator.

**Theorem (Lehmann-Scheffe)**: Let $T = T(X)$ be a sufficient and complete statistic. Then, $\delta(X) = h(T(X))$ is the only function of $T$ that is unbiased and is also the unique UMVUE, so long as $E_\theta[\delta(X)] = g(\theta)$.

_Proof_: If $\tilde h(T(X))$ is also unbiased,
$$
    E[h(T(X)) - \tilde h(T(X))] = g(\theta) - g(\theta) = 0
$$
and by completeness $h = \tilde h$ almost surely. Suppose $\tilde \delta (X)$ is unbiased; then
$$
  \delta(X) = E_\theta[\tilde \delta(X) \mid T(X)]
$$
is also unbiased, but is a function of $T$ and has a lower variance than $\tilde \delta(X)$. By uniqueness from the first part, it must be that $\delta(X) = h(T(X))$ and is thus UMVUE. For uniqueness, look at $\Var_\theta(\tilde \delta)$, and note that it is strictly larger than the variance of $\delta(X)$, unless $\tilde \delta = \delta$ almost surely.

**Theorem (SURE)**: For $X \sim N(\theta, \sigma^2 I_p)$, 
$$
    \|\hat \theta - X\|^2 - \sigma^2 p + 2\sigma^2 \sum_{j=1}^p \frac{\partial \hat \theta_j}{ \partial X_j}
$$
is an unbiased estimator of $E_\theta[(\hat \theta - \theta)^2]$. We also call 
$$
\sum_{j=1}^p \frac{\partial \hat \theta_j}{\partial X_j}
$$
the **degrees of freedom**.

_Proof_: Apply Stein's lemma after introducing a $X - X$ term in the square.

## Maximum Likelihood Theory

-------

In this section, $\theta^*$ always denotes the "true" parameter.

_Def_: Given a likelihood $p_\theta(x)$ and $x_1, \dots, x_n$, the **maximum likelihood estimator (MLE)** is given by
$$
\hat \theta = \argmax_{\theta \in \Theta} \sum_{i=1}^n \log p_\theta(x_i).
$$

We wish to establish some sort of consistency (consider the case of $X_1, \dots, X_n$ Cauchy variables for instance). We have that
$$
\hat \theta = \argmin \frac{1}{n} \sum_{i=1}^n \log\left(\frac{p_{\theta^*}(x_i)}{p_\theta(x_i)}\right)
$$
but by the LLN we have
$$
\frac{1}{n} \sum_{i=1}^n \log \left(\frac{p_{\theta^*}(x_i)}{p_\theta(x_i)}\right) \to \int p_{\theta^*} \log \left( \frac{p_{\theta^*}}{p_\theta} \right).
$$

_Def_: We call 
$$
D(P \| Q) = \int p \log\left(\frac{p}{q}\right) = \int p \log(1/q) - \int p \log(1/p)
$$
the **Kullback-Leibler divergence** and $\int p \log(1/p)$ is the (relative) **entropy**.

**Prop**: The KL divergence is always positive and is in fact bounded below by the square Hellinger ditance
$$
D(P \| Q) \geq \frac{1}{2} \int (\sqrt{p} - \sqrt{q})^2.
$$

Then, the above gives the intuition that
$$
\hat \theta = \argmin_{\theta} D(P_{\theta^*} \| P_\theta).
$$

More formally, we require a constraint.

- **Identifiability**: If $\theta \neq \theta^*$, then $D(P_{\theta^*} \| P_\theta) > 0$. In fact, we must have that $\forall \epsilon > 0$, there is $\delta > 0$ such that $|\theta - \theta^*| > \epsilon$ \implies $D(P_{\theta^*} \| P_\theta) > \delta$.

**Theorem**: In this case, $\hat \theta \to \theta^*$ in probability.

_Proof_: We first compute that
$$
P_{\theta^*}(|\hat \theta - \theta^*| > \epsilon) \leq P_{\theta^*}(D(P_{\theta^*} \| P_{\hat \theta}) > \delta) = P_{\theta^*} \left(\int p_{\theta^*}\log(p_{\theta^*}) - \int p_{\theta^*}\log(p_{\hat \theta}) > 0\right).
$$
Then, we need a uniform LLN, e.g.
$$
\sup_{\theta} \left| \frac{1}{n} \sum_{i=1}^n \log(p_\theta(x_i)) - \int p_{\theta^*}\log(p_\theta) \right| \to 0
$$
in probability under $p_{\theta^*}$. Then approximate the two integrals by the suitable LLN and win.

A proof of the ULLN is omitted, but it holds when the log-densities are low-complexity, e.g.
$$
\sup_{f \in \mathcal F} \left|\frac{1}{n} \sum_{i=1}^n f(x_i) - E[f(X)] \right| \to 0
$$
holds given a condition on the VC dimension/covering number/Dudley integral/Rademacher complexity etc.

_Def_: We define the **score** as
$$
S_\theta(x) = \frac{\partial}{\partial \theta}\log(p_\theta(x))
$$
and the **Fisher information** as
$$
I_\theta = E_\theta[(S(X))^2].
$$

**Theorem**: Suppose that

- $\hat \theta \to \theta^*$ in probability;
- we have
$$
\log(p_{\theta^* + t}(x)) = \log(p_{\theta^*}(x)) + tS_{\theta^*}(x) + |t|r(x,t)
$$
where
$$
\sup_{t \in U_n} \frac{|\nu_n r(\cdot, t)|}{1 + \sqrt{n}|t|} \to 0
$$
in probability and $\nu_n f = \frac{1}{\sqrt{n}} \sum_{i=1}^n(f(X_i) - E[f(X_i)])$ is stochastically bounded;
- and finally
$$
\mathcal L(\theta^* + t) = E_{\theta^*}[L_n(\theta^* + t)] = \mathcal L(\theta^*) - \frac{1}{2}I_{\theta^*}t^2 + o(t^2)
$$
where $I_{\theta^*} = \Var_{\theta^*}(S_{\theta^*}(X))$.

Then the MLE $\hat \theta$ is asymtotically normal:
$$
  \sqrt{n}(\hat \theta - \theta^*) \to N(0, I^{-1}_{\theta^*}).
$$

_Proof_: We first establish a quadratic expansion for $\mathcal L_n(\theta) = \frac{1}{n} \sum_{i=1}^n \log(p_\theta(X_i))$. In particular,
$$
\mathcal L_n(\theta^* + t) = \mathcal L(\theta^* + t) + \frac{1}{n} \nu_n \log(p_{\theta^*}) = L_n(\theta^*) + \frac{t}{\sqrt{n}}\nu_n S_{\theta^*} - \frac{1}{2}I_{\theta^*}t^2 + o(t^2) + \frac{|t|}{\sqrt{n}}r(\cdot, t).
$$
Now we want that $|\hat \theta - \theta^*| = O_P(n^{-1/2})$. Set $\hat t = \hat \theta - \theta^*$; by definition of the MLE, it must be that
$$
\mathcal L_n(\theta^* + \hat t) \geq \mathcal L_n(\theta^*)
$$
and by the previous expansion, we get that
$$
\frac{1}{2}I_{\theta^*} \hat t^2 - o(\hat t^2) \leq O_p(|\hat t|n^{-1/2}) +|\hat t|n^{-1/2}o_P(1 + \sqrt{n}|\hat t|) = O_p(|\hat t|n^{-1/2})
$$
so $|\hat \theta - \theta^*| = O_P(n^{-1/2})$ holds. Then we may control the remainder to be on the order of $n^{-1}$, and use
$$
\mathcal L_n(\theta^* + \hat t) \geq \mathcal L_n\left(\theta^* + n^{-1/2}\nu_n \frac{S_{\theta^*}}{I_{\theta^*}}\right)
$$
to finally get
$$
\frac{1}{2}I_{\theta^*} \left|\sqrt{n}\hat t - \nu_n \frac{S_{\theta^*}}{I_{\theta^*}}\right| = o_P(1)
$$
and so
$$
\sqrt{n} \hat t = \nu_n \frac{S_{\theta^*}}{I_{\theta^*}} + o_P(1) \to N(0, I^{-1}_{\theta^*})
$$
by the CLT.

_Def_: We say a distribution $P_\theta$ is **LAN (local asymptotically normal)** if
$$
\log \left(\prod_{i=1}^n \frac{p_{\theta^*+hn^{-1/2}}(X_i)}{p_{\theta^*}(X_i)}\right) = hn^{-1/2}\sum_{i=1}^n S_{\theta^*}(X_i) - h^2I_{\theta^*} + o_P(1)
$$
for all nonrandom constant $h$.

_Def_: We say that a function $p_{\theta}$ is **DQM (differentiable under quadratic mean)** if
$$
\int \left(\frac{\sqrt{p_{\theta+t}(X)} - \sqrt{p_\theta(X)}}{t} - \frac{1}{2}\sqrt{p_{\theta}} S_\theta \right)^2 \to 0
$$

**Theorem (LeCam)**: DQM implies LAN.

**Corollary**: The second condition above implies the third.

### Optimality of the MLE

**Theorem (Cramer-Rao)**: Suppose we have an unbiased estimator $\hat \theta$ for $X_1, \dots, X_n \sim P_{\theta^*}$; then $\Var_\theta(\hat \theta) \geq (nI_\theta)^{-1}$.

_Proof_: We have by definition
$$
\Var_\theta(\hat \theta) I_\theta = E_\theta[(\hat \theta - \theta)^2]E_\theta[S_\theta^2] \geq \left(\int (\hat \theta - \theta) \cdot p_\theta\frac{\partial}{\partial \theta} \log(p_\theta)\right)^2 \geq 1.
$$

Note that the above bound is true in finite sample sizes only; see Hodge's estimator for an example of super-efficiency. Note that James-Stein is actually another example.

_Def_: A **randomized statistic** of $X$ is some random variable $T(X, U)$ where $U$ is independent of $X$.

**Theorem (Hajek & LeCam)**: Suppose $(P_\theta \mid \theta \in \Theta)$ is DQM at $\theta$, and $T_n$ is some statistic for $(P_{\theta + h/\sqrt{n}} \mid h \in \R)$ satisfying
$$
  \sqrt{n} \left(T_n - (\theta + h/\sqrt{n})\right) \to L_{\theta, h}.
$$
Then there is some randomized statistic $T$ for $(N(h, I_\theta^{-1}) \mid h \in R)$ such that
$$
T - h \sim L_{\theta, h}.
$$

_Def_: A randomized statistic $T$ is **equivariant-in-law** if
$$
T -  h \sim L_\theta.
$$

**Theorem:** If $T$ is a randomized statistic for $(N(h, I_\theta^{-1}) \mid h \in \R)$ and is equivariant-in-law then there is a distribution $M_\theta$ such that
$$
T - h \sim N(0, I_\theta^{-1}) * M_\theta
$$
where the above is a convolution.

**Corollary**: Then we may always write such a statistic as $T = T(X, U) = X + U$.

**Corollary**: In the theorem of Hajek & LeCam, if we have
$$
  \sqrt{n} \left(T_n - (\theta + h/\sqrt{n})\right) \to L_{\theta}
$$
then the Cramer-Rao lower bound holds asymptotically as well.

But we can remove even the requirement of asymptotic equivalence-in-law.

**Lemma**: Suppose $P_\theta$ is DQM for all $\theta \in \Theta$ and $T_n$ is some statistic for $P_\theta^n$ such that
$$
\sqrt{n}(T_n - \theta) \to L_\theta.
$$
Then there is a subsequence such that Lebesgue almost everywhere $(\theta, h)$,
$$
\sqrt{n}(T_n - (\theta + h/\sqrt{n})) \to L_\theta.
$$

**Theorem (Almost Everywhere Convolution Theorem)**: Suppose $P_\theta$ is DQM for all $\theta \in \Theta$, $T_n$ is some statistic for $(P_\theta^n \mid \theta \in \Theta)$ such that for all $\theta$,
$$
\sqrt{n}( T_n - \theta) \to L_\theta;
$$
then $L_\theta = N(0, I_{\theta^*}) * M_\theta$ almost everywhere.

**Corollary**: Superefficiency can only hold on a set of measure zero.

**Theorem (Local Asymptotic Minimaxity)**: Suppose $(P_\theta \mid \theta \in \Theta)$ is DQM at $\theta \in \Theta$, $T_n$ is some statistic for $(P_{\theta + h/\sqrt{n}} \mid h \in \R)$, and $\ell(\cdot)$ is convex. Then
$$
\liminf_{n \to \infty} \sup_{h \in [-\delta_n, \delta_n]}  E^n_{\theta + h/\sqrt{n}} \left[ \ell(\sqrt{n}(T_n - (\theta + h/\sqrt{n}))) \right] \geq \int \ell d N(0, I_\theta^{-1}).
$$
for some $\delta_n \to 0$.

Can we achieve both lower bounds? Yes: assuming DQM, LeCam showed that there is a one-step improvement of the MLE that works. If we assume DQM and some empirical process conditions (as seen before) then the MLE works.
