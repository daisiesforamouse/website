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

## Basic Statistics

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
$$
P(X = x) = P(X = x, T(X) = T(x)) = P(X = x \mid T(X) = T(x)) P(T(X) = T(x)).
$$
