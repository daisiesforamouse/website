---
title: 'Brownian Motion'
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

## Definition

-------

Fix a probability space $(\Omega, \mathcal F, P)$; we characterize the Brownian motion $\{B_t\}_{t \geq 0}$ via the following properties:

- **Independent Increments**: If $s < t$, the random variable $B_t - B_s$ is independent of $\sigma\{B_r: r \leq s\}$
- **Stationary Increments**: If $s < t$, then $B_t - B_s$ has the same distribution as $B_{t-s} - B_0$.
- **Continuity**: The map $t \mapsto B_t$ is almost surely continuous. 

**Theorem**: If a process satisfies the above, then there are $\mu, \sigma^2$ (respectively called the drift and the variance parameter, and $\sigma$ is named the volatility) such that there exist $B_t \sim N(\mu t, \sigma^2 t)$.

_Def_: A stochastic process $\{B_t\}_{t \geq 0}$ is called a (one dimensional) **Brownian motion** (or Wiener process) starting from the origin with drift $\mu$ and variance parameter $\sigma^2$ if $B_t = 0$ and the above three conditions are satisfied, with the imposition that 
$$
B_t - B_s \sim N(\mu (t-s), \sigma^2 (t-s)).
$$

**Prop**: If $B_t$ is a Brownian motion with $\mu = 0, \sigma^2 = 1$ (a so-called **standard Brownian motion**), then $Y_t = \sigma B_t + \mu t$ is a Brownian motion with parameters $\mu, \sigma^2$.

_Proof_: Obvious.

## Construction

-------

Pick a probability space $(\Omega, \mathcal F, P)$ that is rich enough to support a countable collection of independent standard normal variables. If you are particular, the unit interval with Lesbegue measure is sufficient here.

The strategy is as follows: we define $B_t$ for a countable dense set (in particular the dyadic rationals) of times using our precession of standard normals. then, we find some $t \mapsto B_t$ that agrees on the dense set and is uniformly continuous and then extend by continuity.

Set $D_n = \left\{ \frac{k}{2^n}, k = 0, 1, \dots, 2^n \right\}$ and $D = \bigcup_{n=0}^\infty D_n$; index our standard normals by $\{N_{q}\}_{q \in D}$, and set $B_0 = 0, B_1 = N_1$, and $B_{1/2} = \frac{B_1 - B_0}{2} + \frac{1}{2}N_{1/2}$. Just continue the same thing for every such dyadic, such that
$$
\{B_{1/2^n} - B_0, B_{2/2^n} - B_{1/2^n}, \dots, B_1- B_{(2^n-1) / 2^n} \}
$$
are all independent $N\left(0, 2^{-n}\right)$. 
