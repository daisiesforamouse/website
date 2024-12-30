---
title: 'Basic Probability'
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
\DeclareMathOperator{\Cov}{Cov}
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

## Definitions

-------

_Def_: We say that a measurable space $(\Omega, \mathcal F, P)$ is a **probability space** if $P(\Omega) = 1$; in turn, a **random variable** $X$ is a measurable function $X: \Omega \to S$, for some measurable space $S$, usually either $\R^n$ or $\C^n$.

From now on, we will implicitly have all random variables which are mentioned together share the same probability space, and we will only consider real valued variables unless otherwise stated.

_Def_: The **expectation** of a random variable $X$ is simply the integral
$$
    E[X] = \int_\Omega X dP.
$$

_Def_: The **covariance** between two random variables $X$ and $Y$ is the expectation
$$
    \Cov(X, Y) = E[(X - E[X])(Y - E[Y])]
$$
and the **variance** is simply
$$
    \Var(X) = \Cov(X, X) = E[(X - E[X])^2].
$$

_Def_: A finite collection of events $E_1, \dots, E_n$ are **independent** if
$$
    P \left( \bigcap_{i=1}^m E_{k_i} \right) = \prod_{i=1}^m P(E_{k_i})
$$
for all $2 \leq m \leq n$ and $1 \leq k_i \leq n$. An infinite collection of events are independent if any finite subcollection of those events are independent.

_Def_: A collection of classes of sets $\{\mathcal F_n\}_{n=1}^\infty$ are **independent** if for any collection of events $\{E_n\}_{n=1}^\infty$ with $E_n \in \mathcal F_n$ are independent. Similarly, a sequence of random variables $\{X_n\}_{n=1}^\infty$ are **independent** if $\{\sigma(X_n)\}_{n=1}^\infty$ are independent.

_Def_: For a sequence of events $\{E_n\}_{n=1}^\infty$, we say that $E$ is the event which happens **infinitely often** (oft abbreviated i.o.) if 
$$
    E = \{E_n \text{ i.o.} \} = \bigcap_{m=1}^\infty \bigcup_{n=m}^\infty A_n = \limsup_n A_n.
$$
The counterpart is **eventually**, where
$$
    E = \{E_n \text{ eventually} \} = \bigcup_{m=1}^\infty \bigcap_{n=m}^\infty A_n = \liminf_n A_n.
$$

_Def_: Given a sequence of random variables $\{X_n\}_{n=1}^\infty$ as well as another random variable $X$, we say that

1. $X_n \to X$ **in probability** if, for all $\epsilon > 0$,
$$
    \lim_{n \to \infty} P(|X_n - X| > \epsilon) = 0.
$$
2. $X_n \to X$ **almost surely**, if
$$
    P(\lim_{n \to \infty} X_n = X) = 1.
$$
3. $X_n \to X$ **in** $L^p$, if 
$$
    E[|X_n - X|^p]^{1/p} \to 0.
$$

## Fundamental Results

-------

We begin with the three big guns for interchanging integration and limits.
