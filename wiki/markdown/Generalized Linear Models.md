---
title: 'Generalized Linear Models'
subtitle: 'UChicago STAT 34700, Winter 2024'
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

## Generalized Linear Models

-------

Suppose that we have data points $(X_1, y_1), \dots, (X_n, y_n)$; we operate under the framework that $X_i$ are always implicitly conditioned on. The generalized linear model will allow us to generalize from continuous real values to more diverse outputs. Particularly common is $y_i$ arising from an exponential family.

_Def_: A **link function** is how $E[y_i]$ (or $E[y_i \mid X_i]$) depends on $X_i$. In particular, we will have 
$$
    g(E[y_i]) = g(\mu_i) = X_i^T \beta.
$$ 
