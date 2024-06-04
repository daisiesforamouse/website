---
title: 'Modern Methods in Applied Statistics'
subtitle: 'UChicago STAT 34800, Spring 2024'
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
\DeclareMathOperator{\Pr}{Pr}

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


## Supervised Learning and Decision Theory

-------

### Binary Classification

Suppose that we have i.i.d. input-output pairs $(X_1, Y_1), \dots, (X_n, Y_n) \sim \Pr(X, Y)$, where $\Pr$ is some (unknown) distribution from the population and $Y \in \{ 0, 1 \}$. We will then use $P(X, Y)$ to denote properties estimated from the sample.

If a new patient has (discrete-valued) symptoms $X = x$, we can compute the conditional estimates $P(X = x \mid Y = 0)$ and $P(X = x \mid Y = 1)$.

_Def_: In this case we define the **likelihood ratio** as
$$
LR = \frac{P(X = x \mid Y = 0)}{P(X = x \mid Y = 1)}.
$$

Moreover, we interpret this as how much more likely the observed symptoms are under $Y = 0$ than $Y = 1$.

In the case of continuous covariates, the idea remains the same, but now we assume that the continuous covariate is a member of some parametric family, e.g.
$$
\Pr(X \mid Y = y) = P_{\theta_y}(X)
$$
where $\theta_y$ is some parameter depending on $y$.

Of course, the base rate is also important; this is $\Pr(Y = 1)$ which is the distribution of $Y$ before seeing any other data.

**Theorem (Bayes)**: The prior and the posterior is related by the base rates:
$$
\Pr(Y\mid X) = \frac{\Pr(X\mid Y)\Pr(Y)}{\Pr(X)}
$$
where we call the denominator a normalizing constant (as it is constant in $Y$).

We can apply this to the likelihood ratio to get
$$
\text{posterior odds} = \frac{P(Y = 1 \mid X = x)}{P(Y = 0 \mid X = x)} = LR \cdot \frac{P(Y = 1)}{P(Y = 0)} = LR \cdot \text{prior odds}.
$$

_Def_. We say a function $\delta(x)$ is a **decision rule** when it maps input data to some label in $\{ 0, 1\}$. To formulate a decision rule, we utilize a **loss function**, which is any map $\ell(y, \hat y)$ sending a ground truth and a prediction to some real value. The **integrated risk** of a decision rule is then
$$
r(\delta) = E_{X,Y} \left[ \ell(Y, \delta(Y))\right].
$$

**Lemma**: The best possible rule is 
$$
\delta^*(x) = \argmin_a E_{Y \mid X = x}[\ell(Y, a)]
$$
and we call this the **Bayes decision rule** and the corresponding risk the **Bayes risk**. Everything pretty clearly extends to the case of continuous $Y$.

### Generative/Discriminative

In general, we need to estimate $\Pr(Y \mid X)$.

1. The **generative** approach is to model $P(Y)$ and $P(X \mid Y)$ and then do Bayesian inference.
2. The **discriminative** approach is to model $P(Y \mid X)$ directly.

#### Naive Bayes

From before, if we have $X$ high-dimensional, then we need some sort of structural assumption to estimate $P(X = x \mid Y)$ efficiently, since otherwise we would have to estimate the full conditional joint distribution, which can be quite bad. 

_Def_: Let us make the assumption of conditional independence, i.e.
$$
P(X_1 = x_1, \dots, X_p = x_p \mid Y = y) = \prod_{i=1}^p P(X_i = x_i \mid Y_i).
$$
We call this approach **naive** Bayes.

#### Logistic Regression

Here, the assumption is that 
$$
\log \left(\frac{\Pr(Y = 1 \mid X)}{\Pr(Y = 0 \mid X)} \right) = \beta^\top X 
$$
for some coefficient vector $\beta$.


#### KNNs

Here, the assumption is that $\Pr(Y = 1 \mid X)$ is the proportion of the nearest neighbors in the training data which have $Y = 1$.

