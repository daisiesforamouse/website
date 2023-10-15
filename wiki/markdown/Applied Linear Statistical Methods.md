---
title: 'Applied Linear Statistical Methods'
subtitle: 'UChicago STAT 34300, Autumn 2023'
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

Whenever it is unclear, vectors are column vectors.

## Simple Linear Regression

-------

[Wikipedia](https://en.wikipedia.org/wiki/Simple_linear_regression)

The setup of the simplest model is as follows.

Take $x_i \in \mathbb R$ as predictors (alternatively, independent variables, features, etc), and $y_i \in \mathbb R$ as responses (alternatively, outcomes, etc.). Then, we hopes that the response is approximately linear in the predictors, e.g.
$$
  y_i = \beta_0 + \beta_1 x_i + \epsilon_i
$$
which we split (semantically) into a systemic component and some error.

Our goal is then to minimize the residual sum of squares, e.g. 
$$
  \RSS(\mathbf \beta) = n^{-1}\sum_{i=1}^n \epsilon_i^2 = n^{-1}\sum_{i=1}^n(y_i - \beta_0 - \beta_1 x_i)^2.
$$

Do this however you want; it doesn't matter. You will arrive at
$$
\begin{align*}
  \hat \beta_1 &= \frac{\sum_{i=1}^n(x_i - \bar x)(y_i - \bar y)}{\sum_{i=1}^n(x_i - \bar x)^2} \\
  \hat \beta_0 &= \bar y  - \hat \beta_1 \bar x.
\end{align*}
$$

Now, none of the above has any particular randomness as defined. However, we can now impose (some of) the following assumptions:

- We have i.i.d. $(x_1, y_1), \dots, (x_n, y_n)$.
- We have independent $y_i \sim P_{y_i \mid x = x_i}$ given fixed $x_i$.
- **Linearity**: $E[y_i \mid x_i] = \beta_0 + \beta_1 x_i$, and thus our residuals obey $E[\epsilon_i \mid x_i] = 0$.
- **Homoskedasticity**: $\Var(y_i \mid x_i) = \sigma^2 > 0 \iff \Var(\epsilon_i \mid x_i) = \sigma^2 > 0$.

Now for example, if we assume linearity, we can compute that $\hat \beta_1 = \beta_1 + \widetilde \epsilon$ with
$$
  \Var(\widetilde \epsilon) = \frac{\sum(x_i - \bar x)^2 \Var(\epsilon_i)}{\left(\sum(x_i - \bar x)^2 \right)^2} = \frac{\sigma^2}{\sum(x_i - \bar x)^2} = \sigma_{\bar x}^2
$$
where the last equality only holds under homoskedasticity, and we call the last quantity the standard error.

In this case, we can compute that by the CLT the asymptotic distribution of $\hat \beta_{1}$ is $N(\beta_1, \sigma_{\bar x}^2)$. In particular, we can form a confidence interval by doing the usual stuff.

Testing proceeds the same way: set some null $H_0: \beta_1 = c$ and some alternative $H_A: \beta_1 \neq c$, and set 
$$
  \delta = \begin{cases}
    1 & \text{if } H_0 \text{ is rejected} \\
    0 & \text{otherwise} \\
  \end{cases}.
$$

Now take a test statistic $T \sim F_0$ under the null, with $p = 1 - F_0(T)$; this immediately yields that under the null, $p$ is uniformly distributed under the null (as long as $F_0$ is continuous), basically by definition.
Then, we set some threshold for $\alpha$, the probability of a Type-I error, and set $\delta = 1 \iff p \leq \alpha$.

Even under misspecification, we can see that if our samples are i.i.d. distributed $x_1, \dots, x_n \sim X$ and $y_1, \dots, y_n \sim Y$,
$$
  \hat \beta_1 = \frac{n^{-1} \sum_{i=1}^n (x_i - \bar x)(y_i - \bar y)}{n^{-1}\sum_{i=1}^n (x_i - \bar x)^2} \overset{D}{\to} \frac{\Cov(X, Y)}{\Var(X)} = \rho_{XY} \cdot \frac{\sigma_X}{\sigma_Y}
$$
where $\rho_{XY}$ is the correlation and $\sigma_X, \sigma_Y$ are standard deviations.

## Multiple Linear Regression

-------

[Wikipedia](https://en.wikipedia.org/wiki/Linear_regression)

Now we take $x_i \in \R^p$ and $y_i \in \R$, and set our model to be 
$$
  y_i = \sum_{j=1}^p\beta_j (x_i)_j + \epsilon_i.
$$
There are many different features that one could use.

- Taking a feature to be constantly $1$ is equivalent to adding a constant in our regression.
- One can take features that are functions of predictors (e.g. polynomials).
- You can take $\max \{ x_i - c, 0 \}$ to set a cutoff for a feature.

Notationally, set $y = \bmat{y_1 & \dots & y_n}^T \in \R^n$, $X = \bmat{x_1 & \dots & x_n}^T \in \R^{n \times p}$ where each $x_i$ is considered as column vector, $\beta = \bmat{ \beta_1 & \dots & \beta_n }$ and $\epsilon = \bmat{ \epsilon_1 & \dots & \epsilon_n }^T$ so that our model reduces to
$$
  y = X\beta + \epsilon.
$$

Then, the OLS estimator is given by
$$
  \hat \beta = \argmin\{ \RSS(\beta) = \| y - X \beta \|^2\}
$$
which we can compute directly taking a derivative (so long as $n \geq p$, e.g. we have more observations than features):
$$
  \nabla_\beta \RSS(\beta) = -2X^T(Y - X\beta) = 0 \implies X^TX\beta = X^Ty
$$
and if $X$ is invertible, we get
$$
  \hat \beta = (X^TX)^{-1}X^Ty.
$$

Why might $X$ not be invertible? Some examples are

- duplicated features,
- unit conversions / identical measurements,
- batch effects (e.g. one feature is if a patient was given an medication $A$, and another is if a patient was given it by technician $A'$, but it was exactly $A'$ who handed out $A$).

### OLS as a Projection

From a geometric perspective, the OLS predictions are just the projections of $y$ into the image of $X$ as a linear map $\R^p \to \R^n$.

_Def_: $X_{\cdot, 1}, \dots, X_{\cdot p}$ are **orthogonal** if $\left\langle X_{\cdot j}, X_{\cdot k}\right\rangle = 0$ for all $j, k$, and **orthonormal** if they are all of length 1.

**Theorem**: Set $V \subset \R^n$ a linear subspace and $y \in \R^n$; there exists a unique vector
$$
  \proj_V(y) = \argmin_{v \in V} \| y - v \|.
$$
Furthermore, $y - \proj_V(y) \in V^{\perp}$, the perpendicular space to $V$ (or equivalently, $\left\langle y - \proj_V(y), v \right\rangle = 0$ for all $v \in V$).

There are a few more facts:

- $y \mapsto \proj_V(y)$ is a linear operator, and the corresponding matrix is called the projection matrix $P_V \in \R^{n \times n} = P_V$;
- $P_V$ is idempotent, e.g. $P_V^2 = P_V$;
- $P_V^T = P_V$;
- $P_{V^\perp} = I - P_V$, which gives a decomposition $y = P_Vy + P_{V^\perp} y$.
- $\rank(P_V) = \dim(V)$, and has a eigenvalue decomposition with eigenvalues $\{ 0, 1 \}$ in the obvious way.

**Corollary**: Set $P_X = P_{\Im(X)}, \hat y = P_X y = X \hat \beta$. This immediately yields that the "fit" $\hat y$ and the residuals $\hat \epsilon$ are unique, and
$$
  \|y\|^2 = \| \hat y \|^2 + \| \hat \epsilon \|^2
$$
since $\hat \epsilon = (I - P_{X})y$.

_Example_: 

- If $X_{k 1} = 1$ for all $k$ (that is, we include a constant in our regression), we immediately get that $\sum_{i=1}^m \epsilon_i = 0$. 
- If $X$ is full rank, then $P_X = X(X^TX)^{-1}X^T$.

### Reducing to Simple Linear Regression
  
If one column $X_{\cdot j}$ is orthogonal to every other feature, then this reduces to simple linear regressions: 
$$
  \hat \beta_j = \frac{\left\langle X_{\cdot j}, y\right\rangle}{\| X_{\cdot j} \|^2}.
$$

Using this, set 

- $Z_{\cdot 1} = X_{ \cdot 1}$,
- $Z_{\cdot n}$ to the residuals of $X_{\cdot n} \sim Z_{\cdot n-1}$, e.g.
$$
  X_{\cdot n} - \sum_{i=1}^{n-1}\frac{\left\langle Z_{\cdot i}, X_{\cdot i}\right\rangle}{\| Z_{\cdot i}\|^2}Z_{\cdot i}.
$$

Then, we see that
$$
  \hat \beta_p = \frac{\left\langle Z_{\cdot p}, y \right\rangle}{\| Z_{\cdot p}^2\|},
$$
so in general each coefficient depends on other variables, e.g.
$$
  \hat \beta_j = \frac{\left\langle \text{residuals of } X_{\cdot j} \sim \sum_{k \neq j}X_{\cdot k}, y\right\rangle}{ \| \text{residuals of } X_{\cdot j} \sim \sum_{k \neq j}X_{\cdot k} \|^2}.
$$

In particular, the above (which is clearly just doing Gram-Schmidt) gives us a $QR$-decomposition of $X$, where $Q = ZD$ is the orthonormalization of $Z$ and $R$ is the upper triangular matrix corresponding to Gram-Schmidt.

This is (for the most part) what software packages do in computing linear models: you set $y = Q\gamma + \epsilon$, and take $\hat \gamma = Q^Ty$; but since the fit is unique, we have that $Q \hat \gamma = X \hat \beta = QR\hat \beta \implies R\hat \beta = \hat \gamma$, which is numerically tractible.

### OLS and SVD ###

_Def_: [Link](https://en.wikipedia.org/wiki/Singular_value_decomposition) The **SVD decomposition** of a real (resp. complex) matrix $X \in \R^{n \times p}$ (resp. $\C^{n \times p}$) is
$$
  X = U \Sigma V^*
$$
where $U, V$ are orthogonal (resp. unitary) and $\Sigma$ is diagonal. Alternatively, we may sometimes use the "skinny" SVD, where $U, \Sigma$ are just $n \times p$ and $p \times p$ matrices, essentially by dropping the zero rows of $\Sigma$.

SVD is related to eigenvalue decomposition by noting that
$$
  X^*X = V\Sigma V^* U \Sigma V^* = V(\Sigma^* \Sigma)V^*
$$
so that columns of $V$ are eigenvectors of $X^*X$, and similarly columns of $U$ are eigenvectors of $XX^*$.

Then, applying to least squares, we have that
$$
\begin{align*}
  \| y - X \beta \|^2 &= \| y - U \Sigma V^T \beta \|^2  \\
  &= \| U^Ty - \Sigma V^T B \| \\
  &= \| U^Ty - \Sigma \beta^* \|^2 \\
  &= \sum_{j=1}^p ((U^Ty)_j - \sigma_j \beta_j^*)^2 \\
  &= \sum_{j=1}^k (U^Ty)_j - \sigma_j \beta_j^*)^2 + \sum_{j=k+1}^p (U^Ty)_j)^2 \\
\end{align*}
$$

So minimizing with respect to $\beta_j^*$, we pick
$$
  \hat \beta_j^* = \begin{cases}
    \frac{(U^Ty)_j}{\sigma_j} &  j = 1, \dots, k \\
    \text{anything} & j = k + 1, \dots, p
  \end{cases}
$$
and in particular if we take the arbitrary values to be 0 we get the "minimal norm" solution and since $\beta^* = V^T \beta$, the "ridgeless" solution is $\hat \beta = V\beta^*$, and is
$$
  \argmin \{ \| \beta \| \mid \beta \in \R^p, X^TX\beta = X^Ty \}.
$$

## Distributions

-------

We covered a bunch of stuff that you can probably find on Wikipedia.

- [Multivariate Gaussian Distribution](https://en.wikipedia.org/wiki/Multivariate_normal_distribution)
- [Chi-Square Distribution](https://en.wikipedia.org/wiki/Chi-squared_distribution)
- [Student-t Distribution](https://en.wikipedia.org/wiki/Student%27s_t-distribution)
- [F Distribution](https://en.wikipedia.org/wiki/F-distribution)

## Inference in the Homoskedastic Normal Linear Model

-------

For this section, set $X \in \R^{n \times p}$ design matrix, which is taken to be deterministic and of full rank. Then, set $y \sim N(X \beta, \sigma^2 I)$ for some $\sigma > 0, \beta \in \R^p$; that is, set
$$
  Y = X\beta + \epsilon, \epsilon \sim N(0, \sigma^2 I).
$$
We have already seen
$$
  \hat \beta = (X^TX)^{-1}X^TY
$$
and so $\hat \beta$ is normal with mean $\beta$ and variance
$$
  \Var(\hat \beta) = (X^TX)^{-1}X^T(\sigma^2 I)X(X^TX)^{-1} = \sigma^2 (X^TX)^{-1}.
$$

But now we need to understand $\sigma$: we may split $Y$ into $\hat Y = P_XY$ and $\hat \epsilon = (I - P_X)Y$, where $P_X$ is the orthogonal projection matrix onto $X$; properties of the multivariate Gaussian imply that these two are independent, and have distributions
$$
  \begin{align*}
    \hat Y &\sim N(X \beta, \sigma^2 P_X) \\
    \hat \epsilon &\sim N(0, \sigma^2(I - P_X))
  \end{align*}
$$

which we can now use $\hat \epsilon$ to learn about $\sigma^2$ and $\hat Y$ to learn about $\beta$.

Specifically, we use the residual sum of squares $RSS = \| \hat \epsilon \|^2 \sim \sigma^2 \chi^2_{n-p}$; thus this has expectation $\sigma^2(n - p)$, so we can set
$$
  \hat \sigma^2 = \frac{\| \epsilon^2 \|}{n - p} \sim \frac{\sigma^2 \chi^2_{n - p}}{n - p}.
$$
Moreover, since $X^T = X^TP_X$, we already know that $\hat \beta$ only depends on $\hat y$, and so is independent of $\hat \sigma$.

### Applications

We can now do inference for $\beta_j$; in particular $\hat \beta_j \sim N(\beta_j, \sigma^2 (X^TX)_{jj}^{-1})$. Then we know that if we take
$$
\begin{align*}
  A &= \frac{\hat \beta_j -  \beta_j}{\sigma \sqrt{(X^TX)^{-1}_jj}} \sim N(0, 1) \\
  B &= \frac{\hat \sigma^2 }{\sigma^2} \sim \frac{\chi^2_{n-p}}{n - p}
\end{align*}
$$
so $A / \sqrt{B} \sim t_{n - p}$; that is,
$$
  \frac{\hat \beta_j -  \beta_j}{\hat \sigma \sqrt{(X^TX)^{-1}_jj}} \sim t_{n - p}
$$
so we can get a concrete handle on the sampling distributions from only our observations.

