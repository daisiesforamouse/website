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

We can also form a confidence set for $\beta$: since $\hat y \sim N(X \beta, \sigma^2 P_X)$,
$$
  \|X(\beta - \hat\beta)\|_2^2 \sim \sigma^2\chi_p^2
$$
and so
$$
  \frac{p^{-1}\|X(\beta - \hat\beta)\|_2^2}{\hat \sigma^2} \sim F_{p, n-p}.
$$
Then, we can take the confidence set to be the set of $\beta$ that falls below some quantile of the $F$ distribution.

For predictions, if we have some extra data point $x_{n+1}$,
$$
  x_{n+1}^T \hat \beta \sim N(x_{n+1}^T\beta, \sigma^2 x_{n+1}^T(X^TX)^{-1}x_{n+1})
$$
so we have a pivot
$$
  \frac{x_{n+1}^T\beta}{\sigma \sqrt{x_{n+1}^T(X^TX)^{-1}x_{n+1}}} \sim t_{n-p}.
$$

Then, since
$$
  y_{n+1} - x_{n+1}^T \hat \beta = x_{n+1}^T(\beta - \hat \beta) + \epsilon_{n+1} \sim N(0, \sigma^2(x_{n+1}^T(X^TX)^{-1}x_{n+1} + 1)
$$
we can form a confidence interval.

### Nested Model Comparisons

For example, consider the subspace $U$ of $X$ spanned by $(1, 1, \dots, 1)$. Then,
$$
  \|(I - P_U)y^2\|^2 = \|(P_X - P_U)y\| + \|(I - P_X)y\|^2
$$
so
$$
  \|y - \bar y\|^2 = \|\hat y - \bar y\|^2 + \|\hat y - y\|^2
$$
and we label the above equation as $\operatorname{SST} = \operatorname{SSReg} + \operatorname{RSS}$ (the total sum of squares, the regression sum of squares, and the residual sum of squares); further
$$
  R^2 = \frac{\operatorname{SSReg}}{\operatorname{SST}} = 1 - \frac{\operatorname{RSS}}{\operatorname{SSReg}}.
$$
Sometimes one sees adjusted $R^2$, e.g.
$$
  R^2 = \frac{\operatorname{SSReg}}{\operatorname{SST}} = 1 - \frac{\operatorname{RSS} / (n-p)}{\operatorname{SSReg}/(n-1)}.
$$
which includes an adjustment for degrees of freedom. Note that the above assumes that we fit a model with an intercept term; if there is no intercept term, $R^2$ becomes mostly useless (though in reality it's a hard statistic to interpret anyway).

Now suppose that we have two models,
$$
  \text{Model A: } y = Z\beta + \epsilon
$$
and
$$
  \text{Model B: } y = X\gamma + \epsilon.
$$
where the columns of $Z \in \R^{n \times k}$ span a subspace of the column space of $X$. Then,
$$
  \|(I - P_Z)y\|^2 = \|(P_X - P_Z)y\|^2 + \|(I - P_X)y \|^2  
$$
which gives $\RSS_A \geq \RSS_B$, so taking a larger model always reduces $R^2$. But if model A was the "true" model, then
$$
  \|(P_X - P_Z)y\|^2 \sim \sigma^2 \chi^2_{p-k}
$$
and
$$
  \frac{\|(P_X - P_Z)y\|^2/(p-k)}{\hat \sigma^2} = \frac{(\RSS_A - \RSS_B)/(p-k)}{\hat \sigma^2} \sim F_{p-k, n-p}.
$$

### Dropping the Normality Assumption

In the above, we have taken $y \sim N(X\beta, \sigma^2 I_n)$. But in general, we might only know the first two moments of $y$ without knowing a distribution. However, if we have some weaker assumptions, e.g. $y$ has reasonable tails and each $x_i$ is reasonably close to the average, then the central limit theorem still yields that $\hat \beta$ is asymptotically $N(\beta, (X^TX)^{-1} \sigma^2)$.

However, in practice we still compare to the $t_{n-p}$ distribution, since it gives more conservative results, even though all we have is asymptotic normality.

## Diagnostics

Let $Y \sim N(X \beta,\sigma^2 I_n)$, and in this case consider $X$ a random matrix. 

If the model is correct, then $E[\hat \epsilon_i \mid \hat y_i] = 0$, so if we look and this is not the case, the linearity assumption is usually wrong.

Moreover, if we look at the distribution of the residuals, we see that they are not necessarily of identical variance, since they have variance $\sigma^2(I - P_X)$, so if we take
$$
    \epsilon^{m}_i = \frac{\epsilon^m_i}{\hat \sigma \sqrt{1 - h_{i}}}
$$
where $h_i$ is the $i$-th diagonal entry of $I - P_X$.

TODO (missed lecture)

## Inference in the Presence of Heteroskedasticity

Let $X \in \R^{n \times p}$ be fixed, and $Y \sim (X\beta, \Sigma)$. Then, we still have that $\hat \beta^{OLS} = (X^TX)^{-1}X^TY$, which still is $\beta$ in expectation, but the variance is not as nice; in particular
$$
    \Var(\hat \beta^{OLS}) = (X^TX)^{-1}X^T\Sigma X(X^TX)^{-1}.
$$
Thus, if you do inference assuming that homoskedasticity holds, we can get very misleading results.

So there are a couple of things that we could do:
- come up with a new estimator,
- or still work with the OLS estimate and come up with new inference conditions.

For a new estimator, suppose that $\Sigma = \sigma^2 \Gamma > 0$ (one might know this from weaker assumptions). Then,
$$
    \Gamma^{-1/2}Y = \Gamma^{-1/2}X\beta + \Gamma^{-1/2}\epsilon
$$
is a homoskedastic model. Now we can do OLS with $\tilde Y = \Gamma^{-1/2}Y$ and $\tilde X = \Gamma^{-1/2}X$; this is called the GLS, or generalized least squares. In fact,
$$
    \hat \beta^{GLS} = (X^T\Gamma^{-1}X)^{-1}X^T\Gamma^{-1}Y
$$
which has first moment $\beta$ and variance $\sigma^2X^T\Gamma^{-1}X^{-1}$. In fact, $\Var(c^T\hat \beta^{GLS}) \leq \Var(c^T\hat\beta^{OLS})$ for any $c \in \R^p$ (this is just Gauss-Markov). In the case that $\Gamma$ is diagonal, this is just the weighted least squares.

If we choose to work with the OLS estimator, then we need to be very careful when computing. Surprisingly, if $\Sigma$ is diagonal, we can still estimate the variance; in particular we need to know 
$$
    X^T\Sigma X = \sum_{i=1}^n x_ix_i^T \sigma_i^2
$$
and we can take a bad temporary estimate, $\hat \sigma_i^2 = \hat \epsilon_i^2$. However, this gives
$$
    \frac{\widehat{X^T \Sigma X}}{X^T \Sigma X} \to I
$$
in probability. So the estimated variance is
$$
    \widehat{\Var(\hat \beta^{OLS})} = (X^TX)^{-1} \left( \sum_{i=1}^n x_ix_i^T\hat \epsilon_i^2 \right)(X^TX)^{-1}
$$
which gives what is called the (heteroskedasticity) robust standard error.

## Assumption-Lean Regression

If we let $(x_i, y_i) \sim P$ be i.i.d. draws from $\R^{p} \times \R$; we do not require that $E[y_i \mid x_i] = \beta^T x_i$, nor that $\Var(y_i \mid x_i) = \sigma^2 I$.

Before, $\hat \beta^{OLS} \in \argmin_{\beta \in \R^p} \{ n^{-1}\sum_{i=1}^n(y_i - \beta^T x_i)^2 \}$, and now we consider
$$
    \beta^*(P) = \argmin \{ E[(y_i - \beta^T x_i] \}.
$$
If $E[x_ix_i^T] > 0$, then
$$
    \beta^*(P)  = E[x_ix_i^T]^{-1}E[x_i^Ty_i].
$$

We can do a computation to see that
$$
\begin{align*}
    \hat \beta_{OLS} - \beta^*(P) &= (X^TX)^{-1}X^T(Y- X\beta^*(P)) \\
    &= \left( \frac{X^TX}{n}\right)^{-1} \cdot \frac{1}{n} \sum_{i=1}^n x_i(y_i - x_i^T\beta^*(P)).
\end{align*}
$$
But $\frac{X^TX}{n} \to E[x_ix_i^T]$ by the law of large numbers, and
$$
    \frac{1}{n}\sum_{i=1}^n x_i(y_i - x_i^T\beta^*(P) \to E[x_i(y_i - x_i^T\beta^*{P})] = 0
$$
so in probability, $\hat \beta_{OLS} \to \beta^*(P)$. Furthermore, the central limit theorem gives that
$$
    \frac{1}{\sqrt{n}} \sum_{i=1}^n x_i(y_i - x_i^T \beta^*(P)) \to N(0, E[x_ix_i^T(y_i - x_i^T\beta^*(P))^2])
$$
in distribution, so we know that
$$
    \sqrt{n}(\hat \beta_{OLS} - \beta^*(P)) \to N(0, E[x_ix_i^T] E[x_ix_i^T(y_i - x_i^T\beta^*(P))^2]E[x_ix_i^T]^{-1}).
$$

The associated finite estimator is called the Eicker-Huber-White standard error, and it looks like the heteroskedastic robust standard error; in fact this shows that it is also non-linearity robust.

Now, the conditional expectation (put here as $\mu^*(X)$) is the minimizer of the $L2$ distance among all functions measurable w.r.t. $X$; then if you define $\ell^*(X) = X^T\beta^*(p)$, then $\ell^*$ is the minimize of the $L2$ distance among all linear functions. Then, writing
$$
    y_i = \ell^*(x_i) + \mu^*(x_i) - \ell^*(x_i) + y_i - \mu^*(x_i) 
$$
we see that the first difference is some non-linearity error, and the last difference is some unexplainable error.

## Bootstrapping

Bootstrapping is about inference by simulation.

Set $y_1, \dots, y_n \sim P$ as i.i.d. draws on $\R$, with mean $\mu$ and variance $\sigma^2$. Set
$$
    \bar y = \frac{1}{y}\sum y_i, \hat \sigma^2 = \frac{1}{n} \sum (y_i - \bar y)^2.
$$

Then for specific distribution we can figure things out analytically, but we may or may not know $P$ in reality. In this case, what one could do is draw
$$
    y_1^{(1)}, \dots, y_n^{(1)}
$$
to
$$
    y_1^{(B)}, \dots, y_n^{(B)}
$$
from $P$ so we can figure out the distribution of $\bar y$; but since we know $P$ approximately by the original distribution of $y_1, \dots, y_n$, we can just simulate by estimating $P$ by $\hat P$, e.g. taking repeated draws by sampling with replacement.

This is equivalent to the distribution
$$
    \hat P = \frac{1}{n} \sum \delta(y_i)
$$
which gives rise to the so-called empirical distribution function. In fact this might seem crude, but over large sample sizes in low dimensions, this can be quite good.

Then, this lets you get the above statistics in obvious ways. You can also get confidence intervals by, say, using asymptotic normality, or by using the $1 - \alpha/2$ and $\alpha/2$ quantiles, or by using the "bootstrap $t$". This last one is slightly more subtle: since we approximate
$$
    \frac{\bar y - \mu}{\sigma^2 / \sqrt{n}} \text{ by } \frac{\bar y^{(B)} - \bar y}{\hat \sigma^{(B)} / \sqrt{n}}
$$
we can set the interval to be $[\bar y - \frac{\hat \sigma}{\sqrt{n}} u,\bar y + \frac{\hat \sigma^2}{\sqrt{n}} l]$ where $u, v$ are the $1 - \alpha/2, \alpha/2$ quantiles of the $\frac{\bar y^{(B)} - \bar y}{\hat \sigma^{(B)} / \sqrt{n}}$ distribution.

### Bootstrapping for Linear Regression

Take again the general setup of $(x_i, y_i) \sim P$ i.i.d. Then, the pairs bootstrap is done by taking the empirical distribution of $(x_i, y_i)$

If we need $x_i$ to be fixed, then we can do the residual bootstrap, which just means fitting least squares and then bootstrapping over the residuals. But note that this assumes both linearity and homoskedasticity.

The wild bootstrap allows us to drop homoskedasticity; to perform it you fix $x_i^{(b)} = x_i$, and $y_i^{(b)} = \hat \beta^T x_i + \hat \epsilon_i \cdot \eta_i^{(\beta)}$ where $\eta_i^{(b)}$ are i.i.d. (often Rademacher); alternatively,
$$
    \eta_i^{(b)} = \begin{cases}
        \phi & \text{ w.p. } \frac{\phi^{-1}}{\phi + \phi^{-1}} \\
        \phi^{-1} & \text{ w.p. } \frac{\phi}{\phi + \phi^{-1}} \\
    \end{cases}
$$
where $\phi$ is the golden ratio fixes the third moment as well.

If we want to dodge the issue in the pairs bootstrap (where the new design matrix is not of full rank) we can do the Bayesian bootstrap.

Suppose that we had dependence in some structured way, such as time dependence; then we can discretize the observations corresponding to some blocks that isolate the dependence (for example, so blocks corresponding to time intervals). This is called the "box bootstrap".

### Permutation

If we permute the responses, e.g. create some new dataset $(x_{1}, y_{\sigma(1)}), \dots, (x_n, y_{\sigma(n)})$ (labeled as $x_i^{(b)}, y_i^{(b)}$), we can test the null hypothesis that $\beta_1 = 0$ in $y_i = \beta_0 + \beta_1 x_i + \epsilon_i$ by computing the $p$-value
$$
    \frac{1 + \sum_{b=1}^B 1_{|\hat \beta_1^{(b)}| \geq |\hat \beta_1|}}{1 + B}
$$
which gives us something valid in finite samples.

If you want to test one specific predictor, you can first bootstrap the residuals with the model with that predictor removed, and then compute as above.

## Cross Validation

As before, assume we have some dataset $D = \{(x_1, y_1), \dots, (x_n, y_n)\}$, to which we fit some $\hat \mu$ that takes new predictors to an estimated response.

Then if we have some test dataset $(x_{n + 1}, y_{n + 1}), \dots, (x_{n + k}, y_{n + k})$, then we can compute the external error
$$
    \frac{1}{k} \sum_{i=1}^k (y_{n+i} - \hat \mu(x_{n+i})).
$$

If we don't have external testing data, then we can cross-validate, e.g. by creating some hold-out blocks $I_1, \dots, I_k$ we don't train on, and taking the cross-validated error
$$
    \frac{1}{n} \sum_{i=1}^k \sum_{j \in I_l} (y_i - \hat \mu^{(i)}(x_j))
$$
where each $\mu^{(i)}$ is fit on everything but $I_i$. This is called $k$-fold CV; if we have $n$ folds, this is called leave one out CV.

This is very cheap in OLS, since
$$
    \hat y_i = h_i y_i + (1 - h_i)\hat y_{i,-i}
$$
where $\hat y_{i, -i}$ is the predicted value when not fit with $(x_i, y_i)$ and $h_i = (P_X)_{ii}$ is the leverage. So we can compute the leave out out error quickly. In fact, the total CV error is just
$$
    \frac{1}{n} \sum_{i=1}^n \frac{\hat \epsilon_i^2}{(1 - h_i)^2}.
$$

We can use this to select a model with good performance, and then we often refit on the whole dataset after selecting our model.

Alternatively, if we assume homoskedasticity, then we can redraw responses $y_i'\sim (\hat \mu(x_i), \sigma^2)$ and compute the error
$$
    \frac{1}{n}E \left[ \sum_{i=1}^n (y_i' - \hat \mu(x_i))^2\right].
$$
But 
$$
    E[(y_i' - \hat \mu(x_i))^2] = 2\sigma^2 +  E[(y_i - \hat y_i)^2] + 2E[(\hat \mu(x_i) - y_i)(y_i - \hat y_i)]
$$
but $E[(\hat \mu(x_i) - y_i)(y_i - \hat y_i)] = -2\sigma^2 -2\Cov(y_i, \hat y_i)$, so this is just
$$
    E[(y_i' - \hat y_i)^2] = E[(y_i - \hat y_i)^2] + 2 \Cov(y_i, \hat y_i).
$$
If we put $df(\hat y) = \sum_{i=1}^n \frac{\Cov(y_i, \hat y_i)}{\sigma^2}$ (this is called the degrees of freedom of our model), then the error above is just
$$
    \frac{1}{n}E \left[\sum_{i=1}^n (y_i - \hat y_i)^2\right] + \frac{2\sigma^2}{n}df(\hat y).
$$

To see that degrees of freedom is a reasonable name, note that if we were doing a OLS model, then $df(\hat y) = \rank(X)$ if $X$ is full column rank. This gives a way to choose between models with different amounts of predictors.

## Discovery Rates

Suppose you have a family of $n$ hypotheses, and you want the probability of a false discovery to be less than some predefined value $p$.

_Def_: The **Bonferroni correction** is when you require every single hypothesis to have error rate bounded by $p / n$.

Popularized by Benjamini and Hochberg, we define
$$
    FDR = E \left[ \frac{\text{number of false discoveries}}{\max\{\text{number of all discoveries}, 1\}} \right].
$$
Then, let $k^* = \max_{k = 11, \dots, n}\{p_{k} \leq \frac{\alpha k}{n}\}$ where $p_1 \leq \dots \leq p_n$. Then we reject $p_1, \dots, p_{k^*}$, where $\alpha$ is the desired false discovery rate.

**Theorem**: Up to non-pathological distributions, the FDR is bounded by $\alpha$; in pathological cases we still have an upper bound of $\alpha \log(n)$.

## Robust Regression

OLS is very sensitive to outliers, which is something that we sometimes want to fix. There are a few options.

- Least Absolute Deviation: minimize $|Y - X\beta|_1$ instead of $|Y - X\beta|_2$. Under Gaussian linear models, these things converge to the same thing and in fact OLS goes a little faster; under misspecification, this is not true. This is a convex optimization problem.
- Quantile Regression: Let 
$$
    q_\alpha(u) = \begin{cases}
        \alpha |u| & u > 0\\
        (1 - \alpha)|u| & \text{otherwise}
    \end{cases}
$$
and minimize $\sum_{i=1}^n q_\alpha(y_i - x_i^T\beta)$.
- Least trimmed squares.
- Huber regression.

We don't have closed forms for the estimators, and asymptotics are unwieldy for these options, so we do inference via bootstrap and get estimators via iterative optimization.


##  OLS Modifications

### Regression With Intercept

As usual, set $y_i \in \R$, $x_i \in \R^n$ with the model
$$
    y_i = \beta_0 + x_i^T\beta + \epsilon_o
$$
so that we can, by rescaling $\tilde y_i = y - \bar y, \tilde x_ = x_i - \bar x$, force the OLS line through the origin.

### Multivariate Analysis

Suppose we have $x_1, \dots, x_n \sim P$, and $\Sigma = \Var(x_i)$. Suppose that we wish to create new features $z_i = c^Tx_i$ for some unit vector $c$ so that we maximize the variance of $z_i$. Since $\Var(z_i) = c^T \Sigma c$, the maximal choice is the top eigenvector of $\Sigma$; going down the list of eigenvectors gives us the principal components.

In practice, $\Sigma$ is estimated by $X^TX$ for a design matrix $X$, so this reduces to taking the SVD. Principal component regression is just dropping the bottom components. Since PCA is sensitive to the scale of the predictors, we often normalize the centered predictors to have norm 1.

### Ridge/Lasso Regression

We can add a $L^2$ penalty by minimizing
$$
    \|Y  - X\beta\|_2^2 + \lambda \|\beta\|_2^2.
$$
Sometimes we drop the intercept coefficient in the penalty so that we do not force the model to pass through the origin at large $\lambda$; this is the same as centered first. An $L^1$ penalty is given by minimizing
$$
    \|Y  - X\beta\|_2^2 + \lambda \|\beta\|_1
$$
which gives us sparser models. Lastly, we might impose an $L^0$ penalty
$$
    \|Y  - X\beta\|_2^2 + \lambda 1_{\beta \neq 0}.
$$
This is nonconvex, but we can do a greedy search (which is not guaranteed to be optimal).

## Causal Inference

As usual, let $y_i, x_i$ be the response and predictors. We introduce the notation $y_i(0)$ for the response without treatment and $y_i(1)$ for the response with treatment, and the individual fit effect is defined as $ITE_i = y_i(1) - y_i(0)$. Now by assumption $y_i = y_i(w_i)$, and we cannot observe both $y_i(0)$ and $y_i(1)$, so of course $ITE_i$ is unknowable. But we might want to control the moments of $ITE_i$ if we let $(w_i, y_i(0), y_i(1)) \sim P$.

In a randomized control trial (RCT), $E[w_i] = p$ and is independent of everything else. Even without RCT
$$
    \widehat{TDM} \to E[y_i \mid w_i = 1] - E[y_i \mid w_i = 0]
$$
in probability. With RCT,
$$
    \widehat{TDM} \to E[y_i(1)] - E[y_i(0)].
$$

We can show two central limit theorems, where
$$
    \sqrt{n}(\widehat{TDM} - ITE) \to N(0, \sigma^2)
$$
for some variance depending on $P$.
