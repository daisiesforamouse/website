---
title: 'Diffusion Models'
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
\DeclareMathOperator{\tr}{tr}

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


## Basics of Diffusion Models

------------------------

We are considering the problem of sampling from some unknown distribution $P_0$, potentially of very high dimension (such as images).

A diffusion model is composed of a forwards and a backwards process. In the forwards case, a clean sample from the data distribution is progressively contaminated by Gaussian noise, while in the backwards case we sequentially remove noise in hopes of restoring $P_{\text{data}}$. That is, we run two stochastic processes:
the forwards process
$$
    \text{Data} \to \text{Data + Noise} \to \text{Noise}
$$
and the backwards process
$$
    \text{Noise} \to \text{New Data + Noise} \to \text{New Data}
$$
where we learn the backwards process during the forwards pass.

### Forwards Processes

_Def_: For some $X_0 \sim P_{0}$, the **forwards process** is given by an Ornstein-Ulhenbeck process
\begin{equation*}
    dX_t = -\frac{1}{2} g(t)X_t dt + \sqrt{g(t)} dW_t
\end{equation*}
where $W_t$ is a standard Brownian motion and $g(t)$ is a positive nondecreasing weighting function called the **variance schedule**.

While the marginal distribution of $X_t$, (written as $P_t$) is generally intractible, we do know the conditional distribution $X_t \mid X_0$.

_Lemma_: We have that
$$
    X_t \mid X_0 \sim N \left(\alpha(t) X_0, h(t) I\right)
$$
where 
$$
    \alpha(t) = \exp \left(-\frac{1}{2}\int_0^t g(s)ds \right) \text{ and } h(t) = 1 - \alpha(t)^2.
$$

_Proof_: Consider Ito's lemma with $f(t, x) = \frac{X_t}{\alpha(t)}$ applied to $X_t$, which gives that
$$
    d\left(\frac{X_t}{\alpha(t)}\right) = \sqrt{g(t)} \alpha(t)dW_s
$$
which gives that
$$
    X_t = \alpha(t) \left(X_0 + \int_0^t \sqrt{g(s)} \alpha(s)dW_s \right).
$$
The above has the required distribution. $\square$

This means that if $\alpha(t) \to 0$ as $t \to \infty$ (under regularity conditions such that the above stochastic calculus goes through - for example, Lipschitz and bounded is sufficient), then the limiting distribution is $N(0, I)$.

### Backwards Processes

_Def_: Take a general diffusion
$$
    dX_t = \mu(t, X_t)dt + \sigma(t, X_t)dW_t.
$$
Then, the corresponding **backwards SDE** is
$$
    d\widetilde X_t = \widetilde \mu(t, X_t)dt + \widetilde \sigma(t, X_t) d\widetilde W_t
$$
for suitable $\widetilde \mu, \widetilde \sigma$ and Brownian $\widetilde W_t$ such that $\widetilde X_{t} = X_{T - t}$ holds in distribution for some later time $T$, start time $t_0$, and all $t_0 \leq t \leq T$.

**Theorem (Anderson 1982)**: For a general diffusion as above, let $\mu, \sigma$ be such that they guarantee the existence of the probability density of $X_t$, denoted $p_t(x)$, as a smooth and unique solution to the associated Kolmogorov forwards/backwards equations for $t_0 \leq t \leq T$, where $t_0, T$ are starting and end times. Then, we have that the backwards equation governed by 
$$
    d\widetilde X_t = \widetilde \mu(T - t, \widetilde X_t)dt + \sigma(T - t, \widetilde X_t) d\widetilde W_t
$$
where
$$
    \widetilde \mu^i(t, x) = -\mu^i(t, x) + \frac{1}{p_{t}(x)} \sum_{j, k} \frac{\partial}{\partial x^j} \left[p_{t}(x)\sigma^{ik}(t, x)\sigma^{jk}(t, x)\right].
$$
Here, superscripts denote the relevant indices for vector/matrix valued functions (e.g. $\mu^{1}$ is the first coordinate of $\mu$).

_Proof_: See Anderson 1982 for the proof and sufficient conditions on $\mu, \sigma$; note that he uses a slightly different convention so that there is a difference in sign. Essentially, consider the joint density and use Kolmogorov forwards/backwards equations. $\square$

Note that Anderson actually gives $d\widetilde W_t$ explicitly from $dW_t$, which we do not use here, as well more sophisticated regularity conditions on $\mu, \sigma$.

**Corollary**: In the case of the OU diffusion above, we get that the reverse time model is simply 
$$
    d\widetilde X_t = \left(\frac{1}{2} g(T-t) X_t + g(T-t) \nabla \log(p_{T-t}(X_t)) \right)dt + \sqrt{g(T-t)}d\widetilde W_t.
$$

_Def_: We call $\nabla \log(p_t(x))$ the **score** of $p_t$. Note that this is different from the usual meaning of "score" in statistics, where the derivative is with respect to the parameters and not the data.

### Conditional Models

These are very similar in spirit to the unconditional models above; here, however, instead of trying to sample from the unknown $P_0$, we wish to general samples from a conditional data distribution $P_0(\cdot \mid y)$. The conditional forward process is exactly the same:
$$
    dX_t^y = -\frac{1}{2} g(t) X^y_t dt + \sqrt{g(t)}dW_t
$$
with $X_0^y \sim P_0(\cdot \mid y)$. However, the backwards process differs insofaras the score is replaced by the conditional score, i.e.
$$
    d\widetilde X_t^y = \left(\frac{1}{2} g(T-t) X_t^y + g(T-t) \nabla \log(p_{T-t}(X_t^y \mid y)) \right)dt + \sqrt{g(T-t)}d\widetilde W_t.
$$

## The Score

------------------------

Of course, we don't have access to the score $\nabla \log(p_t(x))$, so we need to estimate it from the forwards process.

### Score Matching

The most prominent method for score estimation is score matching, introduced in Hyvarinen 2005 for a different problem, i.e. estimating parameters for nonnormalized statistical models where the partition function is intractible.

The procedure is as follows: for some distribution $P_0$ with corresponding density $p_0(x)$, we wish to solve the minimization problem
$$
    \arg\min_{s \in \mathcal S} E \left[ \| s(x) - \nabla \log(p_{0}(x)) \|_{2}^2 \right]
$$
for some suitable class $\mathcal S$. Of course, this objective is intractible since $p_0$ is generally unknown. However, it turns out to be equivalent to a tractible loss.

**Theorem** (Hyvarinen 2005): Suppose that everything in sight is differentiable and has finite second moment. Then, we have the equivalence 
$$
    \frac{1}{2} E_{x \sim P_0} \left[ \| s(x) - \nabla \log(p_{0}(x)) \|_{2}^2 \right] = E_{x \sim P_0} \left[ \tr(\nabla s(x)) + \frac{1}{2} \|s(x)\|_2^2 \right] + C
$$
for some constant $C$.

However, the computation of $\tr(\nabla s(x))$ does not lend itself to high dimensional problems nor deep networks; thus we have to turn to a slightly different approach.
