---
title: 'Brownian Motion and Stochastic Calculus'
subtitle: 'UChicago STAT 38510, Autumn 2023'
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

##  Brownian Motion 

-------

As a style preference, I'm going to drop all the arguments that are from the probability space.

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

###  Construction 

Pick a probability space $(\Omega, \mathcal F, P)$ that is rich enough to support a countable collection of independent standard normal variables. If you are particular, the unit interval with Lesbegue measure is sufficient here.

The strategy is as follows: we define $B_t$ for a countable dense set (in particular the dyadic rationals) of times using our precession of standard normals. then, we find some $t \mapsto B_t$ that agrees on the dense set and is uniformly continuous and then extend by continuity.

Set $D_n = \left\{ \frac{k}{2^n}, k = 0, 1, \dots, 2^n \right\}$ and $D = \bigcup_{n=0}^\infty D_n$; index our standard normals by $\{N_{q}\}_{q \in D}$, and set $B_0 = 0, B_1 = N_1$, and $B_{1/2} = \frac{B_1 - B_0}{2} + \frac{1}{2}N_{1/2}$. Just continue the same thing for every such dyadic, such that
$$
\{B_{1/2^n} - B_0, B_{2/2^n} - B_{1/2^n}, \dots, B_1- B_{(2^n-1) / 2^n} \}
$$
are all independent $N\left(0, 2^{-n}\right)$. 

**Theorem**: Almost surely, $t \mapsto B_t$, $t \in D$ is uniformly continuous.

_Proof_: Set $K_n = \sup \{ |B_s - B_t| \mid s,t \in D, |s - t| \leq 2^{-n}\}$. We just need to show that $K_n \to 0$ as $n \to \infty$. In fact, something even stronger is true: for $\alpha < \frac{1}{2}$, $\lim_{n \to \infty} 2^{\alpha n}K_n = 0$. Morally speaking, just think that each Brownian increment is about its standard deviation, which is $|t - s|^{1/2}$.

Technically, however, we proceed as follows: set
$$
  Y_n = \max\{B_{1/2^n} - B_0, B_{2/2^n} - B_{1/2^n}, \dots, B_1- B_{(2^n-1) / 2^n} \}
$$
and note that the union bound yields
$$ 
\begin{align*}
  P(Y_n \geq x) &\leq \sum_{j=1}^{2^n} P(|B_{j/2^n} - B_{(j-1)/2^n}| \geq x) \\
      &= 2^nP(B_{1/2^n} \geq x) \\
      &= 2^{n+1}P(B_1 \geq 2^{n/2} x).
\end{align*}
$$
If we choose $x_n$ such that $\sum_{n=1}^\infty 2^{n+1}P(B_1 \geq 2^{n/2}x_n) < \infty$, then by Borel-Cantelli shows that for sufficiently large $n$, $Y_n \leq x_n$ almost surely. Do any reasonable bound you like on the tail of the normal distribution and take a sufficiently large $x$ and call it a day. In particular, if you choose the easiest bound $P(N \geq x) \leq Ce^{-x^2/2}$, you can eventually get the bound:

**Prop**:
$$
  \limsup_{n \to \infty} \frac{2^{n/2}}{\sqrt{n}}Y_n \leq \sqrt{2 \log 2}.
$$

_Proof_: Look at the sum
$$
  \sum_{n=1}^\infty P\left(Y_n > \sqrt{n} \cdot 2 ^{-n/2} \cdot \sqrt{2 \log 2 (1 + \epsilon)}\right)
$$
and apply Borel-Cantelli. In particular, we have
$$
\begin{align*}
  P(Y_n > x_n) &\leq \sum_{j=1^{2^n}} P(|B_{j/2^n} - B_{(j-1)/2^n}| > x_n) \\
  &= 2^{n+1}P(B_{1/2^n} > x_n) \\
  &= 2^{n+1}P\left(B_1 > \sqrt{n} \sqrt{2(\log 2 (1 + \epsilon)}\right) \\
  &\leq C 2^n e^{-\frac{\sqrt{2\log 2 n(1+\epsilon)^2}}{2}} \\
  &\leq C e^{-n\epsilon}.
\end{align*}
$$

**Prop**: Set $K_n = \sup \{ |B_s - B_t| \mid s,t \in D, |s - t| \leq 2^{-n}\}$; then, there is $C$ such that with almost surely,
$$
  \limsup_{n \to \infty} \frac{2^{n/2}}{\sqrt{n}}K_n \leq C.
$$

_Proof_: It's easy to see that $K_n \leq 2 \sum_{j=n+1}^\infty Y_j$ (this is just the triangle inequality). Then for sufficiently large $n$, we get
$$
  K_n \leq 2 \cdot 2 \sum_{j=n+1}^\infty 2^{-j/2}\sqrt{j}
$$
with full probability, and so
$$
  \sup_{\substack{s, t \in D \\ s < t}} \frac{|B_t - B_s|}{\sqrt{(t-s)|\log((t-s)^{-1})|}} < \infty.
$$

Now we may set $B_t$ for $t \in [0, 1]$ by $B_t = \lim_{\substack{s \to t \\ s \in D}}B_s$, and check that this is in fact a genuine Brownian motion, which is not bad; and of course this construction can extend to $[0, \infty)$ easily as well.

###  Properties of Brownian Motion 

_Def_: A function $f: [0, 1] -> \mathbb R$ is called **Hölder continuous** of order $\beta \geq 0$ if there is some $C < \infty$ such that for all $s, t$, $|f(t) - f(s)| \leq C|t-s|^\beta$. Futher, $f$ is **weakly Hölder continuous** of order $\beta$ if it is Hölder continuous of order $\alpha$ for all $\alpha < \beta$. In both cases, we will say Hölder-$\beta$ continuous for short.

**Prop**: Brownian motion paths are weakly Hölder-$\frac{1}{2}$ continuous.

_Proof_: Omitted.

**Theorem**: The function $t \mapsto B_t$ is nowhere differentiable almost surely.

_Proof_: Assume $|f'(t)| < K$; there exists $\delta > 0$ such that if $|s - t| \leq \delta$, $|f(t) - f(s)| \leq 2K|s-t|$; in particular there is an $N$ such that for all $n > N$, $|s - t| \leq n^{-1}, |r - t| \leq n^{-1}$, $|f(s) - f(r)| \leq 4Kn^{-1}$. Then, set 
$$
  Z_{k, n} = \max \left\{ |B_{k/n} - B_{(k-1)/n}|, |B_{(k+1)/n} - B_{k/n}|, |B_{(k+2)/n} - B_{(k+1)/n}|\right\}
$$
and
$$
  Z_n = \min \{ Z_{k, n} \mid k = 1, \dots, n \}.
$$

If $B$ is differentiable, then there is some $M$ such that $Z_n \leq Mn^{-1}$ for all $n$. Now set $E_M$ to be the event that $Z_n \leq Mn^{-1}$ for all sufficiently large; our theorem reduces to showing that $P(E_M)  = 0$ for all $M$. In fact, we will show that
$$
  \lim_{n \to \infty} P(Z_n \leq Mn^{-1}) = 0.
$$

Consider the union bound
$$
\begin{align*}
  P(Z_n \leq Mn^{-1}) &\leq \sum_{j=1}^{n}P(Z(n,k) \leq Mn^{-1}) \\
  &\leq nP\left( \max \left\{ |B_{1/n}|, |B_{2/n} - B_{1/n}|, |B_{3/n} - B_{2/n}|\right\}
\right) \\
&\leq nP(|B_{1/n}| \leq Mn^{-1})^3 \\
&\leq nP(|B_1| \leq Mn^{-1/2})^3
\end{align*}
$$
and just do literally the stupidest estimate you can, e.g. just look at the density and say that the probability is bounded by $2CMn^{-1/2}$, so that the above is sent to zero as $n \to \infty$.


###  Filtrations 

_Def_: A filtration $\{ \mathcal F \}_{t \geq 0}$ is an incerasing collection of sub $\sigma$-algebras. Further, we put
$$
  \mathcal F_{\infty} = \bigcup_{t \geq 0} \mathcal F_t.
$$

_Def_: A stochastic process $\{X_t \}_{t \geq 0}$ is adapted to $\{ \mathcal F_t \}_{t \geq 0}$ if for each $t$, $X_t$ is $\mathcal F_t$-measurable.

_Def_: A process $\{ B_t \}_{t \geq 0}$ is a standard Brownian motion start at 0 w.r.t. $\{ \mathcal F_t \}_{t \geq 0}$ if

- $B_0 = 0$,
- $\{ B_t \}_{t \geq 0}$ is adapted to $\{ \mathcal F_t \}_{t \geq 0}$,
- if $s < t$ then $B_t - B_s$ is independent of $\mathcal F_s$,
- $B_t - B_s \sim N(0, t-s)$,
- and with probability 1 $t \mapsto B_t$ is continuous.

_Def_: A random variable $\tau$ taking values in $[0, \infty]$ is called a stopping time with respect to $\{ \mathcal F_t \}_{t \geq 0}$ if for every $t$, the event $\{\tau \leq t\} \in \mathcal F_t$.

_Examples_: The following are all stopping times:

- constants;
- $\tau = \inf\{ t \mid B_t \in V \}$ where $V$ is Borel;
- $\tau_1 \land \tau_2, \tau_1 \lor \tau_2$, where $\tau_1, \tau_2$ are both stopping times.

_Def_: If $\tau$ is a stopping time, then $\mathcal F_\tau$ is the $\sigma$-algebra corresponding to the collection of events $A$ such that for each $t$, $A \cap \{ \tau \leq t \} \in \mathcal F_t$.

###  The Markov Property of Brownian Motion 

_Def_: For some stochastic process $\{X_t\}_{t \geq 0}$ (or any other indexed set) with filtration $\{F_t\}_{t \geq 0}$, $X_t$ has the **Markov property** if it saitisfies that
$$
  E[f(X_t) \mid \mathcal F_s] = E[f(X_t) \mid \sigma(X_s)].
$$

_Def_: In general, if $X_t$ is a stochastic process and $\tau$ is a stopping time, both adapted to $\{ \mathcal F_t \}_{t \geq 0}$ with $P(\tau < \infty) = 1$, then $X_t$ has the **strong Markov property** if $X_{\tau + t}$ is independent of $\mathcal F_{\tau}$.

**Prop**: Suppose $B_t$ is a Brownian motion and $\tau$ is a stopping time, both with respect to $\{ \mathcal F_t \}$, and assume $P(\tau < \infty) = 1$. Set
$$
  Y_t = B_{\tau + t} - B_{\tau}.
$$
Then $Y_t$  is a Brownian motion independent of $\mathcal F_t$, e.g. Brownian motion has the strong Markov property and the new process is also a Brownian motion.

_Proof_: You proceed by doing successive approximations.

- First, let $\tau$ take a finite amount of values, and use the normal Markov property.
- Then, approximate any $\tau$ by stopping times taking a finite amount of values, such as by
$$
  \tau_n = \begin{cases}
    \frac{k}{2^n} & \frac{k - 1}{2^n} \leq \tau \leq \frac{k}{2^n} \leq n \\
    n & \tau > n
  \end{cases}
$$
  for example.
- Take a limit by continuity.

In particular, the following is clear for Brownian motion:

**Prop**: If $\{B_t\}_{t \geq 0}$ is a Brownian motion, $t$ a fixed time, and $Y_s = B_{t+s} - B_t$, then $\{Y_s\}_{s \geq 0}$, is a BM and independent of $\mathcal F_t = \sigma\{B_s \mid s \leq t\}$.

**Theorem**: Set $B_t$ to a Brownian motion with drift zero; **the reflection principle** is that
$$
  P\left(\max_{0 \leq s \leq t} B_s \geq a\right) = 2P(B_t \geq a).
$$
Set $\tau_a = \min \{ s \mid B_s = a \}$; we also have
$$
  P(\tau_a \leq t) = 2P(B_t \geq a)
$$
or equivalently
$$
  P(B_t \geq a \mid \tau_a \leq t) = \frac{1}{2}.
$$

**Prop**: If $0 < r < s < \infty$,
$$
  q(r, s) = P(B_t = 0 \text{ for some } r \leq t \leq s) = 1 - \frac{2}{\pi} \arctan\left(\sqrt{\frac{r}{s - s}}\right).
$$

_Proof_: First, I claim that $q(r, s) = q(1, s/r)$ just by a change of variables, so we only need to compute $q(t) = q(1, 1 + t)$. Set $A = \{ B_s = 0 \text{ for some } 1 \leq s \leq 1 + t \}$, so that
$$
\begin{align*}
  q(t) &= \frac{2}{\sqrt{2 \pi}}\int_0^\infty P(A \mid B_1 = x) e^{-x^2 / 2}dx.
\end{align*}
$$
However, by the reflection principle, we have that
$$
\begin{align*}
  P(A \mid B_1 = x) &= P\left(\max_{0 \leq s \leq t} B_s \geq x \right) \\
  &= P \left( \min_{0 \leq s \leq t} B_s \leq -x \right) \\
  &= 2P(B_t \geq x) = 2P\left(B_1 \geq \frac{x}{\sqrt{t}}\right)
\end{align*}
$$
upon which we can just compute the integral.

**Corollary**: One dimensional standard Brownian motion is (pointwise) recurrent (that is, the zero set of Brownian motion is unbounded).

**Corollary**: Since $Y_t = t^{-1}B_{1/t}$ is a standard Brownian motion, this shows that for any $\epsilon > 0$, $Z_\epsilon = \{ t \mid B_t = 0, 0 \leq t \leq \epsilon \}$ has more elements that just $0$.

###  Martingales 

_Def_: A process $\{ M_t \}_{t \geq 0}$ is a **supermartingale** (resp. **submartingale**) w.r.t. $\{ \mathcal F_t \}$ if
- $E[|M_t|] < \infty$,
- $M_t$ is $\{ \mathcal F_t \}$ adapted,
- and if $E[M_t \mid \mathcal F_s] \leq M_s$ (resp. $\geq M_s$) almost surely for $s \leq t$.
A process which is both a submartingale and supermartingale is just a martingale, it is continuous if $M_t$ is a continuous function of $t$ almost surely, and is square integrable (or simply $L^2$) if it has finite second moment for all $t$.

As an aside, I'm going to stop saying with probability 1 or almost surely because it's annoying!

**Prop**: Brownian motion (without drift) is an $L^2$ continuous martingale.

**Theorem (Kolmogorov Zero-One)**: Tail events of independent $\sigma$-algebras happen with probability 0 or 1.

_Proof_: Set $\mathcal F_n = \sigma(X_1, \dots, X_n)$, and $\mathcal T_n = \sigma(X_{n+1}, \dots)$, $\mathcal F_\infty = \bigcup_{n} \mathcal F_n$, and $\mathcal T_\infty = \bigcap_n \mathcal T_n$.

Then, one can see that if $A \in \mathcal F_\infty$ and $\epsilon > 0$, there is $n$ and $A_n \in \mathcal F_n$ such that $P(A_n \Delta A) < \epsilon$; even better, there exists $A_n$ independent of $A \in \mathcal T_\infty$ such that $P(A \Delta A_n) < \epsilon$, and so we conclude that $P(A) = P(A)P(A)$.

**Theorem (Blumenthal Zero-One)**: Let $B_t$ be a Brownian motion with the standard filtration, and set $\mathcal F_{0+} = \bigcap_{\epsilon > 0} \mathcal F_{\epsilon}$; then, if $A \in \mathcal F_{0+}$, either $P(A) = 0$ or $1$.

###  Quadratic Variation 

Let $B_t$ be a standard Brownian motion; a partition $\Pi$ of $[0, 1]$ is a sequence $0 = t_0 < t_1 < \dots < t_k = 1$, and the mesh of the partition is just
$$
  \| \Pi \| = \max_{i = 1, \dots, k} t_i - t_{i - 1}.
$$
Now take a sequence of partitions $\Pi_n$, e.g. $0 = t_{0, n} < \dots < t_{k_n, n}$, and define the quantity
$$
  Q(t, \Pi) = \sum_{t_{j} <= t} (B_{t_j} - B_{t_{j-1}})^2
$$
and $Q_n(t) = Q(t, \Pi_n)$ and $Q_n = Q_n(1)$. 

**Theorem**: If $\| \Pi_n \| \to 0$, then $Q_n \to 1$ in probability. Furthermore, if $\sum_{n=1}^\infty \| \Pi_n \| < \infty$, then almost surely $\lim_{n \to \infty} Q_n = 1$.

_Proof_: A simple computation gives us that 
$$
  E(Q_n) = \sum_{j=1}^{k_n}E[(B_{t_j} - B_{t_{j-1}})^2] = \sum_{j=1}^{k_n}(t_j - t_{j-1}) = 1
$$
and
$$
  \Var(Q_n) = \sum_{i=1}^{k_n} \Var((B_{t_j}  - B_{t_{j-1}})^2) = \sum_{j=1}^{k_n} (t_j - t_{j-1})^2 \Var(B_1^2) = 2\sum_{j=1}^{k_n}(t_j - t_{j-1})^2.
$$
Then,
$$
  \Var(Q_n) \leq \| \Pi_n \| 2 \sum_{j=1}^{k_n}(t_j - t_{j-1}) = 2 \| \Pi_n \|
$$
and
$$
  P(|Q_n - 1| \geq \epsilon) \leq \frac{\Var(Q_n)}{\epsilon^2} \leq \frac{2 \| \Pi_n \|}{\epsilon^2}.
$$

The latter half of the theorem follows from Borel-Cantelli.

**Theorem**: In general, if $\sum_{n=1}^\infty \| \Pi_n \| < \infty$, then almost surely we have $Q_n(t) \to t$.

_Proof_: With probability 1 this holds for rational $t$, but by construction $Q_n(t)$ is monotone, so it holds everywhere.

_Def_: In general,  if $X_t$ is a process, then its quadratic variation is
$$
  \left\langle X\right\rangle_t = \lim_{n \to \infty} \sum_{t_{j,n} \leq t} (X_{t_{j, n}} - X_{t_{j-1, n}})^2
$$
(sort of, since sometimes this depends on the partition).

_Def_: Alternatively, if $M_t$ is a continuous $L^2$-martingale, its quadratic variation is the unique predictable process $\left\langle M \right\rangle_t$ that makes $M_t^2 - \left\langle M \right\rangle_t$ a martingale (in particular, this exists by Doob decomposition).

We showed above that $\left\langle B_t \right\rangle_t = t$.

**Prop**: If $B_t$ is a standard Brownian motion and $Y_t = \mu t + \sigma B_t$, then
$$
  \left\langle Y \right\rangle_t = \sigma^2 t.
$$

_Proof_: Just check directly.

###  Law of the Iterated Logarithm 

**Lemma (Relaxed Borel-Cantelli)**: Let $A_1, A_2, \dots$ be a sequence of events, and set $\mathcal F_n = \sigma(A_1, \dots, A_n)$; if there is $q_n$ with
$$
  \sum_{n=1}^\infty q_n = \infty
$$
such that $P(A_n \mid \mathcal F_{n-1}) \geq q_n$, then $A_n$ happens infinitely often almost surely.

_Proof_: Same as usual, but just with a little more caution.

Recall that the second Borel-Cantelli lemma tells us that if they are independent and $\sum_{n=1}^\infty P(A_n) = \infty$, then $A_n$ occurs infinitely often with probability 1; this lemma is stronger than that.


**Theorem**: If $B_t$ is a standard Brownian motion, then
$$
  \limsup_{t \to \infty} \frac{B_t}{\sqrt{2 t \log \log t}} = 1.
$$

**Corollary**: By symmetry,
$$
  \liminf_{t \to \infty} \frac{B_t}{\sqrt{2 t \log \log t}} = -1.
$$

_Proof_: Let $\mathcal T_t = \sigma\{ B_{s + t} - B_t \mid s \geq 0\}$, and $\mathcal T_\infty = \bigcap_t \mathcal T_t$;  one can adapt the arguments from the Kolmogorov 0-1 law to show that everything in $\mathcal T_\infty$ happens with probability 0 or 1. Then,
$$
  A_\epsilon = \left\{ \omega \mid \limsup_{t \to \infty} \frac{B_t}{\sqrt{2 t (1 + \epsilon) \log \log t}} \leq 1\right\}
$$
is a tail event (e.g. in $\mathcal T_\infty$) and thus $P(A_\epsilon) = 0$ or $1$. In fact, if $\epsilon < 0$ then it's 0, and if $\epsilon > 0$ then it's 1. In fact, by this 0-1 law and symmetry one may see
$$
  P \left( \limsup_{t \to \infty} \frac{|B_t|}{\sqrt{2 t (1 + \epsilon) \log \log t}} \leq 1 \right) = P \left(\limsup_{t \to \infty} \frac{B_t}{\sqrt{2 t (1 + \epsilon) \log \log t}} \leq 1\right)
$$
as well.

First take $\epsilon > 0$ and take some $\rho > 1$ to be specified later. Then, let 
$$
  V_n = \left\{ |B_{\rho^n}| \geq \sqrt{2 \rho^n(1-\epsilon) \log \log \rho^n} \right\};
$$
we want to to show that $V_n$ occurs infinitely often. I claim that
$$
  P(V_{n+1} \mid V_1, \dots, V_n) \geq P \left( B_{\rho^{n+1}} - B_{\rho^n} \geq \sqrt{2 \rho^{n+1}(1-\epsilon) \log \log \rho^n} \right).
$$
To see this, compute
$$
\begin{align*}
  &P \left( \frac{B_{\rho^{n+1}} - B_{\rho^n}}{\sqrt{\rho^{n+1} - \rho^n}} \geq \frac{\sqrt{2\rho^{n+1}(1-\epsilon)\log \log \rho^n}}{\sqrt{\rho^n (\rho - 1)}} \right) \\
  &= P \left( \frac{B_{\rho^{n+1}} - B_{\rho^n}}{\sqrt{\rho^{n+1} - \rho^n}} \geq \sqrt{\frac{2 \rho}{\rho - 1}(1-\epsilon)(\log n + \log \log \rho)} \right) \\
\end{align*}
$$
and choose $\rho$ large enough so that $\frac{2 \rho}{\rho - 1}(1-\epsilon) < 1$, and use the estimate $P(B_1 > x) \sim \exp(-x^2 / 2)$ and conclude by the earlier lemma. The other direction for $\epsilon < 0$ is similar (in fact, easier since we may conclude from the first Borel-Cantelli lemma).

###  Zero Sets of Brownian Motion 

_Def_: Set $B_t$ a standard Brownian motion, $Z = \{ t \mid B_t = 0 \}$, and $Z_t = Z \cap [0, t]$. Then, $t \in Z$ is right-isolated if $t \in Z$, and $\exists \epsilon > 0$ such that $(t, t + \epsilon) \cap Z = \emptyset$; similar for left-isolated. A point which is both left and right isolated is just isolated.

With probability $1$, $0$ is not right-isolated; further, $Z_1$ is homeomorphic to the Cantor set.

With probability 1, the sets of left and right isolated points are countable, and there are no isolated points.

We can show that if $q \in \Q_{\geq 0}$, $P(B_q  = 0) = 0$, and $\tau_q = \min \{ t \geq q \mid B_t = 0\}$ is a stopping time; further, the left-isolated points are just $\{ \tau_q \mid q \in \Q_{\geq 0}\}$. And by strong Markov property, no $\tau_q$ is a right-isolated point, so there are no isolated points.

Set $\sigma_q = \max_t \{ t < q \mid B_t = 0 \}$ (this is well-defined, but not a stopping time). Then associate every $q$ to an interval $(\sigma_q, \tau_q)$; then $Z = [0, \infty) \setminus \bigcup_q (\sigma_q, \tau_q)$ and has Lebesgue measure 0. To see this, we just interchange integrals:
$$
  E[ \lambda(Z_1) ] = E\left[ \int_0^1 1_{\{B_s = 0\}} ds \right] = \int_0^1 P\{ B_s = 0 \}ds = 0.
$$

Look at $Z \cap [1, 2]$; now cover $[1, 2]$ by intervals of length $n^{-1}$; put the number of these intervals that intersect $Z$ as $X_n$. Then,
$$
  E[X_n] = \sum_{j=1}^n P\left(Z \cap \left[1 + \frac{j-1}{n}, 1 + \frac{j}{n}  \right]\neq \emptyset \right) \sim Cn^{1/2}.
$$

The box or Minkowski dimension of a set is given by the exponent above (sort of).

_Def_: The **Hausdorff dimension** comes from the **Hausdorff measure**: for $\epsilon > 0$ and $\alpha$, set
$$
  \mathcal H^\alpha_\epsilon = \left\{ \inf \sum_{j=1}^\infty (\operatorname{diam} U_j)^\alpha \right\}
$$
where the infimum is over all coverings $\bigcup_j=1^\infty U_j$ of $V$ with each $\operatorname{diam} (U_j) < \epsilon$. Then the Hausdorff measure is given by
$$
  \mathcal H^\alpha(V) = \lim_{\epsilon \to 0} \mathcal H_{\epsilon}^\alpha (V).
$$
Then, there is some number $D$ such that
$$
  \mathcal H^{\alpha}(V) = \begin{cases}
    \infty & \alpha < D \\
    0 & \alpha > D
  \end{cases}
$$
and we call this $D$ the Hausdorff dimension. In general, the Hausdorff dimension is at most the Minkowski dimension, but in this case we actually do have equality.

###  Local Time 

The local time is the amount of time the Brownian motion spends at $0$ by a certain time. It is sort of like the Cantor measure, insofar as if $s < t$, then $L_t - L_s > 0 \iff (s, t) \cap Z \neq \emptyset$. 

_Def_: We define the **local time** of Brownian motion at $x \in \R$, $L_t$, by first setting
$$
  L_{t, \epsilon}(x) = \frac{1}{2\epsilon} \int_0^t 1_{\{|B_s - x| \leq \epsilon\}} ds
$$
and then letting $L_t(x) = \lim_{\epsilon \to 0} L_{t, \epsilon}$. We abbreviate $L_t = L_t(0)$.

We can compute the expectation
$$
  \begin{align*}
    E[L_{t, \epsilon}] &= \frac{1}{2\epsilon} \int_0^t P(|B_s| \leq \epsilon)ds \\
    &\sim \frac{1}{2\epsilon} \int_0^t 2e \frac{1}{2 \pi s}ds \\
    &= \sqrt{\frac{2}{\pi}}t^{1/2}
  \end{align*}.
$$

**Theorem**: With probability 1, $L_t$ exists for all $t$, and this holds in $L^2$ as well.

There are more facts: $L_t$ is continuous in $t$ and non-decreasing, $L_t - L_s > 0 \iff (s, t) \cap Z \neq \emptyset$, and $L_t$ is weakly $\frac{1}{2}$-Hölder continuous.

**Theorem (Scaling Rule)**: $L_t$ has the same distribution as $t^{1/2} L_1$. Further, $M_t = \max_{ 0 \leq s \leq t}B_s$ has the same distribution as $L_t$.

##  Brownian Motion in Several Dimensions 

-------

_Def_: If $B_t^1, \dots, B_t^d$  are independent standard Brownian motions, then
$$
  B_t = (B_t^1, \dots, B_t^d)
$$
is a standard Brownian motion in $\R^d$.

**Prop**: $B_0 = \mathbf 0$, has independent increments, if $s < t, B_t - B_s \sim N(0, (t-s)I)$, and $B_t$ has continuous paths.

In particular, if we have the above with $B_t - B_s \sim N(\mu, \Gamma)$, then $B_t$ is a Brownian motion with dift $\mu \in \R^d$ and covariance matrix $\Gamma$. Further, if $AA^T = \Gamma$, then we may write $B_t = AY_t + t \mu$ for a standard Brownian motion $Y_t$.

Consider the open annulus with inner radius $r$ and outer radius $R$, e.g. $D(r, R) = \{ y \in \R^d \mid r < |y| < R \}$. Start a Brownian motion at $x$, and let $\tau = \tau(r, R)$ be the first time with $|B_t| = r$ or $|B_t| = R$. What is the probability that $|B_\tau| = R$?

In one dimension, this is easy: stop the Brownian motion at $\tau$ and look at what happens at infinity.

In higher dimensions, we need a quick detour.

###  Harmonic Functions in $\R^d$ 

For this section, a domain is an open connected subset of $\R^d$.

_Def_: For a domain $D$, we say $f: D \to \R$ is harmonic if it is continuous (or merely locally integrable) and satisfies the mean value property
$$
  f(x) = MV(f, x, \epsilon) = \int_{|x - y| = \epsilon} f(y) ds
$$
where $s$ is the surface measure, normalized so that $\int_{|x - y| = \epsilon}ds = 1$.

In a probabilistic vein, if $B_t$ is a standard $d$-dimensional Brownian motion starting at $x$, and $\tau = \min \{ t \mid |B_t - x| = \epsilon \}$, since $B_t$ is rotation invariant, the above is just
$$
  MV(f, x, \epsilon) = E[f(B_\tau)].
$$

_Def_: The Laplacian of $f$ is 
$$
  \nabla f = \sum_{j=1}^d \frac{\partial ^2 f }{\partial x_j^2}.
$$

**Prop**: If $f$ is $C^2$ in $D$, then 
$$
  \frac{1}{2d}\nabla f = \lim_{\epsilon \to 0} \frac{MV(f, x, \epsilon) - f(x)}{\epsilon^2}.
$$

_Proof_: Look at the Taylor expansion.

**Theorem**: $f$ is harmonic in $D$ if and only if it is in $C^2$ and $\nabla f = 0$ everywhere.

###  Hitting Probabilities for Brownian Motion 

Let $\tau = \tau_{r, R} = \min \{ |B_t| = r \text{ or } R \} = \min \{ t \mid B_t \in \partial D \}$. We will let a superscript $x$ denote that $B_0 = x$, and let
$$
  \phi(x) = P^x(|B_\tau| = R).
$$
By rotational invariance, $\phi(x) = \phi(|x|)$, and it is continuous at the boundary; we also clearly have $\phi(x) = 0$ if $|x| = r$ and $\phi(x) = 1$ if $|x| = R$. The strong Markov property furnishes that $\phi$ is harmonic. Set $\tau_\epsilon = \min \{ t \mid |B_t - x| = \epsilon\}$; then
$$
  P^X(|B_\tau| = R \mid \mathcal F_{\tau_\epsilon}) = \phi(B_{\tau_\epsilon})
$$
and
$$
  P^X(|B_\tau| = R) = E^x[P(|B_\tau| = R \mid \mathcal F_{\tau_\epsilon})] = E^x[\phi(B_{\tau_\epsilon})] = MV(\phi, x, \epsilon).
$$
This list of properties gives a unique function (solve an ODE in polar coordinates), and is given by
$$
  \phi(x) = \frac{|x|^{2-d} - r^{2-d}}{R^{2-d} - |r|^{2-d}} 
$$
for $d \neq 2$ and
$$
  \phi(x) = \frac{\log|x| - \log r}{\log R - \log r}
$$
for $d = 2$.


###  Recurrence and Transience of Brownian Motion 

Let $B_t$ be a standard Brownian motion in $\R^d$.

**Theorem**: With probability 1, Brownian motion is transient in dimensions at least 3; that is,
$$
  \lim_{t \to \infty} |B_t| = \infty.
$$

_Proof_: If $d \geq 3$, and we start at $x$ with $|x| > r$, if we set $T_r = \min \{ t \mid |B_t| = r \}$,
$$
  P^x(T_r < \infty)  = \lim_{R \to \infty} P^X(|B_{\tau_{r, R}}| = r) = \left( \frac{r}{|x|} \right)^{d-2} < 1
$$

**Theorem**: With probability 1, Brownian motion is neighborhood recurrent; that is, $\forall z \in \R^2, \epsilon > 0$, Brownian motion visits the disk of radius $\epsilon$ about $z$ infinitely often. However, it is not point recurrent, e.g. it hits $z \neq 0$ with probability 0.

_Proof_:  If $d = 2$, then
$$
  P^x(T_r < \infty) = \lim_{R \to \infty} P(|B_{\tau_{r, R}}| = r) = 1
$$
But if $T = \min\{ t \mid B_t = 0 \}$,
$$
  P^x(T < \infty) \leq \lim_{R \to \infty} \lim_{r \to 0} P^x(|B_{\tau_{r, R}}| = r) = 0
$$

A fun fact is that if $d \geq 2$, then $\{ B_t, t \geq 0 \}$ has Hausdorff dimension $2$ but also zero Hausdorff-2 measure.

###  The Dirchlet Problem 

Take a bounded domain $D \subset \R^d$, and a continuous function $F: \partial D \to \R$; the Dirichlet problem is to find the unique continuous $f: \overline D \to \R$ that agrees with $F$ on $\partial D$ and is harmonic on $D$.

In fact, uniqueness follows from the maximum principle, which says that the maximum of $f$ is attained on the boundary (think about the mean value principle). Then subtract two solutions and see that it is 0.

Let $T = \min\{t \mid B_t \in \partial D \}$, and $f(x) = E^x[F(B_T)]$. This satisfies the mean value principle, and in fact continuous (which is kinda hard) but it is obviously locally integrable, so $f$ here is harmonic.

In general, such a harmonic function is not necessarily continuous on the boundary: take the example of the punctured unit disk, with $F(x) = 1$ for $|x| = 1$ and $F(0) = 0$. And so everything is fine except at the origin, since $f(x) = P^x(|B_T| = 1) = 1$.

_Def_: If $x \in \partial D$, let $\sigma = \inf \{ t > 0 \mid B_t \in \partial D\}$; $x$ is regular if $P^x(\sigma = 0) = 1$.

**Prop**: $f$ as defined above is continuous at every regular boundary point. If we relax the boundary to be only bounded and measurable, then the above is continuous at every regular point at which $F$ is continuous.

Therefore the Dirichlet problem has a solution for every continuous $F$ if and only if every point on $\partial D$ is regular.

_Def_: If $D$ is a domain, then **harmonic measure** on $D$ (or on $\partial D$) at $z \in D$, is the hitting measure of Brownian motion, starting at $z$ and stopped at $\partial D$. That is, if $V \subset \partial D$,
$$
  \hm_D(V, z)  = P^z(B_T \in V).
$$

**Prop**: Note that $\hm_D(\partial D, z) = P^z(T < \infty)$ so if $D$ is bounded, this is a probability measure. Then,
$$
  E^z[F(B_T)] = \int_{\partial D} F(w) \hm_D(dw, z)
$$
and if $\partial D$ is smooth, then $\hm_D(\cdot, z)$ is absolutely continuous with respect to surface measure; if $H_D(z, w)$ is the Poisson kernel,
$$
  \hm_D(v, z) = \int_V H_D(z, w)S(dw)
$$
where $S$ is the surface measure.

A nice example of all of the above is the unit ball. Then, $H_D(z, w) = c_d^{-1}\frac{1 - |z|^2}{|z - w|^d}$ where $c_d$ is the surface measure of the unit sphere.

**Prop**: The following are true.

- If $V \subset \partial D$ and $h(z) = \hm_D(V, z) = P^z(B_T \in V)$, then $h$ is the unique harmonic function in $D$ with boundary condition 
$$
  F(w) = \begin{cases}
    1 & w \in V \\
    0 & w \notin V
  \end{cases}
$$
In particular, one could define the harmonic measure this way and show existence, but this is unreasonably hard compared to just using the Brownian motion.
- If $D$ contains the closed unit ball $B$, and $f$ is harmonic in $D$, then for every $|z| < 1$,
$$
  f(z) = \int_{|w| = 1} f(w) H_B(z, w)S(dw).
$$

**Corollary (Harnack Inequality)**: For every $r \in (0, 1)$, there is $C = C(r, d)$ such that if $f: B \to (0, \infty)$ is harmonic, then for all $|z|, |z'| \leq r$,
$$
  f(z) \leq Cf(z')
$$
and
$$
  C = \max_{\substack{|z|, |z'| \leq r \\ |w| = 1}} \frac{H_B(z, w)}{H_B(z', w)}
$$

**Corollary**: For every $k$ there exists $c = c(k, d)$ such that if $f: B \to \R$ is harmonic, then for any $k$-th order partial derivative,
$$
  |\partial^k f(0)| \leq c\| f \|_\infty.
$$
Moreover, there is $C = C(k, d)$ such that if $f: D \to \R$ is harmonic and $z \in D$, then
$$
  |\partial^k f(z)| \leq \frac{C_k}{(\operatorname{dist}(z, \partial D)^k} \left( \max_{|z - w| \leq \operatorname{dist}(z, \partial D) / 2} |f(w)| \right).
$$

**Theorem (Harnack principle)**: if $D$ is a domain and $K \subset D$ is compact, then there exists some $C = C(K, D)$ such that if $f: D \to (0, \infty)$ is harmonic and $z, w \in K$, then
$$
  f(z) \leq Cf(w).
$$

Now let $D$ and $\R^{d} \setminus D = K$; if $F: K \to \R$ is continuous, we can ask about the existence of $f: \R^d \to \R$ which is harmonic on $D$, coincides with $F$ on $K$, and is continuous on $\R^d$. Such a thing is obviously not unique (think!).

Let both $F$ and $f$ bounded but otherwise as above. Set $T = \min\{ t \geq 0 \mid B_t \in K \}$, and assume for every $z \in D$, $P^z(T < \infty) > 0$. If $d = 1$ or $2$, there is a unique $f$, given by $f(z) = E^z[f(B_T)]$. 

Now if $d \geq 3$ and $K$ is bounded, then $g(z) = P^z(T = \infty) > 0$; $g$ is harmonic and bounded, and if $g = 0$ on $K$ then $g$ is continuous.

**Theorem**: All solutions to the above problem are given by
$$
  f(z) = E^z[(F(B_t) 1_{\{T < \infty\}}] + cP^z(T = \infty).
$$

Intuitively, we just add a point at infinity which has value $c$ and get this formula; but this only works on $\R^d$ and $\Z^d$.

Consider a random walk on an infinite binary tree; such a walk is clearly recurrent, since it moves away from the starting point with probability $2/3$. But now there are lots of infinities.

## Differential Equations

The heat equation is described by some initial function $u_0: D \to \R$, some boundary condition $u(t,x) = F(x)$ for $x \in \partial D$, and $\dot i(t,x) = \frac{1}{2} \Delta u(t,x)$. Then,
$$
    u_t(x) = E^x[1_{\{T \leq t\}}F(B_T) + u_0(B_t)1_{\{T > t\}}].
$$

We can also handle Green's functions. Let $B_t$ be a Brownian motion in $\R^d$, $d \geq 3$; $G(x,y)$ is the "expected amount of time spent at $y$ starting at $x$; that is,
$$
    G(x, y) = \int_0^\infty P_t(x,y)dt = G(y,x) = G(0, y-x).
$$
Further,
$$
    G(x) = G(0, x) = \int_0^\infty \frac{1}{(2\pi t)^{d/2}}\exp(-|x|^2/2t)dt;
$$
and when one computes this integral, we get $C_d |x|^{2-d}$ away from the origin. This is no coincidence: in $d=2$, you get $a(x) = C_2\log(|x|)$, which is also the unique radially symmetric harmonic function.

In a domain $D$ (in any dimension), let
$$
    G_D(x,y) = \int_0^\infty P_t^D(x,y)dt
$$
so in $d \geq 3$,
$$
    G(x,y) - E^x[G(B_T, y)]
$$
and in $d = 2$,
$$
    E^x[a(B_T, y)] - a(x,y)
$$
where $a(x,y) = a(y - x)$.

For fixed $y$, the function $h(x) = G_D(x,y)$ is harmonic in $D \setminus \{ y \}$ as $x \to y$.

## Stochastic Integration

Let $\frac{d}{dt} F(t) = C(t, F(t))$; we then write $dF(t) = C(t, F(t))dt$, so $F(t) = F(0) + \int_0^t C(s, F(s))ds$. Here we are interested in the case where we allow things to be random, e.g.
$$
    dX_t = R_tdt + A_tdB_t.
$$
Intuitively, we require that $X_t$ looks like a BM with some drift $R_t$ and variance $A_t^2$. Then,
$$
    X_t = X_0 + \int_0^tR_sds + \int_0^tA_sdB_s.
$$
Of course, we still need to define the (Itô) stochastic integral.

Let $B_t$ be a Brownian motion with a filtration $\mathcal F_t$. 

_Def_: A process $A_t$ is a **simple process** (with respect to $\mathcal F_t$) if there exists a finite number of times $0 = t_0 < t_1 < \dots < t_n < t_{n+1} = \infty$ and random variables $Y_0, Y_1, \dots, Y_n$ such that $Y_j$ is $\mathcal F_{t_j}$-measurable, and $A_t = Y_j$ on $t_j \leq t < t_{j+1}$. We can have $L^2$ or bounded simple processes, and those are just requirements on the $Y_i$.

_Def_: If $A_t$ is a simple process, we define the stochastic integral of $A_t$ to be 
$$
    Z_t = \int_0^t A_sdB_s = \sum_{k=0}^{j-1}Y_k[B_{t_{k+1}} - B_{t_k}] + Y_j[B_t -  B_{t_j}].
$$

We have some properties immediately:

- If $A_t$ and $C_t$ are simple, and $a, c \in \R$, then $aA_t + cC_t$ is simple and
$$
    \int_0^t (aA_t + cC_t)dB_s = a\int_0^tA_tdB_s + c\int_0^tC_tdB_s.
$$
- If $A_t$ is $L^1$, then $Z_t$ is a martingale (just compute).
- The Itô isometry is, if $A_t$ is in $L^2$,
$$
    \Var(Z_t) = E[Z_t^2] = \int_0^t E[A_s^2]ds = E \left[ \int_0^t A_s^2ds \right].
$$
Again, just compute and use the orthogonality of martingale increments.
- With probability 1, $t \mapsto Z_t$ is a continuous function.

Now let $A_t$ be a bounded, continuous process, adapted to $\mathcal F_t$.
_Def_: We define
$$
    \int_0^t A_sdB_s = \lim_{n \to \infty}\int_0^t A_s^{(n)} dB_s
$$
where $A_s^{(n)}$ is a sequence of simple processes approaching $A_s$.

**Lemma**: For every $t$, we can find a sequence of simple processes $A_t^{(n)}$ with $|A_t^{(n)} \leq K$ and such that
$$
    \lim_{n \to \infty} \int_0^t E[(A_s^{(n)} - A_s)^2]ds = 0.
$$

_Proof_: For ease, take $t = 1$. Then,
$$
    A_t^{(n)} = n \int_{\frac{k-1}{n}}^{\frac{k}{n}}A_sds
$$
does the job. Scale appropriately.

Now if $Z_t^{(n)} = \int_0^t A_s^{(n)}dB_s$, then $E[(Z_t^{(n)} - Z_t^{(m)})^2] = \int_0^tE[(A_s^{(n)} - A_s^{(m)})^2]ds$, so $Z_t^{(n)}$ is a Cauchy sequence in $L^2(\Omega)$, and there is therefore an $L^2$ limit that we call $Z_t$. Further, there is a continuous modification of $Z_t$, which we can show exists by defining it on dyadics and using continuity.

**Prop**: Let $Z_t = \int_0^t A_sdB_s$, with $A_s$ continuous in $L^2$. Then the above properties all hold.

We do not necessarily need continuity. We do need progressive measurability, e.g. $\{A_s(\omega)\}$ is measurable on $\Omega \times [0, t]$.

If $A_s$ is continuous, however, without any other boundedness assumptions, let $\tau_n = \min \{ |A_s| \geq n \}$. Then let $A_s^{(n)} = A_{s \land \tau_n}$ so that $Z_t^{(n)} = \int_0^t A_s^{(n)}dB_s$ is well-defined, and moreover if $n > \sup_{0 \leq s \leq t} |A_s|$, then $Z_t^{(n)} = \int_0^t A_sdB_s$. Thus, we define
$$
    Z_t = \int_0^t A_sdB_s = \lim_{n \to \infty}Z_s^{(n)}
$$
which is a pointwise limit in $\Omega$ (so we have to be careful). With this definition, linearity, continuity, and the Itô isometry holds (allowing for $\infty$ as a possible value), but $Z_t$ is merely a local martingale.

**Prop**: If $M_t$ is a continuous square integrable martingale starting from 0, the for all $\epsilon > 0$,
$$
    P \left(\max_{0 \leq s \leq t}|M_s| \geq \epsilon\right) \leq \frac{E(M_t^2)}{\epsilon^2}.
$$

Now, set $M_t = Z_t - Z_t^{(n)}$  is a continuous square-integrable martingale; thus by the above proposition,
$$
    P(\max_{0 \leq t \leq 1} |M_t| \geq \epsilon ) \leq \frac{E(M_t^2)}{\epsilon^2}.
$$

**Theorem**: Suppose $A_s$ is a continuous process with $\int_0^1 E[A_s^2)ds < \infty$, and $\Pi^{(n)}$ is a sequence of partitions of $[0, 1]$ such that
$$
    \sum_{n=1}^\infty \int_0^1 E[|A_t - A_t^{(n)}|^2]dt < \infty.
$$
If $Y^{(n)} = \max_{0 \leq t \leq 1}|Z_t - Z_t^{(n)}|$, then $Y^{(n)} \to 0$ with probability 1.

_Proof_: Apply Borel-Cantelli to the above.

If $Z_t = \int_0^t A_sdB_s$, then the quadratic variation is
$$
    \left\langle Z\right\rangle_t = \int_0^t A_s^2ds
$$
and $Z_t - \left\langle Z \right\rangle_t$ is a martingale. In fact, the quadratic variation is the unique increasing process such that it is a martingale.

### Ito's Formula

**Theorem**: Suppose that $f: \R \to \R$ is $C^2$, and $B_t$ is a standard Brownian motion. Then 
$$
    f(B_t) - f(B_0) = \int_0^t f'(B_s)dB_s + \frac{1}{2}\int_0^t f''(B_s)ds.
$$

_Proof_: Taylor approximate; prove for a countable dense set of $t$, and conclude for $f$ with compact support. Otherwise, take something that agrees with $f$ on $[-n. n]$ send $n \to \infty$.

People often write this in differential form:
$$
    dX_t = R_tdt + A_tdB_t.
$$

Suppose that $f(t, x)$ is a function from $[0, \infty) \times \to \R$ that is $C^1$ in $t$ and $C^2$ in $x$. Then,
$$
f(t, B_t) - f(0, B_0) = \int_0^t \partial_s f(s,B_s)ds + \int_0^t \partial_x f(s, B_s) + \frac{1}{2} \int_0^t\partial^2_{xx}f(s, B_s)ds
$$

Suppose $f$ is $C^2$, but not all of $\R$, merely on an open interval $(a, b)$, and $B_t$ is a Brownian motion starting in the interval. Then, if $T = \inf \{ t \mid B_t = a, b\}$, then if $t < T$, then Ito's formula holds.

_Def_: Suppose that $dX_t = R_tdt + A_tdB_t$; that is,
$$
    X_t = X_0 + \int_0^t R_sds + \int_0^t A_sdB_s.
$$
Then, we define
$$
    \int_0^t Y_sdX_s = \int_0^t Y_sR_sds + \int_0^t Y_sA_sdB_s
$$
and
$$
    \left\langle X \right\rangle_t = \lim_{n \to \infty} \sum_{j < t/n} (X_{j/n} - X_{(j-1)/n})^2.
$$

If you look closely at $X$, you can see that the quadratic variation becomes
$$
    \left\langle X \right\rangle_t = \lim_{n \to \infty} (A_{j/n} - A_{(j-1)/n})^2 = \int_0^t A_s^2ds.
$$

**Prop**: If $f(t,x)$ is $C^1$ in $t$ and $C^2$ in $x$, and $X_t$ is as above,
$$
    f(t, X_t) - f(0, X_0) = \int_0^t \partial_t f(s,X_s)ds + \int_0^t \partial_x f(s, X_s)dX_s + \frac{1}{2} \int_0^t \partial_{xx}f(s, X_s)A^2_sds.
$$

## Diffusions

_Def_: **Diffusions** are SDEs of the form
$$
    dX_t = m(t, X_t)dt + \sigma(t, X_t)dB_t
$$
where $m, \sigma$ are deterministic. Something satisfies the above if
$$
    X_t = y_0 + \int_0^t m(s, X_s)ds + \int_0^t \sigma(s, X_s)dB_s.
$$
for some initial condition $y_0$.

**Prop**: When $m, \sigma$ are uniformly Lipschitz in the latter arguments, a solution to the above exists.

_Proof_: Let $m, \sigma$ be uniformly Lipschitz with constant $\beta$. For ease, assume that there is no $s$-dependence (you don't need this, but writing it is annoying). We proceed by Picard iteration.

Let $X_t^{(0)} = y_0$. Define $X_t^{(n)}$ by
$$
    X_t^{(n)} = y_0 + \int_{0}^t m(X_t^{(n-1)})ds + \int_0^t \sigma(X_s^{(n-1)})dB_s.
$$

We wish to show $X_t = \lim_{n \to \infty} X_t^{(n)}$ exists (in $L^2$). Look at
$$
\begin{align*}
    E[|X_t^{(k+1)} - X_t^{(k)}|^2] &= E \left[ \left| \int_0^t (m(X_s^{(k)}) - m(X_s^{(k-1)}))ds + \int_0^t (\sigma(X_s^{(k)}) - \sigma(X_s^{(k-1)}))dB_s \right|^2 \right] \\
    &\leq 2E \left[ \left| \int_0^t (m(X_s^{(k)}) - m(X_s^{(k-1)}))ds \right|^2 \right] + E\left[ \left| \int_0^t (\sigma(X_s^{(k)}) - \sigma(X_s^{(k-1)}))dB_s \right|^2 \right] \\
    &\leq E \left[ \left(\int_0^t \beta^2 |X_s^{(k)} - X_s^{(k-1)}|\right)^2 \right] + \int_0^t E[\sigma(X_s^{(k)}) - \sigma(X_s^{(k-1)})]^2dB_s \\
    &\leq \beta^2 E[t \int_0^t |X_s^{(k)} - X_s^{(k-1)}|^2ds] + \beta^2 \int_0^t E[|X_s^{(k)} - X_s^{(k-1)}|^2]ds \\
    &\leq Ct \int_0^t E[|X_s^{(k)} - X_s^{(k-1)}|^2]
\end{align*}
$$
and specifically, one can now show
$$
    E[|X_t^{k+1} - X_t^{(k)}|] \leq \frac{\lambda^{k+1}t^{k+1}}{(k+1)!}
$$
which vanishes as $k \to \infty$.

_Def_ Let $X_t$ be a diffusion as above. The **generator** of $X_t$ is
$$
    Lf(x) = \lim_{t \to 0} \frac{E^x[f(X_t)] - E^x[f(X_0)]}{t}.
$$
Let $f$ be $C^2$; then
$$
\begin{align*}
    f(X_t) - f(X_0) &= \int_0^t f'(X_s)dX_s + \frac{1}{2} \int_0^t f''(X_s)d \left\langle X \right\rangle_s \\
    &= \int_0^t f'(X_s)\sigma(X_s)dB_s + \int_0^t (f'(X_s)m(X_s) + \frac{1}{2}f''(X_s)\sigma^2(X_s))ds.
\end{align*}
$$
Taking the expectation, under suitable regularity conditions (we will take $\sigma, m$ bounded), we get that
$$
    E[f(X_t)] - E[f(X_0)] = 0 + E[t(f'(X_0)m(X_0) + \frac{1}{2}f''(X_0)\sigma^2(X_0) + o(t^2)]
$$
and
$$
    Lf(x) = m(x)f'(x) + \frac{\sigma^2(x)}{2} f''(x).
$$

**Prop**: Let
$$
    dX_t = R_tdt + A_tdB_s
$$
and
$$
    dY_t = S_tdt + C_tdB_s.
$$
Then if we define the covariation as
$$
    \left\langle X, Y \right\rangle_t = \lim_{n \to \infty} \sum_{j < tn}(X_{j/n} - X_{(j-1)/n})(Y_{j/n} - Y_{(j-1)/n}) = \int_0^t A_sC_sds
$$
we have
$$
    dX_tY_t = X_tdY_t + Y_tdX_t + d\left\langle X, Y\right\rangle_t.
$$

Let $B_t = (B_t^1, \dots, B_t^d)$ be a Brownian motion. Then,
$$
    \left\langle B^j, B^k \right\rangle_t = 0
$$
for $k \neq j$. You can then write multidimensional stochastic integrals as stochastic integrals in each of the different dimensions; that is,
$$
    dX_t^{i} = R_t^{i}dt + \sum_{j=1}^d A_t^{(i,j)} dB_t^j.
$$
The covariation is then
$$
    \left\langle X^{(j)}, X^{(k)} \right\rangle_t = \sum_{i=1}^d A^{(j, i)}_tA^{(k, i)}_t.
$$

**Theorem (Multivariate Itö)**: Suppose $f(t, x)$ is a map from $\R^{n + 1} \to \R$, and is $C^1$ in $t$ and $C^2$ in $x$. Then,
$$
    X_t = (X_1^{(1)}, \dots, X_t^{(n)})
$$
satisfies 
$$
    df(t, X_t) = \partial_t f(t, X_t)dt + \sum_{j=1}^{n} \partial_{x_j}f(t, X_t)dX_t^{(j)} + \frac{1}{2}\sum_{j=1}^n\sum_{k=1}^n \partial_{x_j x_k}f(t, X_t)d \left\langle X^{(j)}, X^{(k)} \right\rangle_t.
$$

If $Z_t = \int_0^t A_sdB_s$, then $Z_t$ is not necessarily a martingale, but it is a local martingale.

_Def_: A process $Z_t$ is a **continuous local martingale** on $[0, \tau)$ (where $\tau$ could be $\infty$) if there is a sequence of stopping times $\tau_1 \leq \tau_2 \leq \cdots$ such that a.s. $\lim_{n \to \infty} \tau_n = \tau$ and for each $n$, $Z_{t \land \tau_n}$ is a continuous martingale.

**Prop**: Stochastic integrals are local martingales.

_Proof_: Take hitting times as stopping times.

### Feynman-Kac

Let $X_t$ be a geometric Brownian motion
$$
    dX_t = mX_tdt + \sigma X_t dB_t.
$$

Suppose we have a European option, such that at some fixed time $T > 0$ and fixed price $S$, we can exercise the option to buy some asset at $T$ for $S$. Then the value of an option at time $T$ is $F(X_t)$, where
$$
    F(x) = (x - s)_+ = \max \{x - s, 0\}.
$$

Let $\phi(t, x)$ be the value of this option at time $t < T$; that is,
$$
    \phi(t, x) = E[e^{-r(T - t)}F(X_T) \mid X_t = x].
$$
What PDE does $\phi(t, x)$ satisfy?

More generally, let
$$
\begin{gather*}
    dX_t = m(t, X_t)dt + \sigma(t, X_t)dB_t \\
    dR_t = r(t, X_t)R_tdt \\
    R_t = R_0 \exp\left( \int_0^t r(s, X_s)ds \right) \\
    \phi(t, x) = E \left[ \exp\left(-\int_{t}^Tr(s, X_s)ds\right) F(X_T) \mid X_t = x \right]
\end{gather*}
$$
Now define
$$
    M_t = E[R_T^{-1}F(X_T) \mid \mathcal F_t] = R_t^{-1}E \left[ \exp \left(-\int_t^T r(s,X_s)ds\right) F(X_T) \mid \mathcal F_t\right]
$$
so
$$
    M_t = R_t^{-1}\phi(t, X_t).
$$
is a martingale.

Assume for now that $\phi$ is sufficiently regular to apply Ito; then
$$
\begin{align*}
    d\phi(t, X_t) &= \partial_t \phi dt + \partial_x \phi dX_t + \frac{1}{2}\partial_{xx} \phi d \left\langle X \right\rangle_t \\
    &= \partial_t \phi dt + \partial_x \phi (mdt + \sigma dB_t) + \frac{1}{2}\partial_{xx} \phi \sigma^2 dt
\end{align*}
$$
and if you compute you get the following.

**Theorem**: Let $X_t$ be a geometric Brownian motion as above. Then, if $\phi$ is as above and is $C^1$ in $t$ and $C^2$ in $x$, then $\phi$ satisfies the PDE
$$
    \partial_t \phi(t, x) = -m(t,x) \partial_x \phi(t,x) - \frac{1}{2}\sigma(t,x)^2 \partial_{xx} \phi(t,x) + r(t,x) \phi(t,x)
$$
with terminal condition $\phi(T,x) = F(x)$.

Suppose that we have some SDE
$$
    dX_t = m(X_t)dt + \sigma(X_t)dB_t
$$
with $m, \sigma$ Lipschitz. Then the generator is
$$
    Lf(x) = m(x)f'(x) + \frac{\sigma^2(x)}{2}f''(x)
$$
and, if we have initial condition $F(x)$, then
$$
    u(t, x) = E^x[F(X_t)]
$$
and
$$
    \partial_t u(t,x) = L_x u(t,x), \ \ u(0, x) = F(x).
$$
On the other hand, if we have a terminal condition at $t = T$, we have
$$
    \partial_t\phi(t, x) = L_x(x, t), \ \ \phi(T, x) = F(X)
$$
and
$$
    \phi(t, x) = u(T - t, x), \ \  \partial_t \phi(t, x) = -L_x \phi(t, x).
$$

_Example_: Suppose we have a geometric Brownian motion with $m(x) = mx, \sigma(x) = \sigma x$, alongside some interest rate $r$. Then,
$$
    \phi(t, x) = E[e^{-r(T - t)}F(X_T)]
$$
and
$$
    \partial_t \phi(t, x) = r\phi(t, x) - mx\phi'(t,x) - \frac{\sigma^2}{2}x^2 \phi''(t, x).
$$

## Integrals Against Continuous Martingales

_Def_: If we have $f:[0, 1] \to \R$, its **variation** is
$$
    V_f = \sup_{0 = t_0 < \dots < t_n = 1} \sum_{j=1}^n |f(t_j) - f(t_{j-1})|.
$$
We say that $f$ has **bounded variation** if $V_f < \infty$.

**Prop**: Let $M_t$ be a continuous martingale with respect to the filtration $\{ \mathcal F_t \}$ and $M_0  = 0$. If $M_t$ has paths of bounded variation, then $M_t = 0$ for all $t$ with probability 1.

_Proof_: We will show that $E[M_1^2] = 0$. Then, in the case $V_M(1) \leq K < \infty$, 
$$
    E[M_1^2] = E \left[ \sum_{j=1}^n (M_{j/n} - M_{(j-1)/n})^2 \right].
$$
Bound the inner sum by $\delta_n \sum_{j=1}^n |M_{j/n} - M_{(j-1)/n}| \leq \delta_n K$ where $\delta_n$ is the maximal increment. Since $M_t$ is continuous, $\delta_n K \to 0$.
If we do not have that uniform bound, set $\tau_K = \min \{ t:  V_M(t) = K \}$. Then we have $E[M_{1 \land \tau_K}^2] = 0$ and we conclude by monotone convergence.

_Def_: The quadratic variation $\left\langle M \right\rangle_t$ is the unique increasing process such that $M_t^2 - \left\langle M \right\rangle_t$ is a martingale.

**Theorem**: If $M_t$ is a continuous martingale with respect to $\mathcal F_t$ with quadratic variation $\left\langle M \right\rangle_t$, then $M_t$ is a standard Brownian motion ith respect to $\mathcal F_t$.

_Proof_: Certainly $M_0 = 0$ and we have continuous paths. We only need to show that $E[e^{i\lambda (M_t - M_s)} \mid \mathcal F_s] = e^{-\frac{\lambda^2(t-s)}{2}}$. Since the definitions of the Ito integral and Ito's formula go through with only the quadratic variation assumption, we apply the Ito formula to $f(x) = e^{i\lambda x}$, so
$$
    f(M_t) - f(M_0) = \int_0^t f'(M_s)dB_s + \frac{1}{2} \int_0^t f''(M_s)ds
$$
and in expectation
$$
    E[f(M_t) - f(M_0)] = -\frac{\lambda^2}{2} \int_0^t E[f(M_s)]ds
$$
since $f'' = -\lambda^2f$. Then just solve the differential equation.

Take a standard Brownian motion $B_s, 0 \leq s \leq 1$; the Wiener measure $\mathcal W$ is a measure on the Borel $\sigma$-algebra of $C[0,1]$, the space of continuous functions with the $\sup$ norm. Let $B(f, \epsilon) = \{ g \mid |f - g| < \epsilon\}$. Then,
$$
    \mathcal W[B(f, \epsilon)] = P(|f(t) - B_t| < \epsilon).
$$
Put $\mathcal W_{m, \sigma^2}$ if we have a drift/variance term.

**Prop**: $\mathcal W_{0, 1} \perp \mathcal W_{0, \sigma^2}$ for $\sigma^2 \neq 1$. On the other hand, $\mathcal W_{0, 1}$  is mutually absolutely continuous with $\mathcal W_{m, 1}$.

_Proof_: The first part is easy. Look at $A$, the functions with quadratic variation 1; then $\mathcal W_{0, 1}(A) = 1$, $\mathcal W_{0, \sigma^2}(A) = 0$. In the second case, we go to a weak version of the Girsanov theorem.

**Theorem**: $\mathcal W$ and $\mathcal W_{m, 1}$ are mutually absolutely continuous and
$$
    \frac{d\mathcal W_{m,1}}{d \mathcal W} = \exp\left(mB_1(\omega) - \frac{1}{2}m^2\right).
$$

_Proof_: Let $B_t$ be a standard Brownian motion, and define $M_t = \exp(mB_t - \frac{m^2t}{2})$ so that $dM_t = mM_t dB_t$. On $\mathcal F_t$, define a probability measure $Q_t$ such that
$$
    Q_t[V] = E[1_V M_t]
$$
for all $V \in \mathcal F_t$; by conditioning on $\mathcal F_s$, we can see that if $s < t$ and $V \in \mathcal F_s$, then $Q_s(V) = Q_t(V)$. It is clear that the Radon-Nikodym derivatives are just given by $M_t$ and $1/M_t$.

Now define $Q$ a probability measure on $\mathcal F_\infty$ such that if $\mathcal F_{t}$, $Q(V) = E[1_V M_t]$, which is well defined because of the above. Moreover, if $X$ is $\mathcal F_t$ measurable, this gives that
$$
    E_Q[X] = E[X M_t].
$$
Now we claim that $B_t$ is a Brownian motion under $Q$ with variance parameter 1 and drift $m$. Certainly it holds that $B_0 = 0$ and $t \mapsto B_t$ is continuous with probability 1 under $Q$ and $P$ (since they are mutually absolutely continuous), and if $s, t > 0$ then the conditional distribution of $B_{t + s} - B_s$ given $\mathcal F_s$ is $N(mt, t)$. In fact, all we need to show then is that
$$
    E_Q[\exp(\lambda(B_{t+s} - B_s)) \mid \mathcal F_s] = \exp\left(\lambda mt + \frac{\lambda^2 m^2 t }{2}\right).
$$
This essentially boils down to the definition of the conditional expectation: write it down and unwrap in terms of $P$-expectations and win. 

In the above omitted computation, we essentially do a simple version of the following theorem.

**Theorem**: Let $B_t$ be a Brownian motion on a probability space $(\Omega, \mathcal F, P)$. Suppose we have a nonnegative martingale $M_t, M_0 = 1$ satisfying 
$$
    dM_t = A_tM_t dB_t.
$$
Then, $M_t = M_0e^{Y_t}$ where
$$
    Y_t = \int_0^t A_sdB_s - \frac{1}{2} \int_0^t A_s^2 ds.
$$
This follows by a simple application of Ito.
Define a measure $Q_t$ on $\mathcal F_t$, given by
$$
    Q_t[V] = E[M_t 1_V]
$$
so that $\frac{dQ_t}{dP} = M_t$ and if $s < t$ and $V \in \mathcal F_s$, $Q_s[V] = Q_t[V]$. Then, we define $Q$ as before. Let $X_t = B_t - \int_0^t A_sds$; then $X_t$ is a standard Brownian motion on $(\Omega, \mathcal F, Q)$.

_Proof_: We just do a heuristic argument. As an approximation,
$$
    B_{t + \Delta t} - B_t = \begin{cases}
        \sqrt{\Delta t} & \text{probability } 1/2 \\
        -\sqrt{\Delta t} & \text{probability } 1/2 \\
    \end{cases}
$$
and similarly $M_{t + \Delta t} = M_t(1 \pm A_t \sqrt{\Delta t})$ equiprobably as well. In the new measure, this is essentially tilting the probabilities by $1 \pm A_t \sqrt{\Delta t}$ and so $E_Q[B_{t + \Delta t} - B_t] = A_t \Delta t$.

Check that $X_t$ is a $Q$-martingale and conclude by noting that it has quadratic variation $t$.

## Conformal Invariance

Let $B_t = (B_t^1, B_t^2)$ be a two dimensional Brownian motion; identify it to $B_t^1 + iB_t^2$.

Let $f(B_t) = u(B_t) + iv(B_t)$. Then
$$
    du(B_t) = \nabla u \cdot dB_t + \frac{1}{2} (\Delta u) dt
$$
and
$$
    dv(B_t) = \nabla v \cdot dB_t + \frac{1}{2} (\Delta v) dt
$$
but if $f$ is holomorphic, then the Laplacians vanish and we use the Cauchy-Riemann equations to simplify to
$$
    du(B_t) = u_x(B_t)dB_t^{1} + u_y(B_t)dB_t^2
$$
and
$$
    dv(B_t) = v_x(B_t)dB_t^{1} + v_y(B_t)dB_t^2
$$
which have quadratic variations
$$
   d\langle u(B_t) \rangle_t = d\langle v(B_t) \rangle = (u_x^2 + u_y^2)dt = |f'(B_t)|^2  dt.
$$

Define $T(t) = \int_0^T |f'(B_s)|^2ds = T$. Then $Y_t = f(B_{T(t)})$ is a standard complex Brownian motion.

## Levy Processes

_Def_: A *Levy process* is a stochastic process that satisfies

- $X_0 = 0$ almost surely;
- stationary increments: $X_t - X_s \sim X_{t - s}$;
- independent increments;
- continuity in probability: $X_t \to 0$ in probability as $t \to 0$.

_Def_: A *compound Poisson process* is composed of a Poisson process $N_t$ with some rate $\lambda > 0$, alongside $Y_1, Y_2, \dots$ mutually independent of themselves and $N_t$ and identically distributed random variables with distribution $\mu^\sharp$. Then,
$$
    X_t = \sum_{i=1}^{N_t} Y_i
$$
is the compound Poisson process.

_Def_: A random variable $X$ has an *infinitely divisible* distribution if for each $n$ we may find $Y_1, \dots Y_n$ iid such that $X \sim \sum_{i=1}^n Y_i$.

**Prop**: Levy processes are infinitely divisible at any time.

Now look at the characteristic functions of a Levy process. In particular define the characteristic exponent $\Psi(s)$ where
$$
    \phi_{X_1}(s) = \exp(\Psi(s))
$$
and note that
$$
    \phi_{X_t}(s) = (\phi_{X_1}(s))^t.
$$

Let $X_t, Y_1, \dots$ be as in a compound Poisson process. Then, if we set
$$
    \phi(s) = \phi_{Y_j}(s) = \int_{-\infty}^\infty \exp(isx)\mu^\sharp (dx)
$$
we have
$$
    \phi_{X_1}(s) = \sum_{n=0}^\infty P(N_1 = n) E[e^{sX_1} \mid N_1 = n] = \sum_{n=0}^\infty \frac{\lambda^n}{n!} e^{-\lambda} (\phi(s))^n = \exp(\lambda(\phi(s) - 1))
$$
or equivalently,
$$
    \Psi(s) = \int_{-\infty}^\infty (e^{isx} - 1)\mu(dx)
$$
where $\mu = \lambda \mu^\sharp$. We call $\mu$ the Levy measure for $X_t$.

**Prop**: There are a few properties that we can say: as $t \to 0$,

- $P(N_t = 0) = 1 - \lambda t + o(t)$;
- $P(N_t \geq 2) = o(t)$;
- $P(X_t \in (a, b)) = \mu(a,b)t + o(t)$ for $0 \notin (a, b)$.

_Def_: We define the generator of a Levy process as
$$
    Lf(x) = \lim_{t \to 0}\frac{E^x[f(X_t) - f(x)]}{t} =  \int_{-\infty}^\infty (f(x + y) - f(x)) \mu(dy).
$$

**Prop**: If $X_t, Y_t$ are independent Levy processes, then $X_t + Y_t$ is a Levy process with characteristic exponent $\Psi_X + \Psi_Y$ and generator $L_x + L_y$.

_Def_: Let $X_t$ be a Levy process with measure $\mu$, and put
$$
    m = \int_{-\infty}^\infty x\mu(dx), \ \sigma^2 = \int_{-\infty}^\infty x^2 \mu(dx).
$$
Then $E[X_t] = mt$ and $\Var(X_t) = \sigma^2 t$.

_Def_: A function $f:[0, \infty) \to \R$ is called cadlad (continue a droite, limite a gauche) if the paths are right-continuous and for all $t$, $f(t-)$ exists.

_Def_: If $X_t$ is compound Poisson, then
$$
    \int_0^t A_s dX_s = \sum A_s (X_s - X_{s-})
$$
where the sum is over jump times.

The issue is that if $A_t$ is adapted, then $A_t = X_t$ shows that $\int A_s dX_s = 1$ at the first jump, and is not a martingale; in fact we need $A_t$ adapted to $\mathcal F_{t-}$ instead.
