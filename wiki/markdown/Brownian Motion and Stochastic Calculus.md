#  Brownian Motion and Stochastic Calculus 
##  UChicago STAT 38510, Autumn 2023 

\newcommand{\N}{\mathbb{N}}
\newcommand{\R}{\mathbb{R}}
\newcommand{\Q}{\mathbb{Q}}
\newcommand{\Z}{\mathbb{Z}}
\newcommand{\C}{\mathbb{C}}

\DeclareMathOperator{\Jac}{Jac}
\DeclareMathOperator{\Ker}{Ker}
\DeclareMathOperator{\mesh}{mesh}
\DeclareMathOperator{\Var}{Var}

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

###  Brownian Motion 

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

####  Construction 

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

####  Properties of Brownian Motion 

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


####  Filtrations 

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

####  The Markov Property of Brownian Motion 

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

####  Martingales 

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

####  Quadratic Variation 

Let $B_t$ be a standard Brownian motion; a partition $\Pi$ of $[0, 1]$ is a sequence $0 = t_0 < t_1 < \dots < t_k = 1$, and the mesh of the partition is just
$$
  \| \Pi \| = \max_{i = 1, \dots, k} t_i - t_{i - 1}.
$$
Now take a sequence of partitions $\Pi_n$, e.g. $0 = t_{0, n} < \dots < t_{k_n, n}$, and define the quantity
$$
  Q(t, \Pi) = \sum_{t_{j} <= t} (B_{t_j} - B_{t_{j-1}})^2
$$
and $Q_n(t) = Q(t, \Pi_n)$ and $Q_n = Q_n(1)$. 

**Theorem**: If $\| \Pi_n \| \to 0$, then $Q_n \to 1$ in probability. Futhermore, if $\sum_{n=1}^\infty \| \Pi_n \| < \infty$, then almost surely $\lim_{n \to \infty} Q_n = 1$.

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

####  Law of the Iterated Logarithm 

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

####  Zero Sets of Brownian Motion 

_Def_: Set $B_t$ a standard Brownian motion, $Z = \{ t \mid B_t = 0 \}$, and $Z_t = Z \cap [0, t]$. Then, $t \in Z$ is right-isolated if $t \in Z$, and $\exists \epsilon > 0$ such that $(t, t + \epsilon) \cap Z = \emptyset$; similar for left-isolated. A point which is both left and right isolated is just isolated.

With probability $1$, $0$ is not right-isolated; further, $Z_1$ is homeomorphic to the Cantor set.

With probability 1, the sets of left and right isolated points are countable, and there are no isolated points.

We can show that if $q \in \Q_{\geq 0}$, $P(B_q  = 0) = 0$, and $\tau_q = \min \{ t \geq q \mid B_t = 0\}$ is a stopping time; further, the left-isolated points are just $\{ \tau_q \mid q \in \Q_{\geq 0}\}$. And by strong Markov property, no $\tau_q$ is a right-isolated point, so there are no isolated points.

Set $\sigma_q = \max_t \{ t < q \mid B_t = 0 \}$ (this is well-defined, but not a stopping time). Then associate every $q$ to an interval $(\sigma_q, \tau_q)$; then $Z = [0, \infty) \setminus \bigcup_q (\sigma_q, \tau_q)$ and has Lesbegue measure 0. To see this, we just interchange integrals:
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

####  Local Time 

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

There are more facts: $L_t$ is continuous in $t$ and nondecreasing, $L_t - L_s > 0 \iff (s, t) \cap Z \neq \emptyset$, and $L_t$ is weakly $\frac{1}{2}$-Hölder continuous.

**Theorem (Scaling Rule)**: $L_t$ has the same distribution as $t^{1/2} L_1$. Further, $M_t = \max_{ 0 \leq s \leq t}B_s$ has the same distribution as $L_t$.

###  Brownian Motion in Several Dimensions 

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

####  Harmonic Functions in $\R^d$ 

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

####  Hitting Probabilities for Brownian Motion 

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


####  Recurrence and Transcience of Brownian Motion 

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

####  The Dirchlet Problem 

Take a bounded domain $D \subset \R^d$, and a continous function $F: \partial D \to \R$; the Dirichlet problem is to find the unique continuous $f: \overline D \to \R$ that agrees with $F$ on $\partial D$ and is harmonic on $D$.

In fact, uniqueness follows from the maximum principle, which says that the maximum of $f$ is attained on the boundary (think about the mean value principle). Then subtract two solutions and see that it is 0.

Let $T = \min\{t \mid B_t \in \partial D \}$, and $f(x) = E^x[F(B_T)]$. This satisfies the mean value principle, and in fact continuous (which is kinda hard) but it is obviously locally integrable, so $f$ here is harmonic.

In general, such a harmonic function is not necessarily continuous on the boundary: take the example of the punctured unit disk, with $F(x) = 1$ for $|x| = 1$ and $F(0) = 0$. And so everything is fine except at the origin, since $f(x) = P^x(|B_T| = 1) = 1$.

_Def_: If $x \in \partial D$, let $\sigma = \inf \{ t > 0 \mid B_t \in \partial D\}$; $x$ is regular if $P^x(\sigma = 0) = 1$.

**Prop**: $f$ as defined above is continuous at every regular boundary point.

Therefore the Dirichlet problem has a solution for every continuous $F$ if and only if every opint on $\partial D$ is regular.
