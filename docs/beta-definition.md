# Boundary affection and the implemented beta policies

LAPSE computes conditional boundary exposure using exact finite-field relation
spaces. The default `separated` policy excludes consequences that arise only
by adding the public correctness relation. The exposure policy is part of the
experiment specification and is recorded with every estimate.

## Relation spaces

Write the original input and final output shares as a matrix $B$ with $2n$
rows over the independent source variables. For a leakage instance with row
matrix $L$, define

```math
H_L=\{a\in\mathbb F_q^{2n}:aB\in\operatorname{rowspan}(L)\}.
```

These are precisely the boundary coefficient vectors whose right-hand sides
are known from leakage. Public offsets do not change this coefficient space.
The compiler retains original input SSA sites and final output SSA sites;
intermediate versions remain independent physical leakage sites.

The analyzer requires a proper fresh refresh: equal, nonempty input and output
sharings; preserved share sum; and $\operatorname{rank}(B) = 2n - 1$. Thus the only leak-free
boundary relation is

```math
H_0=\operatorname{span}\{(1,\ldots,1,-1,\ldots,-1)\}.
```

This check rejects identity gadgets with two or more shares, incomplete outputs,
and circuits with additional boundary dependencies.
They need a model that accounts for their nontrivial baseline dependencies.

## Why a policy must be specified

The manuscript defines affected shares using support of boundary relations and
also excludes the public correctness relation. Simply discarding correctness
as a single row is insufficient: a known relation can be shifted by any
multiple of correctness, changing its support without changing its information.

For example, leaking only $X_1$ gives the relation $X_1 = c$. Adding correctness
also gives a relation containing the other boundary coordinates. Counting every
such shifted consequence would affect all shares. Selecting a representative
by row-reduction pivots instead would make affection depend on coordinate order.

The code offers two named policies. Numerical output always records the policy.

## Separated policy

Let $I_L$ be the subspace of $H_L$ supported on input coordinates only, and
$O_L$ the subspace supported on output coordinates only. Define the residual
coupling space

```math
C_L=H_L/(H_0+I_L+O_L).
```

The implemented exposure set is

```math
\operatorname{Aff}_{\mathrm{sep}}(L)=
\begin{cases}
\operatorname{supp}(I_L)\cup\operatorname{supp}(O_L), &\dim C_L=0,\\
\mathbf X\cup\mathbf Y, &\dim C_L>0.
\end{cases}
```

Here the support of a subspace is the union of the supports of its vectors;
computing it on any basis gives the same result. All quantities are independent
of elimination pivots and ordering of leakage sites.

$C_L$ removes both correctness and relations already explained separately on
the two sides. Given the native, a zero residual means that the affine boundary
constraints split into input-only and output-only constraints. Exposing their
supports is sufficient to provide their values. A nonzero residual triggers
exposure of every boundary share; this is a conservative choice for that
instance, rather than an attempt to find a minimal set of participating shares.

To see the conditional-distribution interpretation, let $h_X=(1,\ldots,1,0,\ldots,0)$.
Fixing the native changes the constraint space to $H'_L=H_L+\operatorname{span}(h_X)$.
Its one-sided subspaces are $I_L+\operatorname{span}(h_X)$ and
$O_L+\operatorname{span}(h_Y)$, where $h_Y=(0,\ldots,0,1,\ldots,1)$ and
$h_X-h_Y$ generates correctness. The conditional affine space is a Cartesian
product on its input and output coordinates exactly when these one-sided
subspaces span $H'_L$, equivalently when $C_L=0$. Uniform independent sources
induce a uniform distribution on this affine space. This establishes the
per-instance separation test, not independent domination across random leakage
instances.

This distinguishes genuine residual coupling from the sum of two already known
one-sided relations. For example, independently leaking $X_1$ and $Y_2$ affects
only those shares. Their sum is a cross-boundary relation, but it introduces no
additional dependence.

Important limitation: this policy's affected set need not be monotone in the
leakage set. In the two-share refresh $Y_1=X_1+R$, $Y_2=X_2-R$, observing $R$
leaves a coupled residual and exposes all four shares. Observing $R$ and $X_1$
also determines $Y_1$; the residual vanishes and the exposure set becomes
$\{X_1,Y_1\}$. Therefore the manuscript's monotone diagram-covering arguments and
its independence claims for inflated exposures cannot automatically be applied
to this policy. The computed marginal probability is not by itself a proof of
independent domination or of the composition theorem.

## Literal policy

The optional `literal` policy discards the baseline-only case but counts all
additional boundary relations and all their correctness shifts:

```math
\operatorname{Aff}_{\mathrm{lit}}(L)=
\begin{cases}
\emptyset,&H_L=H_0,\\
\mathbf X\cup\mathbf Y,&H_L\ne H_0.
\end{cases}
```

Every nonzero information class has a representative involving each coordinate,
because correctness has a nonzero coefficient in each coordinate. Thus even a
single leaked input share affects every boundary share under this policy. This
allows comparison with an interpretation that counts every correctness-shifted
support. The default policy for local exposure is `separated`.

## Conditional sampling and reported quantities

For a selected boundary SSA wire $U$, the implementation removes exactly that
wire from the sampling universe and samples every remaining allowed site at
rate $p$. This is precisely the conditional distribution given that $U$ does
not leak. Other SSA wires with the same vector, including copied boundaries,
remain sampleable. The selected wire must have been eligible under the chosen
leakage-role filter. At $p=1$ the conditioning event has probability zero, so
the API rejects the request.

The returned estimate is

```math
\beta_U(p)=\Pr[U\in\operatorname{Aff}(L)\mid U\notin L],\qquad
p_{\mathrm{effective}}=p+(1-p)\beta_U(p).
```

Results are per selected share. Symmetry is not assumed or inferred merely
from a circuit name. Exact enumeration, Monte Carlo, and importance sampling
use the same event. MC and likelihood-ratio IS remain valid for a deterministic
nonmonotone event. Confidence diagnostics describe sampling error; they do not
resolve the mathematical policy distinction. Raw IS estimates are not clipped.

As an analytic oracle, for the two-share refresh with exactly five leakage
sites $X_1,X_2,R,Y_1,Y_2$, conditioning on $X_1$ not leaking gives

```math
\beta_{\mathrm{sep}}(p)=p-p^2+2p^3-p^4,\qquad
\beta_{\mathrm{lit}}(p)=1-(1-p)^4.
```

The test suite checks these formulas over GF(2), GF(257), and GF(4), and checks
the exact relation space against an independently enumerated GF(2) oracle.

## Applying the estimates

Report the selected boundary SSA wire, leakage roles, field, circuit parameters,
and exposure policy. Any use of analytical bounds or composition results must
establish that their assumptions apply to that policy. In particular, the
manuscript's essential-bridge definition tests removal by correctness alone;
the separated policy tests residual coupling modulo correctness and both
one-sided spaces. These tests have different mathematical meanings.
