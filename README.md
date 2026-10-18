# An affine sunflower refuting the conjecture in ECCC TR26-007

Alekseev and Gaevoy's report *New Polynomial-Depth Res(⊕) Lower Bounds*
(ECCC TR26-007, January 2026) derives conditional lower bounds for
polynomial-depth Res(⊕) refutations from a combinatorial conjecture about
unions of affine subspaces (Conjecture 1.4, restated as Conjecture 4.2). This note gives an explicit family showing that
the conjecture is false as stated for every $q > 1$ and every $r > 0$, and
also for $q = 1$ under its stated uniformity in $r$.

## The conjecture

Let $r, q > 0$ be constants. Let $\Phi_1, \dots, \Phi_m$ be affine subspaces
of $\mathbb{F}_2^n$, each of codimension at most $(\log n)^q$, whose union has
at least $2^{n-1}$ points. For each $i$, let $\Phi_i' \subseteq \Phi_i$ be any
subset with $|\Phi_i'| \ge (1 - n^{-r})\,|\Phi_i|$. The conjecture asserts a
constant $c = c(q) > 0$, depending only on $q$, such that

$$
\Bigl|\bigcup_{j} \Phi_j'\Bigr| \ \ge\ \bigl(1 - n^{-rc}\bigr)\,\Bigl|\bigcup_{j} \Phi_j\Bigr| .
$$

The number $m$ of subspaces is unrestricted. (Conjecture 1.4 says
"collections of affine subspaces"; Conjecture 4.2 and its use in the proof of
Theorem 4.3 treat each $\Phi_j$ as one affine subspace, which is the reading
used here.) The report explains that any
naive counterexample has the structure of a sunflower, with every $\Phi_j'$
inside a core and the removed parts $\Phi_j \setminus \Phi_j'$ as petals, and
adds: "Intuitively, it seems like it is impossible to construct such a
sunflower using affine subspaces." The family below is such a sunflower.

## The construction

Fix integers $k, d \ge 1$ and $f \ge 0$, and let

$$
V = \bigl(\mathbb{F}_2^{d+1}\bigr)^k \times \mathbb{F}_2^f, \qquad n = k(d+1) + f .
$$

For every tuple $a = (a_1, \dots, a_k)$ of nonzero vectors
$a_j \in \mathbb{F}_2^{d+1}$, define the linear subspace

$$
S_a = \mathrm{span}(a_1) \times \cdots \times \mathrm{span}(a_k) \times \mathbb{F}_2^f .
$$

Let $L$ be the set of points of $V$ whose $k$ blocks are all nonzero, and
delete it from every member: $S_a' = S_a \setminus L$.

## Proof

1. **Codimension.** Each $S_a$ has dimension $k + f$, so codimension $kd$.
2. **Coverage.** The members cover $V$. Given $x = (x_1, \dots, x_k, w)$, take
   $a_j = x_j$ when $x_j \ne 0$ and any nonzero $a_j$ otherwise.
3. **Local loss.** A point of $S_a$ lies in $L$ only if every block equals
   $a_j$, so $S_a \cap L = \lbrace a \rbrace \times \mathbb{F}_2^f$. Each
   member loses exactly a $2^{-k}$ fraction of its points.
4. **Global loss.** The petals $\lbrace a \rbrace \times \mathbb{F}_2^f$ are
   disjoint and together make up $L$, while every point outside $L$ survives
   in some member. Hence the points lost from the union are exactly $L$, a
   fraction
   $$\frac{|L|}{|V|} = \bigl(1 - 2^{-(d+1)}\bigr)^k .$$
5. **Parameters for $q > 1$ and $r > 0$.** Take
   $k = \lceil r \log_2 n \rceil$, $d = \lfloor (\log n)^q / (2k) \rfloor$ and
   $f = n - k(d+1)$. The codimension $kd$ is at most $(\log n)^q$ and the local
   loss $2^{-k}$ is at most $n^{-r}$. Since $d = \Theta((\log n)^{q-1})$ grows
   faster than $\log_2 k = \Theta(\log \log n)$,
   $$\bigl(1 - 2^{-(d+1)}\bigr)^k \ \ge\ 1 - k\,2^{-(d+1)} \ \longrightarrow\ 1,$$
   whereas $n^{-rc} \to 0$ for every $c > 0$. The conclusion fails for all
   large $n$, even if $c$ may depend on both $q$ and $r$.
6. **The case $q = 1$.** Fix a candidate $c(1) > 0$ and choose a constant $d$
   with $\beta_d = -\log_2\bigl(1 - 2^{-(d+1)}\bigr) < c(1)$. For a small enough
   constant $r > 0$, the choice $k = \lceil r \log_2 n \rceil$ keeps the
   codimension $kd$ below $\log n$ and the local loss at most $n^{-r}$, while
   the global loss is at least $\bigl(1 - 2^{-(d+1)}\bigr)\, n^{-r\beta_d}$,
   which exceeds $n^{-r\,c(1)}$ for large $n$. So no $c$ depending only on
   $q = 1$ works for every $r > 0$. The construction does not decide
   $0 < q < 1$, or a $q = 1$ variant in which $c$ may depend on $r$.

## Remarks

- **The family is necessarily large.** Here
  $m = (2^{d+1} - 1)^k = 2^{\Theta((\log n)^q)}$. Some largeness is
  unavoidable: by the union bound the global loss is at most $m\,n^{-r}$, so
  the conclusion holds trivially whenever $m \le n^{r(1-c)}$.
- **Bounding $m$ does not rescue the statement for $q > 2$.** If one also
  requires $m \le 2^{(\log n)^{q/2}}$, take $k = \lceil r \log_2 n \rceil$ and
  $d + 1 = \lfloor (\log n)^{q/2} / k \rfloor$, which tends to infinity when
  $q > 2$; the same argument applies. For $q < 2$ the construction does not
  fit under that bound, and at $q = 2$ it fits only for small $r$, where it
  again refutes only the uniformity of $c$ in $r$. We have not analyzed which
  restriction on $m$, if any,
  would suffice for the application in Theorem 4.3, where the subspaces are
  the pieces $\Phi_j^\rho$ rather than the proof lines.
- **What is and is not affected.** Theorem 4.3 (Theorem 1.1) is an
  implication from the conjecture with $q > 1$ and $r = k - 2$. For $k > 2$
  the family refutes that premise, so the conjecture as stated cannot
  instantiate the theorem. The implication itself, the report's unconditional
  results, and the intended lower bounds are not refuted. A recent preprint
  claims an unconditional exponential lower bound for the bit pigeonhole
  principle in DAG-like Res(⊕) [2].

## Verification

[`verify.py`](verify.py) uses only the Python 3 standard library. It checks
every identity above by enumerating all of $\mathbb{F}_2^n$ for ten small
parameter choices up to $n = 12$. It then evaluates the parameter choices for
$q \in \lbrace 1.5, 2, 3 \rbrace$, for $q = 1$, and for the bounded-$m$
variant, measuring the codimension budget with natural logarithms (base-2
logarithms only enlarge it). For each instance it prints the largest exponent
$c$ that the instance still allows; these tend to $0$.

```sh
python3 verify.py
```

## Provenance

The construction was found by an autonomous AI research agent (Codex) working
in Prespatial Labs' research program. It was then re-derived by hand and
checked with `verify.py`, also with AI assistance (Claude). Corrections are
welcome as GitHub issues.

## References

1. Y. Alekseev and N. Gaevoy. *New Polynomial-Depth Res(⊕) Lower Bounds.*
   ECCC TR26-007, 2026. <https://eccc.weizmann.ac.il/report/2026/007/>
2. K. Braun. *An exponential lower bound for the bit pigeonhole principle in
   resolution over parities.* arXiv:2609.23015, 2026.
   <https://arxiv.org/abs/2609.23015>

## License

The text is licensed under CC BY 4.0 and the code under the MIT License; see
[`LICENSE`](LICENSE).
