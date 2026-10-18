#!/usr/bin/env python3
"""Check the product-of-lines counterexample to Conjecture 1.4/4.2 of ECCC TR26-007.

Part 1 verifies every identity used in the proof by enumerating all of F_2^n
for small parameters.  Part 2 evaluates the asymptotic parameter choices with
exact integer arithmetic and reports the largest exponent c that each instance
still allows; it tends to 0, so no fixed c > 0 works.  Part 3 does the same for
q = 1 with a fixed c, and Part 4 for families restricted to
m <= 2^((ln n)^(q/2)) when q > 2.

Standard library only.  Exit status 0 means every check passed.
"""

from fractions import Fraction
from itertools import product
import math

LN2 = math.log(2)


def check_small(k: int, d: int, f: int) -> None:
    """Enumerate V = (F_2^(d+1))^k x F_2^f and check the identities exactly."""
    n = k * (d + 1) + f
    zero = (0,) * (d + 1)
    blocks = list(product((0, 1), repeat=d + 1))
    nonzero = [b for b in blocks if any(b)]
    free = list(product((0, 1), repeat=f))

    vectors = {(bs, w) for bs in product(blocks, repeat=k) for w in free}
    assert len(vectors) == 2**n

    def in_lost_set(x) -> bool:
        return all(any(b) for b in x[0])

    lost = {x for x in vectors if in_lost_set(x)}
    union, union_after = set(), set()
    members = list(product(nonzero, repeat=k))
    assert len(members) == (2 ** (d + 1) - 1) ** k

    for a in members:
        # S_a = span(a_1) x ... x span(a_k) x F_2^f: a linear subspace.
        member = {
            (bs, w)
            for bs in product(*[(zero, a_j) for a_j in a])
            for w in free
        }
        assert len(member) == 2 ** (k + f)  # dimension k+f, codimension k*d
        removed = {x for x in member if in_lost_set(x)}
        assert removed == {(a, w) for w in free}  # exactly 2^-k of S_a
        union |= member
        union_after |= member - removed

    assert union == vectors  # the family covers all of V
    assert vectors - union_after == lost  # the actual loss is exactly L
    density = Fraction(len(lost), 2**n)
    assert density == (1 - Fraction(1, 2 ** (d + 1))) ** k
    print(
        f"  k={k} d={d} f={f}: n={n}, m={len(members)}, codim={k * d}, "
        f"local loss=2^-{k}, global loss={density} = {float(density):.4f}"
    )


def log2_global_loss(k: int, d: int) -> float:
    """log2 of the lost fraction (1 - 2^-(d+1))^k, accurate for large d."""
    return k * math.log1p(-(2.0 ** -(d + 1))) / LN2


def log2_largest_c(k: int, d: int, r: float, t: int) -> float:
    """log2 of the largest c with (1 - 2^-(d+1))^k <= n^(-r c), n = 2^t.

    That c equals -log2((1 - 2^-(d+1))^k) / (r t).  For d + 1 > 60 the
    expansion -log(1 - x) = x (1 + x/2 + ...) with x = 2^-(d+1) gives its
    logarithm to within 2^-60, avoiding floating-point underflow.
    """
    if d + 1 <= 60:
        return math.log2(-log2_global_loss(k, d)) - math.log2(r * t)
    return math.log2(k) - (d + 1) - math.log2(LN2) - math.log2(r * t)


def show_power(log2_value: float) -> str:
    if log2_value > -1000:
        return f"{2.0 ** log2_value:.3g}"
    return f"2^{log2_value:.0f}"


def asymptotic_rows(q: float, r: float, exponents) -> list[float]:
    """Instances for fixed q > 1 and r > 0 with n = 2^t (natural-log budget)."""
    allowed = []
    for t in exponents:
        budget = (t * LN2) ** q  # (ln n)^q; base-2 logs only enlarge it
        k = math.ceil(r * t)  # local loss 2^-k <= n^-r
        d = math.floor(budget / (2 * k))
        assert d >= 1 and k * d <= budget
        assert 2**t >= k * (d + 1)  # room for the free coordinates
        log2_c = log2_largest_c(k, d, r, t)
        allowed.append(log2_c)
        print(
            f"  q={q} r={r} n=2^{t}: k={k} d={d} codim={k * d} "
            f"<= {budget:.0f}; global loss={2 ** log2_global_loss(k, d):.4g}; "
            f"largest c allowed={show_power(log2_c)}"
        )
    assert all(b < a for a, b in zip(allowed, allowed[1:]))
    return allowed


def q_equals_one(c: float, exponents) -> None:
    """q = 1: a fixed c(1) fails for a suitably small r."""
    d = 1
    while -log2_global_loss(1, d) >= c:
        d += 1
    beta = -log2_global_loss(1, d)
    r = 0.9 * LN2 / d  # keeps d * ceil(r t) <= ln n for large t
    print(f"  c(1)={c}: choose d={d} (beta_d={beta:.4f} < c) and r={r:.4f}")
    for t in exponents:
        k = math.ceil(r * t)
        assert k * d <= t * LN2
        gap = log2_global_loss(k, d) + r * c * t  # log2(global / n^-rc)
        print(
            f"    n=2^{t}: codim={k * d} <= ln n={t * LN2:.0f}; "
            f"log2(global loss / n^(-rc)) = {gap:.2f}"
        )
    assert gap > 0


def bounded_family(q: float, r: float, exponents) -> None:
    """q > 2: the construction also respects m <= 2^((ln n)^(q/2))."""
    for t in exponents:
        k = math.ceil(r * t)
        d = math.floor((t * LN2) ** (q / 2) / k) - 1
        assert d >= 1
        log2_m = k * math.log2(2 ** (d + 1) - 1)
        assert log2_m <= (t * LN2) ** (q / 2) and k * d <= (t * LN2) ** q
        log2_c = log2_largest_c(k, d, r, t)
        print(
            f"  q={q} r={r} n=2^{t}: d={d}, log2 m={log2_m:.0f} <= "
            f"{(t * LN2) ** (q / 2):.0f}; largest c allowed={show_power(log2_c)}"
        )


def main() -> None:
    print("Part 1: exhaustive check of the identities")
    for k, d, f in [(1, 1, 0), (1, 2, 1), (2, 1, 0), (2, 1, 2), (2, 2, 0),
                    (3, 1, 0), (2, 3, 0), (3, 2, 1), (4, 2, 0), (3, 3, 0)]:
        check_small(k, d, f)

    print("Part 2: fixed q > 1 and r > 0; the largest admissible c tends to 0")
    for q, r in [(1.5, 1), (2, 1), (2, 3), (3, 1)]:
        asymptotic_rows(q, r, [40, 160, 640, 2560, 10240])

    print("Part 3: q = 1 with c depending only on q")
    q_equals_one(0.1, [100, 1000, 10000, 100000])

    print("Part 4: q > 2 with the family size bounded by 2^((ln n)^(q/2))")
    bounded_family(3, 1, [160, 640, 2560, 10240])

    print("All checks passed.")


if __name__ == "__main__":
    main()
