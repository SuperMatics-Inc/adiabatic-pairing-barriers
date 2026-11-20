#!/usr/bin/env python3
"""LINEAR7 replay: Holstein site-field Berry-phase determinant zero (Theorem Z, Lemmas 2-4),
gauge-twist power-law bound (Theorem W), static non-log-concavity (Lemma C),
and the cutoff-ceiling corollary (Corollary 10). Numbering follows RIGOR.md.

Default run exits 0 iff every exact / Arb-ball / interval check passes.
--mutant NAME runs the certificate with one guarded defect and must exit 1.
--all-mutants requires every guarded child to exit 1 and print "REJECTED <name>:".
Labels: CERTIFIED = exact rational or Arb ball; NUMERICAL = float sanity check.
"""
import argparse
import itertools
import math
import subprocess
import sys
from fractions import Fraction as Q
from pathlib import Path

import mpmath as mp
import numpy as np
from flint import arb, arb_mat, ctx

MUTANTS = (
    "smooth-rotation-symmetric-K",
    "three-piece-factor",
    "chebyshev-angle-2pi-over-N",
    "neumann-constant-2pi",
    "star-pattern-nonzero-mean",
    "star-pattern-order",
    "lattice-bracket-4.9",
    "lattice-history-mean",
    "twist-exponent-over-4",
    "shell-bond-count",
    "twist-identity-sign",
    "posterior-variance-quarter",
    "dynamic-variance-sixth",
    "isolated-star-bracket",
)


class Rejected(Exception):
    pass


def require(cond, label):
    if not cond:
        raise Rejected(label)


def A(q):
    """Exact rational -> Arb ball."""
    q = Q(q)
    return arb(q.numerator) / arb(q.denominator)


def certainly_pos(x):
    return bool(x > 0)


def certainly_neg(x):
    return bool(x < 0)


# ---------------------------------------------------------------------------
# 1. Two-level effective models (Lemma 1 constant pi, Lemma 2, Lemma 3, Proposition 5 value)
# ---------------------------------------------------------------------------
def mat2_exp(m):
    return arb_mat(m).exp()


def check_two_level(mutant):
    ctx.prec = 200
    pi = arb.pi()
    # (a) smooth rotation: K = [[A, -p],[p, -A]] satisfies K^2 = (A^2 - p^2) I  (exact polynomial identity)
    for a_, p_ in [(Q(1), Q(2)), (Q(3, 7), Q(5, 11)), (Q(-2), Q(9, 4))]:
        if mutant == "smooth-rotation-symmetric-K":
            K = [[a_, p_], [p_, -a_]]
        else:
            K = [[a_, -p_], [p_, -a_]]
        K2 = [[sum(K[i][k] * K[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
        require(K2 == [[a_ * a_ - p_ * p_, 0], [0, a_ * a_ - p_ * p_]], "smooth-rotation K^2=(A^2-p^2)I")
    # det(1+U_smooth) = 2 - 2cos(sqrt(pi^2 - (bA)^2)); Arb check against U = -exp(-beta K), beta=1
    for bA in (Q(1), Q(3), Q(31, 10), Q(32, 10), Q(4)):
        K = [[A(bA), -pi], [pi, -A(bA)]]
        U = mat2_exp([[-K[i][j] for j in range(2)] for i in range(2)]) * arb(-1)
        d = (arb_mat([[1, 0], [0, 1]]) + U).det()
        s = pi * pi - A(bA) * A(bA)
        formula = 2 - 2 * s.sqrt().cos() if certainly_pos(s) else 2 - 2 * (-s).sqrt().cosh()
        require((d - formula).contains(0) and (d - formula).rad() < arb("1e-40"), "smooth closed form")
        require(certainly_pos(d) if bA < Q(3141, 1000) else certainly_neg(d), "smooth sign change at pi")
    # (b) three-piece factorization: 2 - 2 T_3(c/2) = -(c-2)(c+1)^2 (cubic identity checked at 5 points)
    for c in (Q(0), Q(1), Q(2), Q(-3), Q(7, 5)):
        T3 = 4 * (c / 2) ** 3 - 3 * (c / 2)
        rhs = -(c - 2) * (c + 1) if mutant == "three-piece-factor" else -(c - 2) * (c + 1) ** 2
        require(2 - 2 * T3 == rhs, "three-piece factorization")
    # (c) N-piece closed form: tr U_N = -2 T_N(cos(pi/N) cosh(bA/N)), directions 2*pi*j/N in the doubled plane
    for N in (3, 4, 6):
        for bA in (Q(1), Q(3), Q(9, 2)):
            M = arb_mat([[1, 0], [0, 1]])
            s = A(bA) / N
            for j in range(N):
                ang = 2 * pi * j / N
                n = [[ang.cos(), ang.sin()], [ang.sin(), -ang.cos()]]
                M = mat2_exp([[-s * n[i][k] for k in range(2)] for i in range(2)]) * M
            ang0 = (2 * pi / N) if mutant == "chebyshev-angle-2pi-over-N" else (pi / N)
            cN = ang0.cos() * s.cosh()
            TN = (N * cN.acos()).cos() if bool(abs(cN) < 1) else (N * cN.acosh()).cosh()
            require((M.trace() + 2 * TN).contains(0), f"N-piece closed form N={N}")
    # (d) thresholds bA_N = N arccosh(sec(pi/N)) exceed pi and decrease
    prev = None
    for N in range(3, 65):
        th = N * (1 / (pi / N).cos()).acosh()
        require(certainly_pos(th - pi), "threshold > pi")
        if prev is not None:
            require(certainly_pos(prev - th), "threshold decreasing")
        prev = th
    th3 = 3 * arb(2).acosh()
    th6 = 6 * (2 / arb(3).sqrt()).acosh()
    kappa_star_eff = arb(3) / 2 * th6
    require(bool(abs(th3 - arb("3.95087369077445")) < arb("1e-13")), "3 arccosh 2")
    require(bool(abs(kappa_star_eff - arb("4.94375529900650")) < arb("1e-13")), "kappa_eff star")
    # (e) Neumann positive side: kappa_pos = pi must not exceed any proved zero location
    kappa_pos = 2 * pi if mutant == "neumann-constant-2pi" else pi
    require(not bool(kappa_pos > pi), "positive side <= smooth sharp zero at pi")
    require(certainly_neg(kappa_pos - 5), "positive side below certified lattice zero")
    return th3, th6, kappa_star_eff


# ---------------------------------------------------------------------------
# 2. Ideal star (K_{1,3}) zero-mode geometry: exact rationals (Lemma 4)
# ---------------------------------------------------------------------------
PAT6 = [(1, -1, -1), (1, 1, -1), (-1, 1, -1), (-1, 1, 1), (-1, -1, 1), (1, -1, 1)]


def star_pattern(mutant):
    pat = [list(p) for p in PAT6]
    if mutant == "star-pattern-nonzero-mean":
        pat[0] = [1, 1, 1]
    if mutant == "star-pattern-order":
        pat[1], pat[2] = pat[2], pat[1]
    return pat


def check_star(mutant):
    pat = star_pattern(mutant)
    # star adjacency: centre 0, leaves 1..3; zero modes = leaf vectors with zero sum
    P = [[(Q(1) if i == j else Q(0)) - Q(1, 3) for j in range(3)] for i in range(3)]
    # h_star w = 0 on Ran P
    for w in ([1, -1, 0], [1, 1, -2]):
        require(sum(w) == 0, "zero mode: centre row")  # (h w)_centre = -sum w; (h w)_leaf = -w_centre = 0
    # per-leaf time mean zero and sup = 1
    for leaf in range(3):
        require(sum(p[leaf] for p in pat) == 0, "star pattern zero mean")
    T = []
    for p in pat:
        PVP = [[sum(P[i][k] * p[k] * P[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
        tr = sum(PVP[i][i] for i in range(3))
        c = tr / 2  # rank-2 projector: identity part
        T.append([[PVP[i][j] - c * P[i][j] for j in range(3)] for i in range(3)])
    def tr2(X, Y):
        return sum(X[i][k] * Y[k][i] for i in range(3) for k in range(3))
    def comm(X, Y):
        return [[sum(X[i][k] * Y[k][j] - Y[i][k] * X[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    amp2 = [tr2(X, X) for X in T]
    require(all(x == Q(8, 9) for x in amp2), "traceless amplitude 2/3 (tr T^2 = 8/9)")
    cosines = [tr2(T[j], T[(j + 1) % 6]) / Q(8, 9) for j in range(6)]
    require(all(x == Q(1, 2) for x in cosines), "consecutive doubled angle 60 degrees")
    E = comm(T[0], T[1])
    for j in range(6):
        Cj = comm(T[j], T[(j + 1) % 6])
        require(Cj == E, "consistent orientation of rotation")
    require(any(x != 0 for row in E for x in row), "non-commuting pieces")
    return Q(2, 3)


# ---------------------------------------------------------------------------
# 3. Full 4x4 torus witness: Arb determinant signs (Theorem Z)
# ---------------------------------------------------------------------------
LT = 4
NS = LT * LT


def idx(x, y):
    return (x % LT) * LT + (y % LT)


CENTRE, LEAF_A, LEAF_B, LEAF_C = idx(0, 0), idx(0, 1), idx(3, 0), idx(0, 3)
STAR = (CENTRE, LEAF_A, LEAF_B, LEAF_C)
ADJ = [[0] * NS for _ in range(NS)]
for _x in range(LT):
    for _y in range(LT):
        for _dx, _dy in ((1, 0), (0, 1)):
            _i, _j = idx(_x, _y), idx(_x + _dx, _y + _dy)
            ADJ[_i][_j] = ADJ[_j][_i] = 1
L_FROZEN = Q(400)
W_STAR = Q(1, 200)


def lattice_history(kappa, beta, mutant):
    """Six pieces of length beta/6; leaves A,B,C get V*PAT6[j], V = kappa/beta; all other sites 0."""
    pat = [list(p) for p in PAT6]
    if mutant == "lattice-history-mean":
        pat[3] = [1, 1, 1]
    V = Q(kappa) / Q(beta)
    hist = []
    for p in pat:
        v = [Q(0)] * NS
        v[LEAF_A], v[LEAF_B], v[LEAF_C] = V * p[0], V * p[1], V * p[2]
        hist.append(v)
    return hist


def onebody_U(beta, hist, Wref):
    n = len(hist)
    M = arb_mat([[1 if i == j else 0 for j in range(NS)] for i in range(NS)])
    for v in hist:
        H = [[A(-ADJ[i][j]) for j in range(NS)] for i in range(NS)]
        for i in range(NS):
            H[i][i] = A(Wref[i] + v[i])
        M = (arb_mat(H) * (-A(Q(beta) / n))).exp() * M
    return M


def eye():
    return arb_mat([[1 if i == j else 0 for j in range(NS)] for i in range(NS)])


def check_lattice(mutant, verbose):
    Wref = [W_STAR if i in STAR else L_FROZEN for i in range(NS)]
    results = []
    lo_hi = {Q(20): (Q(3), Q(51, 10)), Q(40): (Q(3), Q(5)), Q(100): (Q(3), Q(5))}
    if mutant == "lattice-bracket-4.9":
        lo_hi[Q(100)] = (Q(3), Q(49, 10))
    for beta, (klo, khi) in lo_hi.items():
        ctx.prec = int(float(beta) * 5) + 256
        hist_lo = lattice_history(klo, beta, mutant)
        hist_hi = lattice_history(khi, beta, mutant)
        for hist in (hist_lo, hist_hi):
            for i in range(NS):
                require(sum(v[i] for v in hist) == 0, "lattice history time-mean zero")
        d0 = (eye() + onebody_U(beta, [[Q(0)] * NS], Wref)).det()
        dlo = (eye() + onebody_U(beta, hist_lo, Wref)).det()
        dhi = (eye() + onebody_U(beta, hist_hi, Wref)).det()
        require(certainly_pos(d0) and certainly_pos(dlo), f"det(1+U)>0 at v=0 and kappa={klo}, beta={beta}")
        require(certainly_neg(dhi), f"det(1+U)<0 at kappa={khi}, beta={beta}")
        results.append((beta, klo, khi, dlo, dhi))
    # robustness: deterministic rational perturbation of the reference (|dW| <= 1/1000), beta = 40
    beta = Q(40)
    ctx.prec = 456
    rng = np.random.default_rng(20261119)
    Wp = [w + Q(int(rng.integers(-1000, 1001)), 10 ** 6) for w in Wref]
    dpert = (eye() + onebody_U(beta, lattice_history(Q(5), beta, mutant), Wp)).det()
    require(certainly_neg(dpert), "perturbed reference still negative at kappa=5")
    # near-zero blow-up of the normalized conditional Green function (beta = 40): bisection in rationals
    lo, hi = Q(49, 10), Q(5)
    for _ in range(30):
        mid = (lo + hi) / 2
        dm = (eye() + onebody_U(beta, lattice_history(mid, beta, mutant), Wref)).det()
        if certainly_pos(dm):
            lo = mid
        elif certainly_neg(dm):
            hi = mid
        else:
            break
    U = onebody_U(beta, lattice_history(lo, beta, mutant), Wref)
    G = (eye() + U).inv()
    g_ab = G[LEAF_B, LEAF_A]
    require(bool(abs(g_ab) > 10 ** 9), "conditional |G_ba| > 1e9 near the zero")
    pair_corr = g_ab * g_ab
    require(bool(pair_corr > 10 ** 19), "conditional pair correlation > 1e19")
    return results, dpert, (lo, hi), g_ab


def check_isolated(mutant):
    """Isolated K_{1,3} star (the L_frozen -> infinity limit, RIGOR Lemma 6 and Theorem Z (Z2)), t = 1, W = 0.
    Certified sign change of det(1+U) along the same six-piece +-V leaf family."""
    adj = [[0, 1, 1, 1], [1, 0, 0, 0], [1, 0, 0, 0], [1, 0, 0, 0]]
    brackets = {Q(1): (Q(3), Q(12)), Q(10): (Q(3), Q(11, 2)), Q(20): (Q(3), Q(51, 10)),
                Q(40): (Q(3), Q(5)), Q(100): (Q(49, 10), Q(5)), Q(1000): (Q(49, 10), Q(5))}
    if mutant == "isolated-star-bracket":
        brackets[Q(1000)] = (Q(3), Q(49, 10))
    i4 = arb_mat([[1 if i == j else 0 for j in range(4)] for i in range(4)])
    out = []
    for beta, (klo, khi) in brackets.items():
        ctx.prec = int(float(beta) * 3) + 256
        dets = []
        for kappa in (Q(0), klo, khi):
            V = kappa / beta
            M = i4
            for p in PAT6:
                H = [[A(-adj[i][j]) for j in range(4)] for i in range(4)]
                for leaf in (1, 2, 3):
                    H[leaf][leaf] = A(V * p[leaf - 1])
                M = (arb_mat(H) * (-A(beta / 6))).exp() * M
            dets.append((i4 + M).det())
        require(certainly_pos(dets[0]) and certainly_pos(dets[1]),
                f"isolated star det(1+U)>0 at v=0 and kappa={klo}, beta={beta}")
        require(certainly_neg(dets[2]), f"isolated star det(1+U)<0 at kappa={khi}, beta={beta}")
        out.append((beta, klo, khi, dets[1], dets[2]))
    return out


# ---------------------------------------------------------------------------
# 4. Theorem W constants and a float sanity check of the twist identity
# ---------------------------------------------------------------------------
def torus_shell_bonds(side, r):
    def dist(p):
        dx = min(p[0] % side, side - p[0] % side)
        dy = min(p[1] % side, side - p[1] % side)
        return dx + dy
    counts = {}
    within = 0
    for x in range(side):
        for y in range(side):
            for dx, dy in ((1, 0), (0, 1)):
                k1, k2 = dist((x, y)), dist((x + dx, y + dy))
                lo = min(k1, k2)
                if lo < r:
                    if k1 == k2:
                        within += 1
                    else:
                        require(abs(k1 - k2) == 1, "bond changes distance by 1")
                        counts[lo] = counts.get(lo, 0) + 1
    return counts, within


def check_twist(mutant):
    ctx.prec = 200
    for r in range(1, 8):
        counts, within = torus_shell_bonds(2 * r + 2, r)
        require(within == 0, "no intra-shell bonds")
        for k in range(r):
            claim = 4 * (2 * k + 2) if mutant == "shell-bond-count" else 4 * (2 * k + 1)
            require(counts.get(k, 0) == claim, "shell bond count 4(2k+1)")
    # cosh s - 1 <= (s^2/2) cosh s : coefficientwise 2 (2m-2)! <= (2m)!
    for m in range(1, 60):
        require(2 * math.factorial(2 * m - 2) <= math.factorial(2 * m), "cosh coefficient inequality")
    # g(r) = sum_{j<r} 1/(2j+1) >= log(2r+1)/2
    g = arb(0)
    for r in range(1, 3001):
        g += arb(1) / (2 * r - 1)
        if r in (1, 2, 5, 10, 100, 1000, 3000):
            require(certainly_pos(g - arb(2 * r + 1).log() / 2), "harmonic-odd lower bound")
    # exponent: eta = 1/(4 max(bt,1/4)); actual exponent eta(1 - 2 bt eta cosh eta) >= claimed theta
    c1 = arb(1).cosh()
    for bt in (Q(1, 100), Q(1, 4), Q(1, 2), Q(1), Q(10), Q(1000)):
        m = max(bt, Q(1, 4))
        eta = 1 / (4 * A(m))
        actual = eta * (1 - 2 * A(bt) * eta * eta.cosh())
        denom = 4 if mutant == "twist-exponent-over-4" else 8
        theta = (2 - c1) / (denom * A(m))
        require(not bool(theta > actual) and certainly_pos(actual - theta + arb("1e-30")), "Theorem W exponent")
    theta_coeff = (2 - c1) / 8
    require(bool(abs(theta_coeff - arb("0.0571149206")) < arb("1e-9")), "theta coefficient value")
    # NUMERICAL: twist identity and bound on a 4-site spinful Hubbard ring (stand-in for any
    # interaction commuting with all n_x); x = 0, y = 2.
    nsite, t, Uh, beta, mu = 4, 1.0, -4.0, 1.5, -1.0
    nm = 2 * nsite
    dim = 2 ** nm
    def ann(mode):
        Cm = np.zeros((dim, dim))
        for s in range(dim):
            if s >> mode & 1:
                sign = (-1) ** bin(s & ((1 << mode) - 1)).count("1")
                Cm[s ^ (1 << mode), s] = sign
        return Cm
    c = [ann(m) for m in range(nm)]  # mode = 2*site + spin
    cd = [x.T for x in c]
    num = [cd[m] @ c[m] for m in range(nm)]
    def build(phi):
        H = np.zeros((dim, dim))
        for i in range(nsite):
            j = (i + 1) % nsite
            for s in range(2):
                a_, b_ = 2 * i + s, 2 * j + s
                H += -t * (np.exp(phi[i] - phi[j]) * cd[a_] @ c[b_] + np.exp(phi[j] - phi[i]) * cd[b_] @ c[a_])
            H += Uh * num[2 * i] @ num[2 * i + 1]
        for m in range(nm):
            H += -mu * num[m]
        return H
    def expm_sym(H):
        w, V = np.linalg.eigh((H + H.T) / 2)
        return (V * np.exp(-beta * w)) @ V.T
    from scipy.linalg import expm
    H0 = build([0.0] * nsite)
    rho = expm_sym(H0)
    Z = np.trace(rho)
    Dx = c[1] @ c[0]  # Delta_0 = c_{0 down} c_{0 up}
    Dy = c[5] @ c[4]
    B = Dx.T @ Dy
    lhs = np.trace(B @ rho) / Z
    eta = 0.3
    phi = [0.0, eta, 2 * eta, eta]
    Hphi = build(phi)
    sgn = 1.0 if mutant == "twist-identity-sign" else -1.0
    twisted = math.exp(sgn * 2 * (phi[2] - phi[0])) * np.trace(B @ expm(-beta * Hphi)) / Z
    require(abs(lhs - twisted) < 1e-10 * max(1.0, abs(lhs)), "twist identity (NUMERICAL)")
    best = min(math.exp(-4 * e + 2 * beta * t * 4 * (math.cosh(e) - 1)) for e in np.linspace(0.01, 3, 300))
    require(abs(lhs) <= best + 1e-12, "twist bound (NUMERICAL)")
    return theta_coeff, lhs, best


# ---------------------------------------------------------------------------
# 5. Lemma C: static marginal is not log-concave at t = 0, mu = -U
# ---------------------------------------------------------------------------
def check_static(mutant):
    # posterior weights at b = -U: P(n) exp(-beta U (n-1)^2 / 2) with P(n) ~ (1, 2e^{-beta U/2}, 1)
    w = [Q(1), Q(2), Q(1)]  # common factor e^{-beta U/2} removed exactly
    Zw = sum(w)
    mean = sum(n * w[n] for n in range(3)) / Zw
    var = sum(n * n * w[n] for n in range(3)) / Zw - mean ** 2
    var_claim = Q(1, 4) if mutant == "posterior-variance-quarter" else var
    require(var == Q(1, 2), "posterior variance 1/2")
    mp.mp.dps = 40
    for bU in (mp.mpf(4), mp.mpf(3), mp.mpf(1)):
        U = mp.mpf(1)
        beta = bU / U
        P = [1, 2 * mp.e ** (-bU / 2), 1]
        p = lambda b: sum(P[n] * mp.e ** (-beta * (b + U * n) ** 2 / (2 * U)) for n in range(3))
        hess = -mp.diff(lambda b: mp.log(p(b)), -U, 2)
        formula = (beta / U) - beta ** 2 * var_claim
        require(abs(hess - formula) < mp.mpf(10) ** -20, "Hessian identity at b=-U")
        require((hess < 0) == (bU > 2), "non-log-concave iff beta U > 2")
    return var


# ---------------------------------------------------------------------------
# 6. Cutoff-ceiling corollary: dynamic variance bracket and exclusion probability
# ---------------------------------------------------------------------------
def check_ceiling(mutant):
    ctx.prec = 200
    denom = 6 if mutant == "dynamic-variance-sixth" else 12
    for a in ("0.001", "0.1", "0.5", "1", "2", "4", "7"):
        aa = arb(a)
        x = aa / 2
        ratio = (x * (1 / x.tanh()) - 1) / (aa * aa / denom)  # sigma^2 / (U Omega^2 beta / denom)
        require(not bool(ratio > 1) and not bool(ratio < 1 - aa * aa / 60), "sigma^2 bracket")
    # exclusion probability per site, kappa = 5: P >= erfc(5 sqrt(12) sqrt(z) / sqrt(2 (1 - a^2/60)))
    probs = []
    for z in ("1/300", "1/1000", "1/30"):
        zz = arb(z.split("/")[0]) / arb(z.split("/")[1])
        arg = 5 * arb(12).sqrt() * zz.sqrt() / (2 * (1 - arb(1) / 60)).sqrt()
        probs.append((z, arg.erfc()))
    require(bool(probs[0][1] > arb("0.3")), "exclusion probability at z=1/300, a<=1")
    return probs


def certificate(mutant=None, verbose=True):
    require(mutant is None or mutant in MUTANTS, "known-mutant")
    th3, th6, keff = check_two_level(mutant)
    amp = check_star(mutant)
    theta, lhs, best = check_twist(mutant)
    var = check_static(mutant)
    probs = check_ceiling(mutant)
    res, dpert, bracket, g_ab = check_lattice(mutant, verbose)
    iso = check_isolated(mutant)
    if verbose:
        print(f"TWO_LEVEL: smooth zero at beta*A=pi; N-piece zeros N*arccosh(sec(pi/N)) > pi, decreasing; "
              f"3arccosh2={th3.str(12)}, 6arccosh(2/sqrt3)={th6.str(12)}")
        print(f"STAR: exact traceless amplitude {amp}*V, consecutive doubled angle 60 deg, fixed orientation; "
              f"effective kappa_* = {keff.str(12)}")
        for beta, klo, khi, dlo, dhi in res:
            print(f"LATTICE beta={beta} L={L_FROZEN} W_star={W_STAR}: det(1+U)>0 at beta*V={klo} "
                  f"({dlo.str(6, radius=False)}), <0 at beta*V={khi} ({dhi.str(6, radius=False)})")
        print(f"PERTURBED_REFERENCE beta=40: det(1+U) at beta*V=5 = {dpert.str(6, radius=False)} < 0")
        print(f"NEAR_ZERO beta=40: zero in beta*V in [{float(bracket[0]):.12f}, {float(bracket[1]):.12f}]; "
              f"|G_ba| = {abs(g_ab).str(6, radius=False)}; conditional pair corr = {(g_ab*g_ab).str(6, radius=False)}")
        for beta, klo, khi, dlo, dhi in iso:
            print(f"ISOLATED_STAR beta*t={beta}: det(1+U)>0 at beta*V={klo} ({dlo.str(4, radius=False)}), "
                  f"<0 at beta*V={khi} ({dhi.str(4, radius=False)})")
        print(f"THEOREM_W: theta = (2-cosh1)/(8 max(beta t,1/4)), (2-cosh1)/8 = {theta.str(10)}; "
              f"4-site Hubbard sanity |<D0+D2>| = {abs(lhs):.4e} <= {best:.4e}")
        print(f"LEMMA_C: posterior Var(n) at b=-U, t=0: {var}; Hessian (beta/U)(1-beta U/2) < 0 iff beta U > 2")
        print("CEILING: sigma^2 in (U Omega^2 beta/12)[1-a^2/60, 1]; exclusion prob lower bounds "
              + ", ".join(f"z={z}: {p.str(6, radius=False)}" for z, p in probs))
        print("SCOPE: finite exact/Arb witnesses; the analytic lemmas are in RIGOR.md")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mutant", choices=MUTANTS)
    parser.add_argument("--all-mutants", action="store_true")
    args = parser.parse_args()

    if args.mutant:
        try:
            certificate(mutant=args.mutant, verbose=False)
        except Rejected as exc:
            print(f"REJECTED {args.mutant}: {exc}")
            return 1
        print(f"UNCAUGHT {args.mutant}")
        return 0

    try:
        certificate(verbose=True)
    except Rejected as exc:
        print(f"BASELINE FAILURE: {exc}")
        return 1

    if not args.all_mutants:
        return 0

    script = str(Path(__file__).resolve())
    for mutant in MUTANTS:
        child = subprocess.run(
            [sys.executable, "-B", script, "--mutant", mutant],
            capture_output=True,
            text=True,
            check=False,
        )
        print(child.stdout.strip())
        if child.stderr:
            print(child.stderr.strip())
        print(f"mutant {mutant}: exit {child.returncode}")
        if child.returncode != 1 or not child.stdout.startswith(f"REJECTED {mutant}:"):
            return 1
    print(f"CERTIFIED: all {len(MUTANTS)} mutants rejected with exit 1")
    return 0


if __name__ == "__main__":
    sys.exit(main())
