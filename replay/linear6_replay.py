#!/usr/bin/env python3
"""LINEAR6 replay: exact (fractions.Fraction) and outward-interval (mpmath.iv) checks
of every CERTIFIED constant and identity in RIGOR.md / RESULT.md.

Usage
  replay.py                  run all exact/interval checks; exit 0 iff all pass
  replay.py --mutant NAME    run the checks with one deliberate error injected;
                             prints REJECTED and exits 1 when a check catches it
  replay.py --all-mutants    run every mutant in a child process; exit 0 iff every
                             child printed REJECTED and exited 1
  replay.py --list-mutants   list mutant names
  replay.py --numerical      additionally rerun the floating-point (NUMERICAL) evidence
                             quoted in RIGOR.md (needs numpy, scipy; takes minutes)

Run with
  UV_OFFLINE=1 UV_CACHE_DIR=D:/dev/.cache/uv uv run --no-project --with numpy \
      --with scipy --with mpmath python -B replay.py
"""
import subprocess
import sys
from fractions import Fraction as Fr
from math import comb

from mpmath import iv

iv.dps = 60

MUTANTS = {
    "c1_one": "A2: sandwich constant c1(a) replaced by 1 (equality-case weight w_1 = c1 fails)",
    "stage7r_reversed": "(1.4): claims STAGE7R constant (a/2)/sinh(a/2) exceeds c1(a)",
    "a3_coeff": "A3: Var-P0(1-P0) = q2(3-2q1-2q2) (wrong coefficient)",
    "a4_below_c1": "A4: claims c = c1/2 already attains the minimal norm U'",
    "a5_half_a2": "A5: subtracts 1/(2a^2) instead of 1/a^2",
    "a6_sign": "A6: same-sign cross term claimed positive",
    "t_g_over_4": "G3 Holstein row with t = g/4 instead of g/64",
    "qprime_8": "G3: attachment constant Q' = 6e claimed <= 8",
    "slice_1_12": "G3 Holstein slice bound claimed <= 1.12 (true value 1.1205)",
    "delta_ssh_250": "G3 SSH row: delta claimed <= 1/250",
    "two16": "Lemma T: 2^16 in place of 2^17 in hypothesis (ii)",
    "bg_quarter": "Theorem D': beta*g >= 1/4 in place of 4",
    "ionic_gap_noR": "Theorem D: Delta_st >= g+|mu| (static truncation R dropped)",
    "b2_sign": "B2: Gaussian term +(U/2beta)theta^2 entered with the wrong sign",
    "s1_2_18": "S1: exponent 2^18 beta U in place of 2^19 beta U",
    "s2_zero_mode": "S2: eps0 = t/g in place of t^2/g (no exact zero mode)",
    "cp_02": "(7.1) at t=0: claims c_P >= 0.2 at g*=1/4, beta U >= 8",
    "b5_markov_78": "B5: claims s(2beta-s)/beta^2 >= 7/8 on |s/beta-1| <= 1/2",
    "ionic_vs_l4_16": "Comparison: threshold Omega < g^{3/2}/(16 sqrt U) in place of 32",
    "odd_charge": "Admissible observables: single-fermion fields c, c^dag (odd charge) admitted",
    "complement_rate_half": "Complement (4.4): LINEAR4 rate claimed >= 1/2 above the D window (fixed-rate T_pair)",
    "omega1p_m_g": "Complement (4.6): Omega'_1 defined with m* = g instead of U (top of window (4.5))",
}

MUT = None


class Fail(Exception):
    pass


def need(cond, what):
    if cond is not True:
        raise Fail(what)


def lt(x, y):
    """Certified strict x < y for intervals or exact numbers."""
    xb = x.b if hasattr(x, "b") else x
    ya = y.a if hasattr(y, "a") else y
    return bool(xb < ya)


def I(v):
    if isinstance(v, Fr):
        return iv.mpf(v.numerator) / iv.mpf(v.denominator)
    return iv.mpf(v)


def isinh(x):
    return (iv.exp(x) - iv.exp(-x)) / 2


def icosh(x):
    return (iv.exp(x) + iv.exp(-x)) / 2


PI = iv.pi
E = iv.e
A_GRID = ["1e-6", "1e-3", "0.1", "0.5", "1", "2", "5", "10", "50", "100", "1000"]
RATS = [Fr(1, 7), Fr(2, 3), Fr(5, 4), Fr(17, 5), Fr(1, 1000), Fr(31, 2)]


# ---------------------------------------------------------------- section 1
def chk_14():
    """(1.4) and the series fact pi^2 > 6 used in its proof."""
    need(lt(I(6), PI * PI), "pi^2 > 6")
    for s in A_GRID:
        a = iv.mpf(s)
        x = a / 2
        c1 = 4 * PI**2 / (4 * PI**2 + a**2)
        if MUT == "c1_one":
            c1 = iv.mpf(1)
        s7 = x / isinh(x)
        if MUT == "stage7r_reversed":
            need(lt(c1, s7), f"(1.4) reversed at a={s}")
        need(lt(s7, c1), f"(1.4) left at a={s}")
        need(lt(I(1), x * icosh(x) / isinh(x)), f"(1.4) right at a={s}")
    return "(1.4) on 11 values of a in [1e-6,1e3]; pi^2>6"


def chk_A1_A2_weights():
    """U - D_k = U nu^2/(nu^2+Om^2); weights increase in |k|; w_{+-1} = c1."""
    for U in RATS[:3]:
        for nu2 in RATS:
            for om2 in RATS:
                need(U - U * om2 / (nu2 + om2) == U * nu2 / (nu2 + om2), "A1 multiplier")
    # with P := pi^2 (any positive stand-in), w_k = 4P k^2/(4P k^2 + a^2)
    for P in [Fr(98696, 10000), Fr(98697, 10000), Fr(10)]:
        for a2 in RATS:
            w = [Fr(4) * P * k * k / (4 * P * k * k + a2) for k in range(1, 60)]
            need(all(w[i] < w[i + 1] for i in range(len(w) - 1)), "A2 weights increase")
            c1 = 4 * P / (4 * P + a2)
            if MUT == "c1_one":
                c1 = Fr(1)
            need(w[0] == c1, "A2 equality case w_1 = c1")
    return "A1 multiplier identity; A2 weight monotonicity and w_1=c1 (exact)"


def chk_A3():
    def var(q1, q2):
        return q1 + 4 * q2 - (q1 + 2 * q2) ** 2

    grid = [Fr(i, 6) for i in range(-6, 13)]
    for q1 in grid:
        for q2 in grid:
            q0 = 1 - q1 - q2
            P0, P1, P2 = 1 - q0, 1 - q1, 1 - q2
            r0 = q2 * (3 - 2 * q1 - 3 * q2)
            if MUT == "a3_coeff":
                r0 = q2 * (3 - 2 * q1 - 2 * q2)
            need(var(q1, q2) - P0 * (1 - P0) == r0, "A3 line 1")
            need(var(q1, q2) - P1 * (1 - P1) == 4 * q0 * q2, "A3 line 2")
            need(var(q1, q2) - P2 * (1 - P2) == q0 * (3 - 2 * q1 - 3 * q0), "A3 line 3")
            need(3 - 2 * q1 - 3 * q2 == 3 * q0 + q1, "A3 nonnegativity rewrite")
    return "A3 three identities on a 19x19 rational grid (exact)"


def chk_A4():
    for P in [Fr(394784, 10000), Fr(394785, 10000), Fr(40)]:  # P = 4 pi^2 stand-ins
        for a2 in RATS:
            U = Fr(3, 2)
            Up = U * a2 / (P + a2)
            c1 = P / (P + a2)
            need((1 - c1) * U == Up, "A4 balance at c1")
            need(U - Up == c1 * U, "(1.7) static coefficient")
            for j in range(0, 41):
                c = Fr(j, 40)
                norm = max((1 - c) * U, Up)
                need(norm >= Up, "A4 lower bound")
                need((norm == Up) == (c >= c1), "A4 minimizer set is [c1,1]")
            if MUT == "a4_below_c1":
                c = c1 / 2
                need(max((1 - c) * U, Up) == Up, "A4 mutant: c1/2 attains U'")
        # U' = tau  <=>  tau^3 = (U - tau) Om^2 / P
        for U in RATS[:4]:
            for tau in RATS:
                for Om in RATS[:4]:
                    a2 = (Om / tau) ** 2
                    lhs = U * a2 / (P + a2) - tau
                    rhs = (U * Om**2 - P * tau**3 - Om**2 * tau) / (P * tau**2 + Om**2)
                    need(lhs == rhs, "A4 remark identity")
    return "A4 norm formula, minimizer set [c1,1], remark identity (exact, P=4pi^2 stand-ins)"


def chk_A5():
    for x in [Fr(i, 12) for i in range(13)]:
        y = x - Fr(1, 2)
        need(y * y / 2 - Fr(1, 24) == (x * x - x + Fr(1, 6)) / 2, "A5 constant term")
    vals = [abs((x * x - x + Fr(1, 6)) / 2) for x in [Fr(i, 1200) for i in range(1201)]]
    need(max(vals) == Fr(1, 12), "A5 sup = 1/12")
    # outward check of the a -> 0 limit with explicit O(a^2) remainder at a = 1e-3
    a = iv.mpf("1e-3")
    sub = 1 / a**2
    if MUT == "a5_half_a2":
        sub = 1 / (2 * a**2)
    for k in range(11):
        x = iv.mpf(k) / 10
        val = icosh(a * (x - iv.mpf(1) / 2)) / (2 * a * isinh(a / 2)) - sub
        lim = (x * x - x + iv.mpf(1) / 6) / 2
        need(lt(abs(val - lim), a * a), f"A5 limit at x={k/10}")
    return "A5 constant term, sup 1/12, limit at a=1e-3 within a^2 (exact + interval)"


def chk_A6():
    for s in A_GRID:
        a = iv.mpf(s)
        D0 = icosh(a / 2) / isinh(a / 2) / 2  # D(0)/(U Omega)
        Dh = 1 / (2 * isinh(a / 2))  # D(beta/2)/(U Omega)
        need(lt(Dh, D0) and lt(I(0), Dh), f"A6 kernel extrema at a={s}")
        cross_upper = -Dh  # same-sign cross term <= -D(beta/2)|d1||d2| < 0
        if MUT == "a6_sign":
            cross_upper = Dh
        need(lt(cross_upper, I(0)), "A6 same-sign cross term negative")
    return "A6 0 < D(beta/2) < D(0) and same-sign cross term < 0 (interval)"


# ---------------------------------------------------------------- section 2
def slice_bound(q):
    return 1 + q / 2 + q * q / (2 * (1 - q))


def chk_G3():
    # Holstein row: t <= g/64, r = g/8, g = 1
    t = Fr(1, 64) if MUT != "t_g_over_4" else Fr(1, 4)
    r = Fr(1, 8)
    g0 = 1 - r
    q = 4 * E * I(t) / I(g0)
    need(lt(q, I(Fr(1, 4))), "Holstein q < 1/4")
    if MUT != "t_g_over_4":
        need(t / g0 == Fr(1, 56), "q = e/14")
    need(slice_bound(Fr(1, 4)) == 1 + Fr(1, 8) + Fr(1, 24), "(2.3) at q=1/4")
    sb = slice_bound(q)
    cap = I(Fr(1121, 1000)) if MUT != "slice_1_12" else I(Fr(112, 100))
    need(lt(sb, cap), "Holstein slice bound")
    gh = 1 - 4 * t
    need(gh >= Fr(15, 16), "Holstein g_h")
    rho = r / Fr(15, 16)
    need(rho == Fr(2, 15), "Holstein rho_V")
    delta = (rho / 2) / (1 - rho)
    need(delta == Fr(1, 13), "Holstein delta")
    # SSH row: g >= 60(t+R), r_b = g/512, g = 1, worst case t+R = 1/60
    qs = 4 * E * (I(Fr(1, 60)) + I(Fr(1, 512)))
    need(lt(qs, I(Fr(203, 1000))), "SSH q < 0.203")
    need(lt(slice_bound(qs), I(Fr(113, 100))), "SSH slice <= 1.13")
    ghs = 1 - Fr(4, 60)
    need(ghs == Fr(14, 15), "SSH g_h")
    rhos = Fr(4, 512) / ghs
    need(rhos == Fr(15, 1792), "SSH rho_V")
    ds = (rhos / 2) / (1 - rhos)
    need(ds == Fr(15, 3554), "SSH delta value")
    need(ds <= (Fr(1, 200) if MUT != "delta_ssh_250" else Fr(1, 250)), "SSH delta bound")
    Qp = 2 * (1 + 2 // 1) * E
    need(lt(Qp, I(17 if MUT != "qprime_8" else 8)), "Q' = 6e <= 17")
    # attachment: sum_y e^{-d} on Z^2 = coth(1/2)^2 <= (1+2/1)^2, and Cauchy-Schwarz 2*3 = 6
    cth = icosh(I(Fr(1, 2))) / isinh(I(Fr(1, 2)))
    need(lt(cth, I(3)), "sum_n e^{-|n|} = coth(1/2) <= 1+2/a at a=1")
    for k in range(0, 41):
        need(sum(comb(k, j) ** 2 for j in range(k + 1)) == comb(2 * k, k) <= 4**k, "minors 4^k")
    return "G3 table: Holstein q<=e/14, slice<1.121, rho=2/15, delta=1/13; SSH q<0.203, slice<1.13, rho=15/1792, delta=15/3554<=1/200; Q'=6e<17; 4^k"


# ---------------------------------------------------------------- section 3
def chk_T():
    c_over_lam_a2 = 64
    two17 = 2**17 if MUT != "two16" else 2**16
    need(2048 * c_over_lam_a2 == two17, "Lemma T: (ii) <=> l >= 2048 c")
    need(2048 * Fr(1, 128) == 16, "Lemma T: 2048 l >= 16")
    need(1024 - 2 * 2048 == -3072, "Lemma T: f' <= -3072")
    need(Fr(8192, 3072) == Fr(8, 3) and 8 + Fr(8, 3) <= 12, "Lemma T: 8+8/3 <= 12")
    need(64 * 2**17 == 2**23, "LINEAR4 consistency: tau^3 >= 2^23 U Om^2")
    need(8 * 2**17 == (2**10) ** 2, "(3.4): 2^10")
    need(512 * 2**17 == (2**13) ** 2, "D': 2^13")
    need(Fr(1, 128) / Fr(1, 8) == Fr(1, 16), "(i) at r=g/8 <=> beta g >= 1/16 <=> tau <= 16 g")
    bg_min = Fr(4) if MUT != "bg_quarter" else Fr(1, 4)
    need(bg_min / 512 == Fr(1, 128), "(i) at r_b=g/512 <=> beta g >= 4")
    need(60 * 4096 == 245760 and 245760 >= 4, "D': beta g >= 245760")
    need(Fr(2048, 8) == 256 and Fr(1024, 8) == 128 and Fr(1024, 512) == 2, "tail exponents")
    # threshold equality at tau = 2^10 sqrt(U/g) Omega with U=1, g=1/4 (sqrt(U/g)=2)
    U, g, Om = Fr(1), Fr(1, 4), Fr(3, 7)
    tau = 2**11 * Om
    need(two17 * U * Om**2 / tau**2 == g / 8, "(3.4) equality at the threshold")
    return "Lemma T algebra, LINEAR4 reduction 2^23, (3.4) 2^10, D' 2^13, beta g thresholds (exact)"


def chk_totals():
    T = 12 * iv.exp(iv.mpf(-128))
    need(lt(T, I(Fr(1, 2))), "12 e^-128 <= 1/2")
    need(4 * 12 == 48 and 17 * 48 == 816 and 17 * 8 == 136, "products")
    need(Fr(4096**2, 8) == 2**21, "static envelope exponent 2^21")
    sC = 816 * iv.exp(iv.mpf(-128))
    need(lt(sC, I(Fr(1, 64))), "Theorem C: 816 e^-128 < 1/64")
    sD = sC + 136 * iv.exp(-iv.mpf(2**21))
    need(lt(sD, I(Fr(1, 64))), "Theorem D total < 1/64")
    need(2 * 245760 == 491520, "2 beta g exponent")
    sDp = 17 * (48 * iv.exp(-iv.mpf(491520)) + 8 * iv.exp(-iv.mpf(2**21)))
    need(lt(sDp, I(Fr(1, 64))), "Theorem D' total < 1/64")
    # Theorem D gap with static truncation R = 4096 U
    U, g, mu = Fr(1), Fr(1, 4), Fr(-1)
    R = 4096 * U
    Dst = g + R + abs(mu) if MUT != "ionic_gap_noR" else g + abs(mu)
    worst_gap = Dst - R - abs(mu)
    need(worst_gap >= g, "Theorem D gap |W| >= g")
    return "Theorem C/D/D' totals < 1/64 (interval); Theorem D gap"


# ---------------------------------------------------------------- sections 5-7
def chk_B2():
    for kap in RATS[:4]:
        for U in RATS[1:4]:
            for beta in RATS[2:5]:
                for Bp, Bm in [(1, 0), (0, 1), (3, 2), (5, 7)]:
                    B = Bp + Bm
                    sg = 1 if MUT != "b2_sign" else -1

                    def ex(th):
                        return (-th * kap * Bp - th * (2 * U + kap) * Bm
                                + sg * U / (2 * beta) * th * th * B + 2 * U * th * Bm)

                    th = beta * kap / U
                    need(ex(th) == -beta * kap**2 * B / (2 * U), "B2 completed square")
                    h = Fr(1, 10**6)
                    need(ex(th + h) >= ex(th) and ex(th - h) >= ex(th), "B2 optimal theta")
    return "B2 exponent at theta=beta kappa/U equals -beta kappa^2|B|/(2U) and is minimal (exact)"


def chk_S1():
    p19 = 2**19 if MUT != "s1_2_18" else 2**18
    for U in RATS[1:5]:
        for w in RATS:
            z = 1024 * U + w
            need(1024 * z - z * z / (2 * U) == p19 * U - w * w / (2 * U), "S1 identity")
    return "S1 exponent identity 1024z - z^2/(2U) = 2^19 U - w^2/(2U) (exact)"


def chk_S2():
    for t in [Fr(1, 100), Fr(1, 64), Fr(1, 10), Fr(1, 2), Fr(1)]:
        g = Fr(1)
        eps0 = t * t / g if MUT != "s2_zero_mode" else t / g
        h = [[eps0, -t], [-t, g]]
        p0 = [g, t]
        p1 = [t, -g]
        hp0 = [h[0][0] * p0[0] + h[0][1] * p0[1], h[1][0] * p0[0] + h[1][1] * p0[1]]
        need(hp0 == [0, 0], "S2 exact zero mode")
        lam = g + t * t / g
        hp1 = [h[0][0] * p1[0] + h[0][1] * p1[1], h[1][0] * p1[0] + h[1][1] * p1[1]]
        need(hp1 == [lam * p1[0], lam * p1[1]], "S2 second eigenpair")
        n2 = g * g + t * t
        need((p0[1] * p1[1]) ** 2 / (n2 * n2) == (g * t / n2) ** 2, "S2 overlap gt/(g^2+t^2)")
    need(Fr(1, 9) + 1 == Fr(10, 9), "S2 (3pi)^-2 + pi^-2 = (10/9) pi^-2")
    return "S2 exact zero mode, eigenpair, overlap, Matsubara factor (exact)"


def chk_S3():
    for U in RATS[1:4]:
        for beta in RATS[1:4]:
            e1 = -beta * (U * U - 0) / (2 * U)
            e2 = -beta * (4 * U * U - U * U) / (2 * U)
            need(e1 == -beta * U / 2 and e2 == -3 * beta * U / 2, "S3 Gaussian ratios")
            need((e1, e1 + 2 * beta * U) == (-beta * U / 2, 3 * beta * U / 2), "S3 first range")
            need((e2, e2 + 2 * beta * U) == (-3 * beta * U / 2, beta * U / 2), "S3 second range")
            # extremal m = delta_{2 beta}: log L(b) = -2 beta b - 2 U beta
            need((-2 * beta * (-U)) - 0 == 2 * beta * U, "S3 upper end attained")
    return "S3 ratio ranges [e^-bU/2, e^3bU/2], [e^-3bU/2, e^bU/2] and extremals (exact)"


def chk_Ploc_t0():
    for lam in [Fr(0), Fr(1, 3), Fr(1), Fr(8), Fr(40)]:
        x = iv.exp(-I(lam) / 2)
        p1 = 2 * x / (2 + 2 * x)
        need(not lt(x, p1), "P(n=1) <= e^{-beta U/2}")
    cP = I(Fr(9, 32)) - iv.log(iv.mpf(2)) / 8
    need(lt(I(Fr(19, 100) if MUT != "cp_02" else Fr(20, 100)), cP), "c_P >= 0.19")
    return "(7.1) at t=0: P(n=1) <= e^{-bU/2}; c_P = 9/32 - log2/8 > 0.19 (interval)"


def chk_B5():
    lo = Fr(3, 4) if MUT != "b5_markov_78" else Fr(7, 8)
    for i in range(-50, 51):
        sig = 1 + Fr(i, 100)
        need(sig * (2 - sig) == 1 - (sig - 1) ** 2, "B5 identity")
        need(sig * (2 - sig) >= lo, "B5 Markov constant")
    need([n * (2 - n) for n in (0, 1, 2)] == [0, 1, 0], "n(2-n) = 1_{n=1}")
    for xs in ["0.01", "0.5", "2", "30"]:
        for ys in ["0.02", "0.7", "3", "100"]:
            x, y = iv.mpf(xs), iv.mpf(ys)
            lm = (x - y) / (iv.log(x) - iv.log(y))
            need(not lt((x + y) / 2, lm), "log mean <= arithmetic mean")
    return "B5 Markov constant 3/4, n(2-n), log-mean <= arithmetic mean (exact + interval)"


def chk_compare():
    k = 32 if MUT != "ionic_vs_l4_16" else 16
    need(Fr(2**30, 2**25) == k, "D beats LINEAR4 iff Omega < g^{3/2}/(32 sqrt U)")
    need(Fr(2**39, 2**25) == 2**14, "D' beats LINEAR4 iff Omega < g^{3/2}/(2^14 sqrt U)")
    # instance: U=1, g=1/4, t=1/256 (= g/64), Omega=1e-6, E = t (Holstein)
    need(Fr(1, 256) <= Fr(1, 4) / 64, "instance satisfies t <= g/64")
    U, g, t, Om = iv.mpf(1), iv.mpf(1) / 4, iv.mpf(1) / 256, iv.mpf("1e-6")
    tauD = 2**10 * iv.sqrt(U / g) * Om
    need(lt(tauD, g) and lt(tauD, U), "instance: tau_D <= min(g,U), so (4.3) is nonempty")
    # Omega <= Omega_1: LINEAR4 criterion holds at m* = min(g,U) = 1/4
    ms = g
    need(lt(2**25 * U * Om**2 * (1 + iv.log(1 + t / ms)), ms**3), "instance: Omega <= Omega_1")
    # LINEAR4 criterion fails at tau = tau_D, so (4.4) is strictly below what LINEAR4 certifies
    lhs = tauD**3
    rhs = 2**25 * U * Om**2 * (1 + iv.log(1 + t / tauD))
    need(lt(lhs, rhs), "LINEAR4 criterion fails at tau_D (D strictly better there)")
    return "comparison constants 32 and 2^14; explicit ionic instance with Omega <= Omega_1 (exact + interval)"


def chk_complement():
    # admissible observables (RIGOR §4): opposite nonzero EVEN charge, as LINEAR4 §§4-5 require
    q = {"c": -1, "cdag": 1}
    A, B = (["c"], ["cdag"]) if MUT == "odd_charge" else (["c", "c"], ["cdag", "cdag"])
    qA, qB = sum(q[f] for f in A), sum(q[f] for f in B)
    need(qA != 0 and qA + qB == 0 and qA % 2 == 0, "pair fields have opposite nonzero even charge")
    # Definition (4.0): above the D window LINEAR4 gives only rate a/2 < 1/2 (Holstein instance, tau=0.3 > m*=1/4)
    U, g, t = iv.mpf(1), iv.mpf(1) / 4, iv.mpf(1) / 256
    tau = iv.mpf(3) / 10
    a = iv.mpf(1) if MUT == "complement_rate_half" else 1 / (1 + 16 * t / tau)
    need(lt(0, a) and lt(a, 1), "LINEAR4 rate a/2 in (0,1/2) at tau=0.3: rate must be tau-dependent")
    need(lt(abs(a - iv.mpf("0.8276")), iv.mpf("1e-4")), "a(0.3) = 0.8276 on the instance")
    # SSH: 16 J/tau > 16*4096, so a < 1/65537 < 1.526e-5 for every tau > 0
    need(Fr(1, 1 + 16 * 4096) < Fr(1526, 10**8), "SSH LINEAR4 rate a < 1/65537 < 1.526e-5")
    # Omega_1 (4.4): m* = min(g,U); closed form solves m*^3 = 2^25 U Om^2 [1+log(1+t/m*)]; window nonempty
    for U, g in [(Fr(1), Fr(1, 4)), (Fr(1, 8), Fr(1)), (Fr(8), Fr(4)), (Fr(1), Fr(1000))]:
        Ui, gi = I(U), I(g)
        ti = gi / 64
        ms = gi if g < U else Ui
        br = 1 + iv.log(1 + ti / ms)
        om1 = ms ** iv.mpf(1.5) / (iv.mpf(2) ** iv.mpf(12.5) * iv.sqrt(Ui) * iv.sqrt(br))
        need(lt(abs(ms**3 - 2**25 * Ui * om1**2 * br), ms**3 * iv.mpf("1e-40")), "Omega_1 closed form")
        lo = 2**10 * iv.sqrt(Ui / gi) * om1
        need(lt(lo, ms * iv.mpf(2) ** iv.mpf(-2.5) * (1 + iv.mpf("1e-40"))), "2^10 sqrt(U/g) Omega_1 <= 2^-5/2 m*")
    need(lt(iv.mpf("1e-6"), (iv.mpf(1) / 4) ** iv.mpf(1.5) / (iv.mpf(2) ** iv.mpf(12.5)
         * iv.sqrt(1 + iv.log(1 + iv.mpf(1) / 64)))), "instance Omega = 1e-6 <= Omega_1")
    # Omega'_1 (4.6): m* = U (top of (4.5)), g >= 60(t+4096U); criterion at tau=U and window nonempty
    for U in [Fr(1, 8), Fr(1), Fr(8)]:
        for t in [U / 100, U, 10 * U]:
            g = 60 * (t + 4096 * U)
            Ui, ti, gi = I(U), I(t), I(g)
            ms = gi if MUT == "omega1p_m_g" else Ui
            om1p = ms ** iv.mpf(1.5) / (iv.mpf(2) ** iv.mpf(12.5) * iv.sqrt(Ui)
                                         * iv.sqrt(1 + iv.log(1 + (ti + Ui) / ms)))
            rhs = 2**25 * Ui * om1p**2 * (1 + iv.log(2 + ti / Ui))
            need(lt(rhs, Ui**3 * (1 + iv.mpf("1e-40"))), "LINEAR4 SSH criterion holds at tau=U for Omega=Omega'_1")
            lo = 2**13 * iv.sqrt(Ui / gi) * om1p
            need(lt(lo, iv.sqrt(2) * Ui / iv.sqrt(245760) * (1 + iv.mpf("1e-40"))) and lt(lo, Ui),
                 "2^13 sqrt(U/g) Omega'_1 <= sqrt2 U/sqrt(245760) < U")
    return "even pair charge; LINEAR4 rate a/2 < 1/2 above windows; Omega_1, Omega'_1 explicit; windows nonempty (interval)"


CHECKS = [chk_14, chk_A1_A2_weights, chk_A3, chk_A4, chk_A5, chk_A6, chk_G3, chk_T,
          chk_totals, chk_B2, chk_S1, chk_S2, chk_S3, chk_Ploc_t0, chk_B5, chk_compare,
          chk_complement]


def run_checks():
    for c in CHECKS:
        try:
            msg = c()
        except Fail as e:
            return False, f"{c.__name__}: {e}"
        print(f"PASS {c.__name__}: {msg}")
    return True, ""


# ---------------------------------------------------------------- numerical evidence
def numerical():
    import numpy as np
    import scipy.linalg as sl
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    import mpmath as mp

    print("NUMERICAL A1: Matsubara series (|k|<=2e5) vs closed form")
    worst = 0.0
    for a in [1e-6, 1e-3, 0.1, 1, 10, 100]:
        for p in [0.01, 0.3, 0.5]:
            k = np.arange(1, 200001)
            nu = 2 * np.pi * k
            ser = 0.5 * 2 * np.sum(4 * np.sin(np.pi * k * p) ** 2 / (nu**2 + a * a))
            cl = float(mp.sinh(a * p / 2) * mp.sinh(a * (1 - p) / 2) / mp.sinh(a / 2)) / a
            worst = max(worst, abs(ser - cl) / cl)
    print(f"  max relative error {worst:.3e} (claim < 1e-4)")
    assert worst < 1e-4

    print("NUMERICAL Lemma T with the LINEAR4 tail majorant at r = 2^17 U a^2")
    mp.mp.dps = 30
    for rr in [1 / 32, 1 / 8]:
        for beta in [1.0, 10.0, 100.0]:
            a2 = rr / 2**17
            c = 64 * beta * a2
            l = beta * rr
            f = lambda y: mp.e ** (1024 * y) * min(1, 8 * mp.e ** (-y**2 / c))
            T = mp.e ** (1024 * l) * min(1, 8 * mp.e ** (-l**2 / c)) + 1024 * mp.quad(
                f, [l, l + 50 * mp.sqrt(c) + 10, mp.inf])
            B = 12 * mp.e ** (-1024 * l)
            print(f"  r/U={rr} betaU={beta}: T={mp.nstr(T, 6)} bound={mp.nstr(B, 6)}")
            assert T <= B

    print("NUMERICAL B1 displacement identity and B3 complete monotonicity (2 sites, N_ph=60)")

    def ops(nm):
        I2 = np.eye(2)
        Z = np.diag([1.0, -1.0])
        an = np.array([[0.0, 1.0], [0.0, 0.0]])
        cs = []
        for j in range(nm):
            m = np.array([[1.0]])
            for kk in range(nm):
                m = np.kron(m, Z if kk < j else (an if kk == j else I2))
            cs.append(m)
        return cs

    cc = ops(4)
    Ie = np.eye(16)

    def model(t, Om, U, beta, mu, Nph):
        F = np.sqrt(U * Om**2)
        an = np.diag(np.sqrt(np.arange(1, Nph)), 1)
        X = (an + an.T) / np.sqrt(2 * Om)
        hop = -t * sum(cc[i].T @ cc[i + 2] + cc[i + 2].T @ cc[i] for i in (0, 1))
        n1 = cc[0].T @ cc[0] + cc[1].T @ cc[1]
        n2 = cc[2].T @ cc[2] + cc[3].T @ cc[3]
        Ip = np.eye(Nph)
        H = np.kron(hop - mu * (n1 + n2), Ip) + np.kron(Ie, Om * an.T @ an) + F * np.kron(n1, X)
        w0 = np.linalg.eigvalsh(H).min()

        def tr(M):
            return np.sum(np.exp(-beta * (np.linalg.eigvalsh(M) - w0)))

        return H, F * np.kron(Ie, X), np.kron(n1, Ip), tr

    for (t, Om, U, beta, mu, th) in [(0.3, 1.0, 2.0, 2.0, -1.0, 0.7), (0.5, 0.5, 1.0, 3.0, -0.5, -1.2),
                                     (0.1, 2.0, 3.0, 1.0, -1.5, 2.0)]:
        H, b1, N1, tr = model(t, Om, U, beta, mu, 60)
        Z = tr(H)
        lhs = tr(H - (th / beta) * b1) / Z
        rhs = np.exp(U * th**2 / (2 * beta)) * tr(H + (U * th / beta) * N1) / Z
        print(f"  (t,Om,U,beta,mu,theta)={(t, Om, U, beta, mu, th)}: rel err {abs(lhs - rhs) / rhs:.2e}")
        assert abs(lhs - rhs) / rhs < 1e-12
    H, b1, N1, tr = model(0.4, 0.7, 2.0, 4.0, -1.0, 60)
    Z = tr(H)
    Psi = lambda e: tr(H + e * N1) / Z
    for e0 in [-1.0, 0.0, 0.8]:
        vals = [Psi(e0 + 0.3 * k) for k in range(10)]
        for k in range(1, 9):
            d = sum((-1) ** j * comb(k, j) * vals[j] for j in range(k + 1))
            assert d >= -1e-12, ("CM", e0, k, d)
    es = np.linspace(-1, 1, 21)
    ps = [Psi(e) for e in es]
    assert all(np.diff(ps) < 0) and all(np.diff([p * np.exp(8.0 * e) for p, e in zip(ps, es)]) > 0)
    print("  alternating forward differences k<=8 at 3 base points; monotonicities hold")

    print("NUMERICAL S2: Matsubara-truncated ||C0 V|| vs (6.2), 800 frequencies")
    for beta in (10.0, 100.0, 1000.0):
        for t in (1 / 100, 1 / 64, 1 / 10):
            g = 1.0
            r = g / 8
            h = np.array([[t * t / g, -t], [-t, g]])
            ns = np.arange(-400, 400)
            w = (2 * ns + 1) * np.pi / beta
            N = len(ns)
            C0 = np.zeros((2 * N, 2 * N), complex)
            V = np.zeros((2 * N, 2 * N), complex)
            for i, wi in enumerate(w):
                C0[2 * i:2 * i + 2, 2 * i:2 * i + 2] = np.linalg.inv(1j * wi * np.eye(2) + h)
            for i in range(N):
                for j in (i - 1, i + 1):
                    if 0 <= j < N:
                        V[2 * i + 1, 2 * j + 1] = r / 2
            nrm = np.linalg.norm(C0 @ V, 2)
            lb = np.sqrt(10) / 3 * r * beta / (2 * np.pi) * g * t / (g * g + t * t)
            print(f"  beta g={beta:g} t/g={t:g}: norm {nrm:.4f} bound {lb:.4f}")
            assert nrm >= lb

    print("NUMERICAL G3: 4x4 torus, exact block-cyclic equal-time Green function")
    rng = np.random.default_rng(7)
    L = 4
    N = L * L
    ix = lambda i, j: (i % L) * L + (j % L)
    A = np.zeros((N, N))
    for i in range(L):
        for j in range(L):
            A[ix(i, j), ix(i + 1, j)] = A[ix(i + 1, j), ix(i, j)] = 1
            A[ix(i, j), ix(i, j + 1)] = A[ix(i, j + 1), ix(i, j)] = 1
    Dm = np.zeros((N, N), int)
    for x in range(N):
        for y in range(N):
            dx, dy = abs(x // L - y // L), abs(x % L - y % L)
            Dm[x, y] = min(dx, L - dx) + min(dy, L - dy)

    def green(h, Vs, beta):
        # block-cyclic antiperiodic Trotter system, solved sparsely (dense N*M square is ~0.5 GB at beta=60)
        M = len(Vs)
        dt = beta / M
        Bs = [sl.expm(-dt * (h + np.diag(v))) for v in Vs]
        blocks = [[None] * M for _ in range(M)]
        for kk in range(M):
            blocks[kk][kk] = sp.identity(N, format="csr")
        blocks[0][M - 1] = sp.csr_matrix(Bs[M - 1])
        for kk in range(1, M):
            blocks[kk][kk - 1] = sp.csr_matrix(-Bs[kk - 1])
        O = sp.bmat(blocks, format="csc")
        Ev = np.zeros((N * M, N))
        Ev[:N] = np.eye(N)
        return spla.spsolve(O, Ev)[:N]

    worst = 0.0
    worstslice = 0.0
    for trial in range(60):
        g, t, r = 1.0, 1 / 64, 1 / 8
        beta = rng.choice([1.0, 5.0, 20.0, 60.0])
        mode = trial % 4
        if mode == 0:
            sg = rng.choice([-1, 1], N)
        elif mode == 1:
            sg = np.array([(-1) ** (x // L + x % L) for x in range(N)])
        elif mode == 2:
            sg = np.ones(N)
        else:
            sg = np.where(np.arange(N) < N // 2, 1, -1)
        W = sg * g * (1 + rng.uniform(0, 0.2, N) * (trial % 2))
        h = np.diag(W) - t * A
        M = int(max(40, 8 * beta))
        if trial % 3 == 0:
            Vs = [rng.uniform(-r, r, N) for _ in range(M)]
        elif trial % 3 == 1:
            Vs = [-sg * r * (1 if (kk // max(1, M // 7)) % 2 == 0 else -1) for kk in range(M)]
        else:
            Vs = [-sg * r * np.ones(N) for kk in range(M)]
        G = green(h, Vs, beta)
        worst = max(worst, (np.abs(G) * np.exp(Dm)).max())
        for y in range(N):
            Ey = np.diag(np.exp(Dm[:, y]))
            worstslice = max(worstslice, np.linalg.norm(Ey @ G @ np.linalg.inv(Ey), 2))
    print(f"  max conjugated slice norm {worstslice:.4f} (bound 1.121); max |G_xy| e^d {worst:.4f}")
    assert worstslice < 1.121


def main(argv):
    global MUT
    if "--list-mutants" in argv:
        for k, v in MUTANTS.items():
            print(f"{k}: {v}")
        return 0
    if "--all-mutants" in argv:
        bad = []
        for name in MUTANTS:
            p = subprocess.run([sys.executable, "-B", __file__, "--mutant", name],
                               capture_output=True, text=True)
            ok = p.returncode == 1 and "REJECTED" in p.stdout
            print(f"{'OK ' if ok else 'BAD'} mutant {name}: exit {p.returncode}; "
                  f"{p.stdout.strip().splitlines()[-1] if p.stdout.strip() else ''}")
            if not ok:
                bad.append(name)
        print(f"{len(MUTANTS) - len(bad)}/{len(MUTANTS)} mutants rejected")
        return 0 if not bad else 1
    if "--mutant" in argv:
        name = argv[argv.index("--mutant") + 1]
        if name not in MUTANTS:
            print(f"unknown mutant {name}")
            return 2
        MUT = name
        ok, why = run_checks()
        if ok:
            print(f"MUTANT SURVIVED {name}")
            return 0
        print(f"REJECTED {name} by {why}")
        return 1
    ok, why = run_checks()
    if not ok:
        print(f"FAIL {why}")
        return 1
    if "--numerical" in argv:
        numerical()
        print("NUMERICAL evidence reproduced")
    print(f"ALL {len(CHECKS)} CHECK GROUPS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
