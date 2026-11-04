* stpaf_cert_shapley.do  --  certification of combined() (average sequential / Shapley PAFs) and
* r(V_paf) in stpaf 0.12.0 (test 20)
* Run with stpaf.ado on the adopath (e.g. -adopath ++ "path/to/folder"-), or from run_stpaf_0120.do.
*
* Model A (data _csh_sim: recurrent events, 2 strata, 30% delayed entry, x time-varying 0-3):
*   stcox i.x##c.age c.z##c.w, strata(stratum) vce(cluster id) breslow
* factors x (set(x=0)), z (eliminate(z)), w (cap(w=0)); z and w interact (continuous interaction).
*
* What is checked
*  20a. K = 2 (z, w), explicit leave-one-cluster-out loop of stcox refits (20 clusters): Shapley
*       PAFs = average over the 2 orders of the subset PAFs; full jackknife covariance of all rows by
*       hand = r(V_paf) of refit (1e-8) and of the Newton jackknife (1e-4, scaled); Shapley SEs
*  20b. K = 3, influence function: every subset row = separate eliminate()/set()/cap() call (1e-10);
*       Shapley PAFs = average over all 3! orders by hand (1e-12); sum = joint (1e-12), r(sumdiff);
*       = r(W_shap) x subset PAFs; diag(r(V_paf)) = se^2; Shapley block of r(V_paf) = W V_sub W';
*       r(V_paf) of an each call = the matching block; labels, r(K), r(factors), r(subsets)
*  20c. grouped factor: indicator model with combined((eliminate(d1) eliminate(d2 d3)) ...) = i.x
*       model with set(x=0) (IF and jackknife, table and V); blanks, abbreviation; one variable per
*       factor outside parentheses (combined(eliminate(z w)) = combined(eliminate(z) set(w=0)))
*  20d. efron ties with tied failure times, if subset: IF and Newton jackknife (400 clusters), Shapley
*       CIs (normal / t(G-1), PAF scale), E_cf and se_logR missing; IF vs jackknife covariances;
*       Newton = refit (20 clusters)
*  20e. r(V_paf) in the other modes: varlevels() (joint row/column = sum of the level rows/columns;
*       IF and jackknife); single scenario diag = se^2; fweights without clusters = expanded data
*  20f. after mi estimate: per-imputation loop and hand Rubin (PAF scale for the Shapley rows,
*       log scale for the others), multivariate Rubin for r(V_paf), r(paf_mi), r(se_mi),
*       per-imputation additivity, r(sumdiff) with citrans(log) and (none); jackknife: BR df by hand
*  20g. errors (198 / 111 / 7) and K limit
*  20h. hand computation (Mata, columns built from the data): model c.w##i.z (continuous variable
*       in an interaction with a factor variable): cap(w=0), set(z=0), both, and their Shapley PAFs
*  20i. timings K = 6 and K = 8 (IF and default jackknife), additivity

version 15.0
clear all
set more off
set type double
set seed 20261105

* ----------------------------------------------------------------------------- helpers
capture program drop _csh_sim
program define _csh_sim
    syntax, n(integer) [ROUND]
    clear
    set obs `n'
    gen long id = _n
    gen double age = rnormal()
    gen byte stratum = 1 + (runiform() > .5)
    gen byte xb = floor(3*runiform())
    gen byte z = runiform() < .4
    gen double w = rnormal()
    forvalues j = 1/9 {
        gen byte v`j' = runiform() < .3
    }
    gen double C = 2 + 3*runiform()
    gen double sw = rexponential(2)
    gen double fr = rgamma(2, .5)
    gen double e0 = cond(runiform() < .3, .5*runiform(), 0)
    gen double len0 = min(sw, C)
    gen double len1 = max(C - sw, 0)
    gen double lp = 0.3*age + 0.5*z + 0.2*w + 0.15*z*w + 0.1*v1
    gen double r0 = (0.6 + 0.4*(stratum == 2))*fr*exp(lp + 0.35*xb + 0.1*xb*age)
    gen double r1 = (0.6 + 0.4*(stratum == 2))*fr*exp(lp + 0.35*(xb + 1) + 0.1*(xb + 1)*age)
    gen long n0 = rpoisson(r0*len0)
    gen long n1 = rpoisson(r1*len1)
    expand n0 + n1 + 1
    bysort id: gen long k = _n
    gen double tev = .
    replace tev = runiform()*len0 if k <= n0
    replace tev = sw + runiform()*len1 if k > n0 & k <= n0 + n1
    replace tev = C if k == n0 + n1 + 1
    gen byte ev = k <= n0 + n1
    if "`round'" != "" replace tev = ceil(4*tev)/4
    gsort id tev -ev
    by id: gen double tstart = cond(_n == 1, 0, tev[_n-1])
    drop if tev <= tstart
    stset tev, id(id) failure(ev) time0(tstart) enter(time e0) exit(time .)
    stsplit blk, after(time = sw) at(0)
    gen byte x = xb + (_t0 >= sw)
    forvalues j = 0/3 {
        gen byte d`j' = (x == `j')
    }
    gen int grp = 1 + mod(id, 20)
end

capture mata: mata drop _csh_*()
mata:
// Stata matrix out = rows r, columns c of Stata matrix in
void _csh_sel(string scalar outn, string scalar inn, real rowvector r, real rowvector c)
{
    real matrix A
    A = st_matrix(inn)
    st_matrix(outn, A[r, c])
}
// largest relative difference |a-b|/(|b|+1); both missing counts as equal
real scalar _csh_mrd(string scalar an, string scalar bn)
{
    real matrix A, B, ma, mb
    A = st_matrix(an)
    B = st_matrix(bn)
    if (rows(A) != rows(B) | cols(A) != cols(B)) return(1e300)
    ma = (A :>= .)
    mb = (B :>= .)
    if (any(ma :!= mb)) return(1e300)
    A = editmissing(A, 0)
    B = editmissing(B, 0)
    return(max(abs(A - B) :/ (abs(B) :+ 1)))
}
// scaled difference of two covariance matrices: max|A-B| / max|B|
real scalar _csh_vrd(string scalar an, string scalar bn)
{
    real matrix A, B
    A = st_matrix(an)
    B = st_matrix(bn)
    if (rows(A) != rows(B) | cols(A) != cols(B)) return(1e300)
    if (hasmissing(A) | hasmissing(B)) return(1e300)
    return(max(abs(A - B))/max(abs(B)))
}
// Shapley values by brute force over all K! orders; P: rows = masks 1..2^K-1 (bit k-1 = factor k),
// columns = horizons; result K x horizons
void _csh_perm(string scalar pn, real scalar K, string scalar outn)
{
    real matrix P, O
    real rowvector dig, used
    real scalar nseq, q, i, v, ok, m, mprev, nperm, k
    P = st_matrix(pn)
    O = J(K, cols(P), 0)
    nseq = K^K
    nperm = 0
    for (q = 0; q < nseq; q++) {
        dig = J(1, K, 0)
        v = q
        for (i = 1; i <= K; i++) {
            dig[i] = mod(v, K) + 1
            v = floor(v/K)
        }
        used = J(1, K, 0)
        ok = 1
        for (i = 1; i <= K; i++) {
            if (used[dig[i]]) ok = 0
            used[dig[i]] = 1
        }
        if (!ok) continue
        nperm++
        m = 0
        for (i = 1; i <= K; i++) {
            k = dig[i]
            mprev = m
            m = m + 2^(k-1)
            O[k, .] = O[k, .] + P[m, .] - (mprev > 0 ? P[mprev, .] : J(1, cols(P), 0))
        }
    }
    st_matrix(outn, O/nperm)
}
// 20a: replicate PAFs PJ (G x 6: joint h1 h2, z h1 h2, w h1 h2) -> rows of the combined() table
// (joint, Shapley z, Shapley w, {z}, {w}; 2 horizons) and their jackknife covariance
void _csh_jkhand(string scalar pn, string scalar vn, string scalar rn)
{
    real matrix P, R
    real rowvector pm
    real scalar G
    P = st_matrix(pn)
    G = rows(P)
    R = (P[., 1..2], (P[., 3..4] + P[., 1..2] - P[., 5..6])/2,
        (P[., 5..6] + P[., 1..2] - P[., 3..4])/2, P[., 3..4], P[., 5..6])
    pm = mean(R)
    st_matrix(vn, (G - 1)/G*quadcrossdev(R, pm, R, pm))
    st_matrix(rn, R)
}
// Shapley block of V = A V A' (A from W, subset rows: scenario 1 and K+2...); scaled difference
real scalar _csh_shapv(string scalar vn, string scalar wn, real scalar K, real scalar nh)
{
    real matrix V, W, A, VS, VA
    real colvector idx
    real scalar ns, s, h, k, R
    V = st_matrix(vn)
    W = st_matrix(wn)
    ns = cols(W)
    R = rows(V)
    A = J(K*nh, R, 0)
    for (k = 1; k <= K; k++) {
        for (h = 1; h <= nh; h++) {
            for (s = 1; s <= ns; s++) {
                A[(k-1)*nh + h, (s == 1 ? h : (K + s - 1)*nh + h)] = W[k, s]
            }
        }
    }
    idx = (nh + 1)::((K + 1)*nh)
    VS = V[idx, idx]
    VA = A*V*A'
    return(max(abs(VS - VA))/max(abs(VS)))
}
// Shapley PAFs from r(W_shap) and the subset rows of the table: scaled difference
real scalar _csh_shapp(string scalar tn, string scalar wn, real scalar K, real scalar nh)
{
    real matrix T, W
    real scalar ns, s, h, k, d, b
    T = st_matrix(tn)
    W = st_matrix(wn)
    ns = cols(W)
    d = 0
    for (k = 1; k <= K; k++) {
        for (h = 1; h <= nh; h++) {
            b = 0
            for (s = 1; s <= ns; s++) b = b + W[k, s]*T[(s == 1 ? h : (K + s - 1)*nh + h), 3]
            d = max((d, abs(b - T[k*nh + h, 3])))
        }
    }
    return(d)
}
// diag(V) vs se^2 (column 4 of the table): max relative difference of the SEs
real scalar _csh_diag(string scalar vn, string scalar tn)
{
    real matrix V, T
    V = st_matrix(vn)
    T = st_matrix(tn)
    return(max(abs(sqrt(diagonal(V)) - T[., 4]) :/ (T[., 4] :+ 1e-300)))
}
// IF vs jackknife covariance: max |var ratio - 1| and max |correlation difference|
void _csh_vcmp(string scalar an, string scalar bn)
{
    real matrix A, B, CA, CB
    real colvector da, db
    A = st_matrix(an)
    B = st_matrix(bn)
    da = sqrt(diagonal(A))
    db = sqrt(diagonal(B))
    CA = A :/ (da*da')
    CB = B :/ (db*db')
    st_numscalar("VRAT", max(abs((da:^2) :/ (db:^2) :- 1)))
    st_numscalar("CDIF", max(abs(CA - CB)))
}
// varlevels(): row and diagonal of the joint scenario = sums over the level scenarios
real scalar _csh_lvl(string scalar vn, real scalar nsc, real scalar nh)
{
    real matrix V
    real colvector L
    real scalar h, s, d, rj
    V = st_matrix(vn)
    d = 0
    for (h = 1; h <= nh; h++) {
        rj = h
        L = J(0, 1, .)
        for (s = 2; s <= nsc; s++) L = L \ ((s - 1)*nh + h)
        d = max((d, max(abs(V[rj, .] - colsum(V[L, .])))/max(abs(V))))
        d = max((d, abs(V[rj, rj] - sum(V[L, L]))/V[rj, rj]))
    }
    return(d)
}
// mi by hand: P (M x k) PAFs, S (M x k) within SEs on the PAF scale, VS (k x k) sum of the
// within covariances; Shapley rows (scenarios 2..K+1, flag f) on the PAF scale, the other rows
// on log(1-PAF) with within SEs SL; returns table columns paf se (k x 2) and V
void _csh_mihand(string scalar pn, string scalar sn, string scalar sln, string scalar vsn,
    string scalar fn, string scalar outn, string scalar von)
{
    real matrix P, S, SL, VS, O
    real colvector f
    real scalar M, k, c, qb, w, b, t
    P = st_matrix(pn)
    S = st_matrix(sn)
    SL = st_matrix(sln)
    VS = st_matrix(vsn)
    f = st_matrix(fn)
    M = rows(P)
    k = cols(P)
    O = J(k, 2, .)
    for (c = 1; c <= k; c++) {
        if (f[c]) {
            qb = mean(P[., c])
            w = mean(S[., c]:^2)
            b = variance(P[., c])
            O[c, .] = (qb, sqrt(w + (1 + 1/M)*b))
        }
        else {
            qb = mean(ln(1 :- P[., c]))
            w = mean(SL[., c]:^2)
            b = variance(ln(1 :- P[., c]))
            t = w + (1 + 1/M)*b
            O[c, .] = (1 - exp(qb), exp(qb)*sqrt(t))
        }
    }
    st_matrix(outn, O)
    st_matrix(von, VS/M + (1 + 1/M)*variance(P))
}
// Barnard-Rubin df (as mi estimate) of column c of Q with variances U, complete df nuc
real scalar _csh_br(string scalar qn, string scalar un, real scalar c, real scalar nuc)
{
    real matrix Q, U
    real scalar M, ub, b, t, r, nul, gam, nuo
    Q = st_matrix(qn)
    U = st_matrix(un)
    M = rows(Q)
    ub = mean(U[., c])
    b = variance(Q[., c])
    t = ub + (1 + 1/M)*b
    r = (1 + 1/M)*b/ub
    nul = (M - 1)*(1 + 1/r)^2
    gam = (1 + 1/M)*b/t
    nuo = nuc*(nuc + 1)*(1 - gam)/(nuc + 3)
    return(1/(1/nul + 1/nuo))
}
// PAF and E_cf by hand (Breslow, strata, (t0, t1] records): observed columns xl, coefficient
// row vector HB (same order), counterfactual column lists separated by |; horizons 1 and 2.5;
// result HAND (scenarios*2 x 4: scenario, horizon, paf, E_cf)
void _csh_hand_m(string scalar xl, string scalar scl, string scalar es)
{
    real matrix X, AC, H
    string rowvector sl
    real colvector t0, t1, dd, st, a, ev, risk
    real rowvector b, E0, tau
    real scalar ns, s, ti, i, r, tt, S0, E
    b = st_matrix("HB")
    X = st_data(., tokens(xl), es)
    sl = tokens(scl, "|")
    sl = select(sl, sl :!= "|")
    ns = cols(sl)
    t0 = st_data(., "_t0", es)
    t1 = st_data(., "_t", es)
    dd = st_data(., "_d", es)
    st = st_data(., "stratum", es)
    a = exp(X*b')
    AC = J(rows(X), ns, .)
    for (s = 1; s <= ns; s++) AC[., s] = exp(st_data(., tokens(sl[s]), es)*b')
    ev = selectindex(dd :== 1)
    H = J(2*ns, 4, .)
    tau = (1, 2.5)
    for (ti = 1; ti <= 2; ti++) {
        E = 0
        E0 = J(1, ns, 0)
        for (i = 1; i <= rows(ev); i++) {
            r = ev[i]
            if (t1[r] > tau[ti]) continue
            tt = t1[r]
            risk = selectindex((st :== st[r]) :& (t0 :< tt) :& (t1 :>= tt))
            S0 = sum(a[risk])
            E = E + 1
            E0 = E0 + colsum(AC[risk, .])/S0
        }
        for (s = 1; s <= ns; s++) H[(s - 1)*2 + ti, .] = (s, tau[ti], 1 - E0[s]/E, E0[s])
    }
    st_matrix("HAND", H)
}
end

capture program drop _csh_cmp
program define _csh_cmp
    args a b tol what
    mata: st_numscalar("CSHD", _csh_mrd("`a'", "`b'"))
    di as txt "   `what': max rel. difference " %9.2e scalar(CSHD)
    assert scalar(CSHD) < `tol'
end
capture program drop _csh_cmpv
program define _csh_cmpv
    args a b tol what
    mata: st_numscalar("CSHD", _csh_vrd("`a'", "`b'"))
    di as txt "   `what': max scaled difference " %9.2e scalar(CSHD)
    assert scalar(CSHD) < `tol'
end
* additivity: sum of the Shapley rows (scenarios 2..K+1) = joint row, each horizon
capture program drop _csh_add
program define _csh_add
    args T nh K what
    local dmax 0
    forvalues h = 1/`nh' {
        scalar SP = 0
        forvalues s = 2/`=`K'+1' {
            scalar SP = SP + `T'[(`s' - 1)*`nh' + `h', 3]
        }
        local dmax = max(`dmax', abs(SP - `T'[`h', 3]))
    }
    di as txt "   `what': sum of Shapley PAFs - joint PAF, max " %9.2e `dmax'
    assert `dmax' < 1e-12
end

local coxA "stcox i.x##c.age c.z##c.w, strata(stratum) vce(cluster id) breslow"
local coxB "stcox d1 d2 d3 age c.d1#c.age c.d2#c.age c.d3#c.age c.z##c.w, strata(stratum) vce(cluster id) breslow"
local f3 "set(x=0) eliminate(z) cap(w=0)"

* =========================================================================== 20a
_csh_sim, n(400)
qui `coxA'
stpaf, combined(eliminate(z) cap(w=0)) subsets horizon(1 2.5) cluster(grp) refit nodots
assert r(K) == 2 & "`r(factors)'" == "z w" & "`r(jk_method)'" == "refit" & r(df) == 19
assert "`r(scen1)'" == "joint" & "`r(scen2)'" == "Shapley: z" & "`r(scen3)'" == "Shapley: w"
assert "`r(scen4)'" == "z" & "`r(scen5)'" == "w"
matrix TR = r(table)
matrix VR = r(V_paf)
assert rowsof(TR) == 10 & rowsof(VR) == 10
* explicit loop: stcox refit without cluster g, subset PAFs from separate calls
matrix PJ = J(20, 6, .)
forvalues g = 1/20 {
    qui stcox i.x##c.age c.z##c.w if grp != `g', strata(stratum) breslow
    qui stpaf, eliminate(z) cap(w=0) each horizon(1 2.5) vce(if)
    matrix T = r(table)
    forvalues r = 1/6 {
        matrix PJ[`g', `r'] = T[`r', 3]
    }
}
mata: _csh_jkhand("PJ", "VH", "RH")
_csh_cmpv VR VH 1e-8 "20a refit r(V_paf) vs hand jackknife covariance"
* point estimates: Shapley = average of the two orders (full data)
qui `coxA'
qui stpaf, eliminate(z) cap(w=0) each horizon(1 2.5) vce(if)
matrix T = r(table)
forvalues h = 1/2 {
    local pj = T[`h', 3]
    local pz = T[2 + `h', 3]
    local pw = T[4 + `h', 3]
    assert abs(TR[2 + `h', 3] - (`pz' + `pj' - `pw')/2) < 1e-12
    assert abs(TR[4 + `h', 3] - (`pw' + `pj' - `pz')/2) < 1e-12
    * Shapley rows: SE from the hand covariance, t(19) CI on the PAF scale, E_cf/se_logR missing
    forvalues k = 1/2 {
        local r = 2*`k' + `h'
        assert reldif(TR[`r', 4], sqrt(VH[`r', `r'])) < 1e-6
        assert reldif(TR[`r', 5], TR[`r', 3] - invttail(19, .025)*TR[`r', 4]) < 1e-10
        assert reldif(TR[`r', 6], TR[`r', 3] + invttail(19, .025)*TR[`r', 4]) < 1e-10
        assert TR[`r', 8] >= . & TR[`r', 9] >= .
        assert reldif(TR[`r', 7], TR[`h', 7]) < 1e-12
    }
}
_csh_add TR 2 2 "20a refit"
* Newton jackknife
qui stpaf, combined(eliminate(z) cap(w=0)) horizon(1 2.5) cluster(grp)
assert "`r(jk_method)'" == "newton"
matrix TN = r(table)
matrix VN = r(V_paf)
_csh_cmp TN TR 1e-5 "20a Newton vs refit, r(table)"
_csh_cmpv VN VH 1e-4 "20a Newton r(V_paf) vs hand jackknife covariance"
_csh_add TN 2 2 "20a Newton"
di as txt "[ok] 20a: K = 2, jackknife covariance and Shapley PAFs = explicit refit loop"

* =========================================================================== 20b
qui `coxA'
stpaf, combined(`f3') subsets horizon(1 2.5) vce(if)
assert r(K) == 3 & "`r(factors)'" == "x z w" & "`r(vce)'" == "if"
assert "`r(scen1)'" == "joint" & "`r(scen2)'" == "Shapley: x" & "`r(scen4)'" == "Shapley: w"
assert "`r(scen5)'" == "x" & "`r(scen7)'" == "w" & "`r(scen8)'" == "x,z" & "`r(scen10)'" == "z,w"
assert r(sumdiff) < 1e-12
matrix T3 = r(table)
matrix V3 = r(V_paf)
matrix W3 = r(W_shap)
matrix S3 = r(subsets)
assert rowsof(T3) == 20 & colsof(T3) == 10 & rowsof(V3) == 20
assert rowsof(W3) == 3 & colsof(W3) == 7 & rowsof(S3) == 7 & colsof(S3) == 3
assert S3[1, 1] == 1 & S3[1, 2] == 1 & S3[1, 3] == 1 & S3[2, 1] == 1 & S3[2, 2] == 0 & S3[7, 1] == 0
assert reldif(W3[1, 1], 1/3) < 1e-15 & reldif(W3[1, 2], 1/3) < 1e-15 & reldif(W3[1, 3], -1/6) < 1e-15
* subset rows = separate calls; PM: subset PAFs by mask (x = 1, z = 2, w = 4), 2 horizons
matrix PM = J(7, 2, .)
local opt1 "set(x=0)"
local opt2 "eliminate(z)"
local opt4 "cap(w=0)"
local opt3 "set(x=0) eliminate(z)"
local opt5 "set(x=0) cap(w=0)"
local opt6 "eliminate(z) cap(w=0)"
local opt7 "`f3'"
* scenario of each mask in the combined() table
local sc1 5
local sc2 6
local sc4 7
local sc3 8
local sc5 9
local sc6 10
local sc7 1
foreach m in 1 2 4 3 5 6 7 {
    qui stpaf, `opt`m'' horizon(1 2.5) vce(if)
    matrix D = r(table)
    matrix D = D[1..2, 2...]
    local r0 = (`sc`m'' - 1)*2 + 1
    matrix C = T3[`r0'..`r0'+1, 2...]
    _csh_cmp C D 1e-10 "20b subset `opt`m''"
    matrix PM[`m', 1] = D[1, 2]
    matrix PM[`m', 2] = D[2, 2]
}
* Shapley PAFs = average over the 3! orders (brute force)
mata: _csh_perm("PM", 3, "SH")
matrix C = (T3[3..4, 3]', T3[5..6, 3]', T3[7..8, 3]')
matrix SHv = (SH[1, 1..2], SH[2, 1..2], SH[3, 1..2])
_csh_cmp C SHv 1e-12 "20b Shapley PAFs vs average over the 3! orders"
_csh_add T3 2 3 "20b IF"
mata: st_numscalar("CSHD", _csh_shapp("T3", "W3", 3, 2))
di as txt "   20b Shapley = W_shap x subset PAFs, max difference " %9.2e scalar(CSHD)
assert scalar(CSHD) < 1e-12
mata: st_numscalar("CSHD", _csh_diag("V3", "T3"))
di as txt "   20b sqrt(diag V_paf) vs se, max rel. difference " %9.2e scalar(CSHD)
assert scalar(CSHD) < 1e-10
mata: st_numscalar("CSHD", _csh_shapv("V3", "W3", 3, 2))
di as txt "   20b Shapley block of V_paf vs W V_sub W', scaled " %9.2e scalar(CSHD)
assert scalar(CSHD) < 1e-10
* r(V_paf) of an each call = matching block of the combined() covariance
qui stpaf, `f3' each horizon(1 2.5) vce(if)
matrix VE = r(V_paf)
matrix rownames VE = r1 r2 r3 r4 r5 r6 r7 r8
matrix colnames VE = r1 r2 r3 r4 r5 r6 r7 r8
* (each: scenarios joint, z, x, w -- eliminate() changes come before set(); combined(): {x} = 5,
*  {z} = 6, {w} = 7)
mata: _csh_sel("VC", "V3", (1, 2, 11, 12, 9, 10, 13, 14), (1, 2, 11, 12, 9, 10, 13, 14))
_csh_cmpv VC VE 1e-10 "20b each-call V_paf = block of combined() V_paf"
* display without subsets: same r(table)
qui stpaf, combined(`f3') horizon(1 2.5) vce(if)
matrix T = r(table)
_csh_cmp T T3 1e-15 "20b without option subsets"
di as txt "[ok] 20b: K = 3 Shapley PAFs, subset rows, additivity, covariance (IF)"

* =========================================================================== 20c
qui `coxB'
stpaf, combined((eliminate(d1) eliminate(d2 d3)) eliminate(z) cap(w=0)) horizon(1 2.5) vce(if)
assert r(K) == 3 & "`r(factors)'" == "d1+d2+d3 z w"
matrix TB = r(table)
matrix VB = r(V_paf)
_csh_cmp TB T3 1e-6 "20c grouped indicators vs i.x set(x=0), IF table"
_csh_cmpv VB V3 1e-6 "20c grouped indicators vs i.x, IF V_paf"
qui stpaf, combined( ( elim(d1 d2)  set(d3 = 0) )  eliminate(z)  cap( w = 0 ) ) horizon(1 2.5) vce(if)
assert "`r(factors)'" == "d1+d2+d3 z w"
matrix T = r(table)
_csh_cmp T TB 1e-12 "20c blanks, abbreviation, mixed group"
qui stpaf, combined((eliminate(d1) eliminate(d2 d3)) eliminate(z) cap(w=0)) horizon(1 2.5) cluster(grp)
assert "`r(vce)'" == "jackknife"
matrix TBj = r(table)
matrix VBj = r(V_paf)
qui `coxA'
qui stpaf, combined(`f3') horizon(1 2.5) cluster(grp)
matrix T = r(table)
matrix V = r(V_paf)
_csh_cmp TBj T 1e-6 "20c grouped indicators vs i.x, jackknife table"
_csh_cmpv VBj V 1e-5 "20c grouped indicators vs i.x, jackknife V_paf"
_csh_add T 2 3 "20c jackknife"
* one factor per variable outside parentheses
qui stpaf, combined(eliminate(z w)) horizon(1 2.5) vce(if)
assert r(K) == 2 & "`r(factors)'" == "z w"
matrix T1 = r(table)
qui stpaf, combined(eliminate(z) set(w=0)) horizon(1 2.5) vce(if)
matrix T2 = r(table)
_csh_cmp T1 T2 1e-14 "20c combined(eliminate(z w)) = combined(eliminate(z) set(w=0))"
di as txt "[ok] 20c: grouped factors = factor variable; parsing"

* =========================================================================== 20d
preserve
_csh_sim, n(400) round
local coxAe "stcox i.x##c.age c.z##c.w, strata(stratum) vce(cluster id) efron"
qui `coxAe'
stpaf if age > -1, combined(`f3') horizon(1 2.5) vce(if)
assert "`r(if_ties)'" == "efron" & r(K) == 3
matrix TI = r(table)
matrix VI = r(V_paf)
_csh_add TI 2 3 "20d efron IF, if age > -1"
forvalues r = 3/8 {
    assert reldif(TI[`r', 5], TI[`r', 3] - invnormal(.975)*TI[`r', 4]) < 1e-10
    assert reldif(TI[`r', 6], TI[`r', 3] + invnormal(.975)*TI[`r', 4]) < 1e-10
    assert TI[`r', 8] >= . & TI[`r', 9] >= .
}
mata: st_numscalar("CSHD", _csh_diag("VI", "TI"))
assert scalar(CSHD) < 1e-10
qui stpaf if age > -1, combined(`f3') horizon(1 2.5) vce(jackknife)
assert "`r(jk_method)'" == "newton" & r(N_clust) >= 390
local G = r(N_clust)
matrix TJ = r(table)
matrix VJ = r(V_paf)
_csh_add TJ 2 3 "20d efron Newton jackknife (all clusters)"
forvalues r = 3/8 {
    assert reldif(TJ[`r', 5], TJ[`r', 3] - invttail(`G' - 1, .025)*TJ[`r', 4]) < 1e-10
}
mata: _csh_vcmp("VI", "VJ")
di as txt "   20d IF vs jackknife covariance: max |var ratio - 1| " %6.3f scalar(VRAT) ///
    ", max |correlation difference| " %6.3f scalar(CDIF)
assert scalar(VRAT) < .15 & scalar(CDIF) < .05
* Newton = refit (20 clusters)
qui stpaf if age > -1, combined(`f3') horizon(1 2.5) cluster(grp)
matrix TN = r(table)
matrix VN = r(V_paf)
qui stpaf if age > -1, combined(`f3') horizon(1 2.5) cluster(grp) refit nodots
matrix TR = r(table)
matrix VR = r(V_paf)
_csh_cmp TN TR 1e-5 "20d Efron Newton vs refit, table"
_csh_cmpv VN VR 1e-4 "20d Efron Newton vs refit, V_paf"
restore
di as txt "[ok] 20d: efron ties, if subset, IF and jackknife"

* =========================================================================== 20e
qui `coxA'
qui stpaf, varlevels(x) horizon(1 2.5) vce(if)
matrix V = r(V_paf)
matrix T = r(table)
mata: st_numscalar("CSHD", _csh_lvl("V", 4, 2))
di as txt "   20e varlevels IF: joint row of V_paf = sum of level rows, scaled " %9.2e scalar(CSHD)
assert scalar(CSHD) < 1e-10
mata: st_numscalar("CSHD", _csh_diag("V", "T"))
assert scalar(CSHD) < 1e-10
qui stpaf, varlevels(x) horizon(1 2.5) cluster(grp)
matrix V = r(V_paf)
mata: st_numscalar("CSHD", _csh_lvl("V", 4, 2))
di as txt "   20e varlevels jackknife: joint row = sum of level rows, scaled " %9.2e scalar(CSHD)
assert scalar(CSHD) < 1e-10
qui stpaf, set(x=0) vce(if)
matrix V = r(V_paf)
matrix T = r(table)
assert rowsof(V) == 1 & reldif(sqrt(V[1, 1]), T[1, 4]) < 1e-10
* fweights without clusters = expanded data
preserve
clear
set obs 300
gen byte z = runiform() < .5
gen double w = rnormal()
gen double t = rexponential(1)/exp(.5*z + .3*w)
gen byte d = t < 2
replace t = min(t, 2)
gen int fw = 1 + floor(3*runiform())
stset t [fw=fw], failure(d)
qui stcox z w, breslow
qui stpaf, eliminate(z) set(w=0) each horizon(.5 1) vce(if)
matrix T1 = r(table)
matrix V1 = r(V_paf)
expand fw
stset t, failure(d)
qui stcox z w, breslow
qui stpaf, eliminate(z) set(w=0) each horizon(.5 1) vce(if)
matrix T2 = r(table)
matrix V2 = r(V_paf)
_csh_cmp T1 T2 1e-8 "20e fweights vs expanded data, table"
_csh_cmpv V1 V2 1e-8 "20e fweights vs expanded data, V_paf"
restore
di as txt "[ok] 20e: r(V_paf) with varlevels(), one scenario, fweights"

* =========================================================================== 20f
preserve
_csh_sim, n(400)
keep if _st == 1
gen double t0b = _t0
gen double t1b = _t
gen byte db = _d
stset, clear
keep id t0b t1b db x d0 d1 d2 d3 age z w stratum grp
gen double u = runiform()
bysort id (t0b): replace u = u[1]
replace w = . if u < .3
drop u
mi set flong
mi register imputed w
mi register regular x d0 d1 d2 d3 age z stratum db t0b t1b grp id
mi impute regress w i.x z age i.stratum db, add(5) rseed(18)
mi stset t1b, id(id) failure(db) time0(t0b) exit(time .)
capture erase f20f.ster
qui mi estimate, saving(f20f, replace): `coxA'
stpaf, combined(`f3') horizon(1 2.5) vce(if)
assert r(M) == 5 & r(mi_nsaved) == 5 & r(K) == 3 & "`r(factors)'" == "x z w"
matrix TA = r(table)
matrix PA = r(paf_mi)
matrix SA = r(se_mi)
matrix VA = r(V_paf)
local sdl = r(sumdiff)
assert colsof(TA) == 13 & rowsof(TA) == 20 & rowsof(VA) == 20
* explicit loop over the imputations (data reloaded from a tempfile: one preserve per level)
matrix PL = J(5, 20, .)
matrix SL = J(5, 20, .)
matrix SR = J(5, 20, .)
tempfile mid
qui save `mid'
forvalues m = 1/5 {
    qui use `mid', clear
    qui mi extract `m', clear
    qui `coxA'
    qui stpaf, combined(`f3') horizon(1 2.5) vce(if)
    matrix T = r(table)
    if `m' == 1 matrix VS = r(V_paf)
    else matrix VS = VS + r(V_paf)
    forvalues c = 1/20 {
        matrix PL[`m', `c'] = T[`c', 3]
        matrix SL[`m', `c'] = T[`c', 4]
        matrix SR[`m', `c'] = T[`c', 9]
    }
}
qui use `mid', clear
_csh_cmp PA PL 1e-10 "20f r(paf_mi) = loop"
matrix FL = J(20, 1, 0)
forvalues r = 3/8 {
    matrix FL[`r', 1] = 1
}
mata: _csh_mihand("PL", "SL", "SR", "VS", "FL", "OH", "VHm")
matrix C = TA[1..., 3..4]
_csh_cmp C OH 1e-10 "20f hand Rubin (Shapley rows PAF scale, others log scale), paf and se"
_csh_cmpv VA VHm 1e-10 "20f r(V_paf) = mean within + (1+1/M) between covariance"
* r(se_mi): PAF-scale SEs for the Shapley rows, se_logR for the others
mata: _csh_sel("SEX1", "SR", (1..5), (1..2))
mata: _csh_sel("SEX2", "SL", (1..5), (3..8))
mata: _csh_sel("SEX3", "SR", (1..5), (9..20))
matrix SEX = SEX1, SEX2, SEX3
_csh_cmp SA SEX 1e-10 "20f r(se_mi)"
* per-imputation additivity
scalar DM = 0
forvalues m = 1/5 {
    forvalues h = 1/2 {
        scalar SP = PA[`m', 2 + `h'] + PA[`m', 4 + `h'] + PA[`m', 6 + `h']
        scalar DM = max(DM, abs(SP - PA[`m', `h']))
    }
}
di as txt "   20f per-imputation additivity, max difference " %9.2e scalar(DM)
assert scalar(DM) < 1e-12
* citrans(log): r(sumdiff) = difference of the combined estimates
scalar DL = 0
forvalues h = 1/2 {
    scalar SP = TA[2 + `h', 3] + TA[4 + `h', 3] + TA[6 + `h', 3]
    scalar DL = max(DL, abs(SP - TA[`h', 3]))
}
di as txt "   20f citrans(log): sum of Shapley PAFs - joint PAF " %9.2e scalar(DL) " (r(sumdiff) " %9.2e `sdl' ")"
assert reldif(`sdl', scalar(DL)) < 1e-6
* citrans(none): exact, and diag(V_paf) = se^2 for all rows
qui mi estimate using f20f
qui stpaf, combined(`f3') horizon(1 2.5) vce(if) citrans(none)
assert r(mi_nsaved) == 5 & r(sumdiff) < 1e-12
matrix T = r(table)
matrix V = r(V_paf)
_csh_add T 2 3 "20f citrans(none)"
mata: st_numscalar("CSHD", _csh_diag("V", "T"))
assert scalar(CSHD) < 1e-10
* jackknife (20 clusters): Barnard-Rubin df of a Shapley row by hand
qui mi estimate, saving(f20f, replace): `coxA'
qui stpaf, combined(`f3') horizon(1 2.5) cluster(grp) citrans(none)
assert "`r(vce)'" == "jackknife" & r(df_c) == 19 & r(sumdiff) < 1e-12
matrix T = r(table)
matrix P = r(paf_mi)
matrix S = r(se_mi)
mata: st_matrix("U", st_matrix("S"):^2)
mata: st_numscalar("BRDF", _csh_br("P", "U", 3, 19))
di as txt "   20f jackknife: df of Shapley row (x, h1) " %8.3f T[3, 11] ", hand Barnard-Rubin " %8.3f scalar(BRDF)
assert reldif(T[3, 11], scalar(BRDF)) < 1e-8
capture erase f20f.ster
restore
di as txt "[ok] 20f: combined() after mi estimate"

* =========================================================================== 20g
qui `coxA'
capture stpaf, combined(eliminate(z) cap(w=0)) set(x=0)
assert _rc == 198
capture stpaf, combined(eliminate(z) cap(w=0)) each
assert _rc == 198
capture stpaf, combined(eliminate(z) cap(w=0)) varlevels(x)
assert _rc == 198
capture stpaf, combined(eliminate(z))
assert _rc == 198
capture stpaf, combined((eliminate(z) cap(w=0)))
assert _rc == 198
capture stpaf, combined(eliminate(z) set(z=1))
assert _rc == 198
capture stpaf, combined(((eliminate(z))) cap(w=0))
assert _rc == 198
capture stpaf, combined(remove(z) cap(w=0))
assert _rc == 198
capture stpaf, combined(eliminate(nosuchvar) cap(w=0))
assert _rc == 111
capture stpaf, combined(set(x=a) eliminate(z))
assert _rc != 0
capture stpaf, combined(eliminate(v1-v9))
assert _rc == 198
qui stpaf, combined(set(x=0 z=0)) vce(if)
assert r(K) == 2 & "`r(factors)'" == "x z"
di as txt "[ok] 20g: errors"

* =========================================================================== 20h
* continuous variable in an interaction with a factor variable (c.w##i.z): PAFs by hand
qui stcox c.w##i.z age, strata(stratum) breslow
matrix HB = (_b[w], _b[1.z], _b[1.z#c.w], _b[age])
gen double hz = (z == 1)
gen double hzw = (z == 1)*w
gen double wc = min(w, 0)
gen double hzwc = (z == 1)*wc
gen double zero = 0
tempvar es
gen byte `es' = e(sample)
* scenarios: 1 both (w capped, z = 0), 2 w capped, 3 z = 0
mata: _csh_hand_m("w hz hzw age", "wc zero zero age|wc hz hzwc age|w zero zero age", "`es'")
qui stpaf, combined(cap(w=0) set(z=0)) subsets horizon(1 2.5) vce(if)
matrix T = r(table)
* rows: 1-2 joint, 3-4 Shapley w, 5-6 Shapley z, 7-8 {w}, 9-10 {z}
matrix C = (T[1..2, 3], T[1..2, 8] \ T[7..8, 3], T[7..8, 8] \ T[9..10, 3], T[9..10, 8])
matrix H = HAND[1..., 3..4]
_csh_cmp C H 1e-10 "20h c.w##i.z: subset PAFs and E_cf vs hand"
qui stpaf, cap(w=0) horizon(1 2.5) vce(if)
matrix T1 = r(table)
matrix C = T1[1..., 3], T1[1..., 8]
matrix H = HAND[3..4, 3..4]
_csh_cmp C H 1e-10 "20h cap(w=0) alone vs hand"
forvalues h = 1/2 {
    local pj = HAND[`h', 3]
    local pw = HAND[2 + `h', 3]
    local pz = HAND[4 + `h', 3]
    assert abs(T[2 + `h', 3] - (`pw' + `pj' - `pz')/2) < 1e-10
    assert abs(T[4 + `h', 3] - (`pz' + `pj' - `pw')/2) < 1e-10
}
drop hz hzw wc hzwc zero
di as txt "[ok] 20h: c.w##i.z by hand; Shapley PAFs from hand subset PAFs"

* =========================================================================== 20i
qui stcox i.x##c.age c.z##c.w v1-v8, strata(stratum) vce(cluster id) breslow
timer clear
timer on 1
qui stpaf, combined(`f3' eliminate(v1 v2 v3)) vce(if)
timer off 1
assert r(K) == 6 & r(sumdiff) < 1e-12 & rowsof(r(table)) == 69
timer on 2
qui stpaf, combined(`f3' eliminate(v1 v2 v3)) cluster(grp)
timer off 2
assert "`r(vce)'" == "jackknife" & r(K) == 6 & r(sumdiff) < 1e-12
timer on 3
qui stpaf, combined(`f3' eliminate(v1-v5)) vce(if)
timer off 3
assert r(K) == 8 & r(sumdiff) < 1e-12 & rowsof(r(table)) == 263 & rowsof(r(V_paf)) == 263
timer on 4
qui stpaf, combined(`f3' eliminate(v1-v5)) cluster(grp)
timer off 4
assert "`r(vce)'" == "jackknife" & r(K) == 8 & r(sumdiff) < 1e-12
qui timer list
di as txt "   20i timings (s): K=6 IF " %6.2f r(t1) ", jackknife " %6.2f r(t2) ///
    "; K=8 IF " %6.2f r(t3) ", jackknife " %6.2f r(t4)
di as txt "[ok] 20i: K = 6 and 8"

di as res _n "stpaf combined() certification finished (20a-i)"
