"""Tradecraft transpositions x periodic substitution, EXACT crib tests.
Families (all perms m[pt]=ct, both directions: perm and its inverse):
  MY  Myszkowski: all weak orderings (rank vectors with ties) of width w; ct read = for each rank,
      row by row, left to right across the tied columns. (Distinct ranks = plain columnar.)
  DV  VIC disrupted columnar. Variant 'gap' (VIC/SECOM): triangle t starts in the column read t-th,
      occupies that col..right edge, shrinking one col per row; after its single last-col cell one
      full row, then triangle t+1. Variant 'nogap': next triangle starts on the very next row.
      Fill: non-disrupted cells row-major (first N cells of the grid), then disrupted cells row-major;
      read columns in key order.
  IR  Irregular columnar: 'right' = short last row right-justified (long columns are the rightmost r);
      'keylong' = the long columns are the first r columns in key (read) order.
Scorer: run3 = new_doublecol.run2 extended: periodic key ST/KA Vig/Beau periods 1..24 (pt- or ct-indexed,
redundancy>=6), arbitrary alphabet periods 1..12 (>=4 strong equalities). Cross-checked vs gen.run.
usage: python3 new_tradecraft.py [budget_w9_mysz_fraction] [w9/w10 disrupted samples]
"""
import sys, time, itertools, random, math
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, '/tmp/claude-0/-home-claude/bd15c051-bac5-5123-8757-e9921faf2557/scratchpad')
from stack import CT0, ANCH, KA
from gen import run, ST

N = 97; PKEY = 24; PANY = 12; CHUNK = 20000
AN = np.array(sorted(ANCH))
PR = np.array(list(itertools.combinations(range(len(AN)), 2))); P0, P1 = PR[:, 0], PR[:, 1]
PT_PAIRS = {P: np.nonzero((AN % P)[P0] == (AN % P)[P1])[0] for P in range(1, PKEY + 1)}
PT_RED = {P: len(AN) - len(set((AN % P).tolist())) for P in range(1, PKEY + 1)}
Pl = np.array([ord(ANCH[i]) for i in AN]); Peq = Pl[P0] == Pl[P1]

DIV = np.array([sum(1 << (P - 1) for P in range(1, PKEY + 1) if d % P == 0) for d in range(N)], dtype=np.int32)
DIVPT = DIV[np.abs(AN[P0] - AN[P1])]
FULLK = (1 << PKEY) - 1; FULLA = (1 << PANY) - 1
def run3(ct, Ms, labs, RTHR=6, ANYSTRONG=4):
    """Bitmask re-implementation of new_doublecol.run2 (keyed periods 1..PKEY, any-alphabet 1..PANY).
    For each perm: invalid-period mask = OR over 'conflicting' anchor pairs of the mask of periods dividing
    the pair's index difference (pt or ct). Same hit tuples as gen.run (cross-checked in campaign())."""
    hits = []
    CI = Ms[:, AN].astype(np.int16)
    Cl = np.frombuffer(ct.encode(), dtype=np.uint8)[CI]
    B = DIV[np.abs(CI[:, P0] - CI[:, P1])]
    for an, A in (('ST', ST), ('KA', KA)):
        lut = np.zeros(128, dtype=np.int8); lut[[ord(c) for c in A]] = np.arange(26)
        C = lut[Cl]; P_ = np.array([A.index(ANCH[i]) for i in AN], dtype=np.int8)
        for rule in (0, 1):
            K = ((C - P_) % 26 if rule == 0 else (C + P_) % 26).astype(np.int8)
            D = K[:, P0] != K[:, P1]
            vpt = ~np.bitwise_or.reduce(np.where(D, DIVPT[None, :], 0), axis=1) & FULLK
            vct = ~np.bitwise_or.reduce(np.where(D, B, 0), axis=1) & FULLK
            for x in np.nonzero(vpt)[0]:
                for P in range(1, PKEY + 1):
                    if vpt[x] >> (P - 1) & 1 and PT_RED[P] >= RTHR: hits.append(('key', an, rule, 'pt', P, PT_RED[P], labs[x]))
            for x in np.nonzero(vct)[0]:
                for P in range(1, PKEY + 1):
                    if vct[x] >> (P - 1) & 1:
                        rv = (CI[x] % P).tolist(); red = len(rv) - len(set(rv))
                        if red >= RTHR: hits.append(('key', an, rule, 'ct', P, red, labs[x]))
    Ceq = Cl[:, P0] == Cl[:, P1]; bad = Peq[None, :] ^ Ceq
    vpt = ~np.bitwise_or.reduce(np.where(bad, DIVPT[None, :], 0), axis=1) & FULLA
    vct = ~np.bitwise_or.reduce(np.where(bad, B, 0), axis=1) & FULLA
    for x in np.nonzero(vpt)[0]:
        for P in range(1, PANY + 1):
            if vpt[x] >> (P - 1) & 1:
                strong = int(Peq[PT_PAIRS[P]].sum())
                if strong >= ANYSTRONG: hits.append(('any', 'pt', P, strong, labs[x]))
    for x in np.nonzero(vct)[0]:
        for P in range(1, PANY + 1):
            if vct[x] >> (P - 1) & 1:
                sm = (CI[x, P0] - CI[x, P1]) % P == 0; st = int((sm & Peq).sum())
                if st >= ANYSTRONG: hits.append(('any', 'ct', P, st, labs[x]))
    return hits

def order_to_perm(order):            # order[k] = pt index read k-th  ->  m[pt]=ct
    m = np.empty(N, dtype=np.int16); m[np.asarray(order)] = np.arange(N); return m

# ------------------------------------------------------------------ Myszkowski
def rgs(w):                          # restricted growth strings = set partitions
    out = []
    def rec(s, k):
        if len(s) == w: out.append((tuple(s), k)); return
        for b in range(k + 1): rec(s + [b], max(k, b + 1))
    rec([0], 1); return out
def weak_orders(w):
    blocks = []
    by_k = {}
    for s, k in rgs(w): by_k.setdefault(k, []).append(s)
    for k, ss in sorted(by_k.items()):
        S = np.array(ss, dtype=np.int8); Pm = np.array(list(itertools.permutations(range(k))), dtype=np.int8)
        blocks.append(Pm[:, S].reshape(-1, w))       # rank = perm[block]
    return np.concatenate(blocks)
IDX = np.arange(N)
def mysz_perms(R):                   # R: (n,w) ranks -> (n,N) m[pt]=ct
    w = R.shape[1]
    key = R[:, IDX % w].astype(np.int32) * 128 + IDX[None, :]
    order = np.argsort(key, axis=1, kind='stable')
    m = np.empty_like(order); np.put_along_axis(m, order, np.arange(N)[None, :].repeat(len(R), 0), axis=1)
    return m.astype(np.int16)
def mysz_ref(rank):                  # slow reference (validation)
    w = len(rank); rows = math.ceil(N / w); order = []
    for r in sorted(set(rank)):
        cols = [c for c in range(w) if rank[c] == r]
        for row in range(rows):
            for c in cols:
                if row * w + c < N: order.append(row * w + c)
    return order_to_perm(order)

# ------------------------------------------------------------------ VIC disrupted
def disrupted_order(order, gap, n=N):
    w = len(order); rows = math.ceil(n / w); D = np.zeros((rows, w), bool)
    s = 0
    for t in range(w):
        c = order[t]
        if s >= rows: break
        for j in range(w - c):
            if s + j < rows: D[s + j, c + j:] = True
        s += (w - c) + (1 if gap else 0)
    cells = [(r, c) for r in range(rows) for c in range(w) if r * w + c < n]
    fill = [x for x in cells if not D[x]] + [x for x in cells if D[x]]
    pos = {x: i for i, x in enumerate(fill)}            # grid cell -> pt index
    rd = [pos[(r, c)] for c in order for r in range(rows) if r * w + c < n]
    return rd, D
def disrupted_perm(order, gap): return order_to_perm(disrupted_order(order, gap)[0])

# ------------------------------------------------------------------ irregular columnar
def irregular_perm(order, var):
    w = len(order); L, r = divmod(N, w)
    if var == 'right': long_cols = set(range(w - r, w))
    else: long_cols = set(order[:r])                    # 'keylong'
    lens = [L + (c in long_cols) for c in range(w)]
    # fill row-major across the cells that exist: row L (last) has cells only in long columns
    cells = [(row, c) for row in range(L + 1) for c in range(w) if row < lens[c]]
    pos = {x: i for i, x in enumerate(cells)}
    return order_to_perm([pos[(row, c)] for c in order for row in range(lens[c])])

# ------------------------------------------------------------------ job generation
# job = (family, width, variant, start, count); worker regenerates perms for that slice
SAMPLES = {}
MY_R = {}
def get_R(w):
    if w not in MY_R: MY_R[w] = weak_orders(w)
    return MY_R[w]
def orders_for(fam, w, var):
    key = (fam, w, var)
    if key in SAMPLES: return SAMPLES[key]
    return None
def build(job):
    fam, w, var, s, cnt = job
    if fam == 'MY':
        R = get_R(w); sel = SAMPLES.get(('MY', w)); idx = np.arange(s, s + cnt) if sel is None else sel[s:s + cnt]
        F = mysz_perms(R[idx]); ids = idx
    else:
        sel = SAMPLES.get((fam, w))
        if sel is None: ords = list(itertools.islice(itertools.permutations(range(w)), s, s + cnt)); ids = np.arange(s, s + len(ords))
        else: ords = [tuple(o) for o in sel[s:s + cnt]]; ids = np.arange(s, s + len(ords))
        fn = (lambda o: disrupted_perm(o, var == 'gap')) if fam == 'DV' else (lambda o: irregular_perm(o, var))
        F = np.array([fn(o) for o in ords], dtype=np.int16)
    Inv = np.argsort(F, axis=1).astype(np.int16)
    M = np.concatenate([F, Inv]); labs = [(fam, w, var, int(i), 0) for i in ids] + [(fam, w, var, int(i), 1) for i in ids]
    return M, labs
def describe(lab):
    fam, w, var, i, inv = lab
    if fam == 'MY': k = tuple(get_R(w)[i])          # i is the true weak-order index (also when sampled)
    else:
        sel = SAMPLES.get((fam, w)); k = tuple(sel[i]) if sel is not None else next(itertools.islice(itertools.permutations(range(w)), i, None))
    return f'{fam}/{var} w{w} key={tuple(int(x) for x in k)}' + (' [inverse]' if inv else '')

CT_CUR = None
def worker(job):
    M, labs = build(job); return run3(CT_CUR, M, labs), len(M)

def campaign(ct, name, jobs, validate=True):
    global CT_CUR; CT_CUR = ct; t = time.time()
    if validate:                      # cross-check run3 vs gen.run on the first job
        M, labs = build(jobs[0]); L = np.arange(len(M))
        a = sorted(repr(h) for h in run(ct, ANCH, M, L, Pmax=PKEY, anyalpha=False))
        b = sorted(repr(h) for h in run(ct, ANCH, M, L, Pmax=PANY, RTHR=99) if h[0] == 'any')
        c = sorted(repr(h) for h in run3(ct, M, list(L)))
        assert sorted(a + b) == c, (name, len(a), len(b), len(c))
    with Pool(2) as pool: res = pool.map(worker, jobs, chunksize=1)
    hits = [h for r, _ in res for h in r]; n = sum(k for _, k in res)
    print(f'  {name}: {n} perms scored, {len(hits)} hits, {time.time()-t:.0f}s', flush=True)
    return hits, n

def jobs_for(fam, w, var, total):
    return [(fam, w, var, s, min(CHUNK // 2, total - s)) for s in range(0, total, CHUNK // 2)]

if __name__ == '__main__':
    T0 = time.time()
    MY9_FRAC = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    DV_SAMP = int(sys.argv[2]) if len(sys.argv) > 2 else 100000
    rng = np.random.default_rng(5)
    # ---------------- validation of implementations
    for w in range(3, 8):
        R = get_R(w); sub = R[rng.choice(len(R), min(300, len(R)), replace=False)]
        assert all(np.array_equal(a, mysz_ref(tuple(r))) for a, r in zip(mysz_perms(sub), sub))
    print('weak-order counts', {w: len(get_R(w)) for w in range(3, 10)})
    msg = 'WEAREDISCOVEREDFLEEATONCEXQZJ'
    for gap in (True, False):
        rd, D = disrupted_order((2, 0, 4, 1, 3), gap, len(msg))
        enc = ''.join(msg[i] for i in rd); dec = [''] * len(msg)
        for k, i in enumerate(rd): dec[i] = enc[k]
        assert ''.join(dec) == msg and sorted(rd) == list(range(len(msg)))
        print(f'DV test gap={gap} key-order (2,0,4,1,3) grid:'); print('\n'.join('    ' + ''.join('#' if d else '.' for d in row) for row in D))
        print(f'    enc {enc}  dec-roundtrip OK')
    for w in range(2, 9):
        for o in list(itertools.permutations(range(w)))[:50]:
            for v in ('right', 'keylong'): assert sorted(irregular_perm(o, v).tolist()) == list(range(N))
    print('IR validation: perms are bijections; example right w5 (1,0,2,3,4) ct-order head',
          np.argsort(irregular_perm((1, 0, 2, 3, 4), 'right'))[:22].tolist())
    # ---------------- job lists
    jobs = []; cover = []
    for w in range(3, 9): n = len(get_R(w)); jobs += jobs_for('MY', w, '-', n); cover.append(f'MY w{w} {n}/{n}')
    n9 = len(get_R(9))
    if MY9_FRAC < 1:
        SAMPLES[('MY', 9)] = np.sort(np.random.default_rng(9).choice(n9, int(n9 * MY9_FRAC), replace=False))
    k9 = len(SAMPLES[('MY', 9)]) if ('MY', 9) in SAMPLES else n9
    jobs += jobs_for('MY', 9, '-', k9); cover.append(f'MY w9 {k9}/{n9}')
    for w in range(5, 11):
        tot = math.factorial(w)
        if tot > 40320:
            r2 = np.random.default_rng(w)
            SAMPLES[('DV', w)] = np.array([r2.permutation(w) for _ in range(min(DV_SAMP, tot))])
        k = len(SAMPLES[('DV', w)]) if ('DV', w) in SAMPLES else tot
        for var in ('gap', 'nogap'): jobs += jobs_for('DV', w, var, k)
        cover.append(f'DV w{w} {k}/{tot} x2var')
    for w in range(2, 9):
        for var in ('right', 'keylong'): jobs += jobs_for('IR', w, var, math.factorial(w))
        cover.append(f'IR w{w} all x2var')
    print('coverage:', '; '.join(cover), flush=True)
    # ---------------- positive controls: one planted per family
    def plant(perm, seed):
        r = random.Random(seed); pt = [r.choice(ST) for _ in range(N)]
        for i, c in ANCH.items(): pt[i] = c
        x = [KA[(KA.index(c) + KA.index('OCHRE'[i % 5])) % 26] for i, c in enumerate(pt)]
        z = [''] * N
        for i in range(N): z[perm[i]] = x[i]
        return ''.join(z)
    ctrl = [('MY', 7, '-', int(np.nonzero((get_R(7) == np.array([2, 0, 1, 0, 3, 2, 1])).all(1))[0][0])),
            ('DV', 8, 'gap', list(itertools.permutations(range(8))).index((3, 0, 6, 1, 7, 2, 5, 4))),
            ('IR', 7, 'keylong', list(itertools.permutations(range(7))).index((4, 1, 6, 0, 3, 5, 2)))]
    ok_all = True
    for fam, w, var, i in ctrl:
        M, labs = build((fam, w, var, i, 1)); ct = plant(M[0], 3)
        fj = [j for j in jobs if j[0] == fam and j[1] == w and j[2] == var]
        h, _ = campaign(ct, f'CONTROL {fam} {describe(labs[0])}', fj)
        good = [x for x in h if x[0] == 'key' and x[1:5] == ('KA', 0, 'pt', 5) and x[-1] == labs[0]]
        others = {x[-1] for x in h} - {labs[0]}
        print(f'    recovered={bool(good)}  other perms hitting in that sub-family: {len(others)}', flush=True)
        ok_all &= bool(good)
    # ---------------- K4 and null
    out = {}
    for name, ct in (('K4', CT0), ('NULL', ''.join(random.Random(11).sample(CT0, N)))):
        h, n = campaign(ct, name, jobs); out[name] = (h, n)
        import pickle; pickle.dump(h, open(f'new_tradecraft_{name}.pkl', 'wb'))
        fams = {}
        for x in h: fams.setdefault(x[-1][0], []).append(x)
        print(f'  {name} hits by family:', {k: len(v) for k, v in fams.items()}, '| key-type hits:',
              sum(x[0] == 'key' for x in h))
        for x in sorted(h, key=lambda x: -x[-2] if x[0] == 'key' else -x[3])[:12]: print('     ', x[:-1], describe(x[-1]))
    print(f'SUMMARY perms={out["K4"][1]} K4_hits={len(out["K4"][0])} NULL_hits={len(out["NULL"][0])} '
          f'controls={"OK" if ok_all else "FAIL"} total {time.time()-T0:.0f}s')
