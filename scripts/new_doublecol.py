"""Double columnar transposition (two successive keyed columnar transpositions,
different widths) composed with a periodic key underneath (ST/KA Vig/Beau + any-alphabet).
Composite m2[m1[pt]] for widths w1,w2 in WIDTHS (forward/encipher direction only),
plus the full inverse of each composite; dedupe by bytes; periodic-key test Pmax=12.
run2() is a faster re-implementation of gen.run (same hit tuples); it is validated against
gen.run on the first VALID chunk of every campaign.
Positive control: random plaintext + cribs, KA Vigenere 'OCHRE', columnar 4 then columnar 5.
Null: shuffled K4.
usage: python3 new_doublecol.py [widths e.g. 2,3,4,5,6] [nproc]
"""
import sys, time, itertools, random
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, '/tmp/claude-0/-home-claude/bd15c051-bac5-5123-8757-e9921faf2557/scratchpad')
from stack import CT0, ANCH, KA
from gen import run, colperms, ST

N = 97
WIDTHS = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [2, 3, 4, 5, 6]
NPROC = int(sys.argv[2]) if len(sys.argv) > 2 else 2
PMAX = 12
CHUNK = 20000
VALID = 2          # chunks per campaign cross-checked against gen.run

# ---------------------------------------------------------------- fast scorer
AN = np.array(sorted(ANCH))
PR = np.array(list(itertools.combinations(range(len(AN)), 2)))
P0, P1 = PR[:, 0], PR[:, 1]
PT_PAIRS = {Per: np.nonzero((AN % Per)[P0] == (AN % Per)[P1])[0] for Per in range(1, PMAX + 1)}
PT_RED = {Per: len(AN) - len(set((AN % Per).tolist())) for Per in range(1, PMAX + 1)}
Pl = np.array([ord(ANCH[i]) for i in AN]); Peq = Pl[P0] == Pl[P1]

def run2(ct, Ms, labs, Pmax=PMAX, RTHR=6, ANYSTRONG=4):
    """Same semantics/hit tuples as gen.run for ct without '?' (all anchors valid)."""
    assert '?' not in ct
    hits = []
    CI = Ms[:, AN].astype(np.int16)
    cl = np.frombuffer(ct.encode(), dtype=np.uint8); Cl = cl[CI]
    Rct = {Per: (CI % Per).astype(np.int8) for Per in range(1, Pmax + 1)}
    NEEDct = {Per: Rct[Per][:, P0] == Rct[Per][:, P1] for Per in range(1, Pmax + 1)}
    for an, A in (('ST', ST), ('KA', KA)):
        lut = np.zeros(128, dtype=np.int8); lut[[ord(c) for c in A]] = np.arange(26)
        C = lut[Cl]; P = np.array([A.index(ANCH[i]) for i in AN], dtype=np.int8)
        for rule in (0, 1):
            K = ((C - P) % 26 if rule == 0 else (C + P) % 26).astype(np.int8)
            same = K[:, P0] == K[:, P1]
            notsame = ~same
            for Per in range(1, Pmax + 1):
                # kind 'pt': residues independent of perm
                ok = same[:, PT_PAIRS[Per]].all(axis=1)
                red = PT_RED[Per]
                if red >= RTHR:
                    for x in np.nonzero(ok)[0]: hits.append(('key', an, rule, 'pt', Per, red, labs[x]))
                # kind 'ct'
                ok = ~(NEEDct[Per] & notsame).any(axis=1)
                for x in np.nonzero(ok)[0]:
                    rv = Rct[Per][x]; red = len(rv) - len(set(rv.tolist()))
                    if red >= RTHR: hits.append(('key', an, rule, 'ct', Per, red, labs[x]))
    # any-alphabet test
    Ceq = Cl[:, P0] == Cl[:, P1]; bad = Peq[None, :] ^ Ceq
    for Per in range(1, Pmax + 1):
        pp = PT_PAIRS[Per]; strong = int(Peq[pp].sum())
        if strong >= ANYSTRONG:
            ok = ~bad[:, pp].any(axis=1)
            for x in np.nonzero(ok)[0]: hits.append(('any', 'pt', Per, strong, labs[x]))
        sm = NEEDct[Per]
        ok = ~(bad & sm).any(axis=1) & ((sm & Peq[None, :]).sum(axis=1) >= ANYSTRONG)
        for x in np.nonzero(ok)[0]: hits.append(('any', 'ct', Per, int((sm[x] & Peq).sum()), labs[x]))
    return hits

# ---------------------------------------------------------------- build perms
t0 = time.time()
FWD = {}
for w in WIDTHS:
    allp = colperms(N, w)
    FWD[w] = allp[:len(allp) // 2]          # forward (encipher) direction only
orders = {w: list(itertools.permutations(range(w))) for w in WIDTHS}
pairs = list(itertools.product(WIDTHS, WIDTHS))
raw = 2 * sum(len(FWD[w1]) * len(FWD[w2]) for w1, w2 in pairs)
M = np.empty((raw, N), dtype=np.int16); L = np.empty(raw, dtype=np.int64)   # preallocated (int16, ~2 bytes/entry)
block_meta = []; pos = 0; IOTA = np.arange(N, dtype=np.int16)
for bid, (w1, w2) in enumerate(pairs):
    F1 = FWD[w1]; F2 = FWD[w2]; n1 = len(F1); n2 = len(F2); n = n1 * n2
    comp = M[pos:pos + n]; comp[...] = F2[:, F1].reshape(n, N)          # comp[k=b*n1+a, i] = F2[b, F1[a, i]]
    inv = M[pos + n:pos + 2 * n]; inv[np.arange(n)[:, None], comp] = IOTA[None, :]   # full inverse of the composite
    k = np.arange(n, dtype=np.int64)
    L[pos:pos + n] = (bid << 40) + k * 2; L[pos + n:pos + 2 * n] = (bid << 40) + k * 2 + 1
    block_meta.append((w1, w2, n1)); pos += 2 * n
assert pos == raw
Mv = M.view(np.dtype((np.void, M.dtype.itemsize * N))).ravel()
first = np.unique(Mv, return_index=True)[1]; first.sort(); del Mv
M = np.ascontiguousarray(M[first]); L = L[first]; del first
print(f'widths {WIDTHS}: raw composites {raw}, unique perms {len(M)}, '
      f'mem {M.nbytes/1e6:.0f} MB, build {time.time()-t0:.1f}s', flush=True)

def decode(lab):
    lab = int(lab); bid = lab >> 40; rest = lab & ((1 << 40) - 1); k = rest >> 1; inv = rest & 1
    w1, w2, n1 = block_meta[bid]; b, a = divmod(k, n1)
    return f'T1 w{w1} {orders[w1][a]} -> T2 w{w2} {orders[w2][b]}' + (' [inverse]' if inv else '')

CT_CUR = None
def worker(s):
    return run2(CT_CUR, M[s:s + CHUNK], L[s:s + CHUNK])

def campaign(ct, name):
    global CT_CUR
    CT_CUR = ct; t = time.time()
    starts = list(range(0, len(M), CHUNK))
    # cross-check fast scorer vs gen.run on the first VALID chunks
    for s in starts[:VALID]:
        a = sorted(map(repr, run(ct, ANCH, M[s:s + CHUNK], L[s:s + CHUNK], Pmax=PMAX)))
        b = sorted(map(repr, run2(ct, M[s:s + CHUNK], L[s:s + CHUNK])))
        assert a == b, (name, s, a[:3], b[:3])
    with Pool(NPROC) as pool:
        res = pool.map(worker, starts, chunksize=1)
    hits = [h for r in res for h in r]
    print(f'{name}: {len(hits)} hits in {time.time()-t:.1f}s (fast scorer == gen.run on first {VALID} chunks)', flush=True)
    for h in hits[:20]: print('   ', h[:-1], decode(h[-1]))
    return hits

if __name__ == '__main__':
    # ---- positive control ---------------------------------------------------
    rng = random.Random(7)
    pt = [rng.choice(ST) for _ in range(N)]
    for i, c in ANCH.items(): pt[i] = c
    key = 'OCHRE'
    x = [KA[(KA.index(c) + KA.index(key[i % 5])) % 26] for i, c in enumerate(pt)]   # KA Vigenere, pt-index period 5
    o1 = (2, 0, 3, 1); o2 = (3, 1, 4, 0, 2)
    m1 = FWD[4][orders[4].index(o1)]; m2 = FWD[5][orders[5].index(o2)]; comp = m2[m1]
    z = [''] * N
    for i in range(N): z[comp[i]] = x[i]
    ctrl_ct = ''.join(z)
    ctrl_hits = campaign(ctrl_ct, 'CONTROL')
    want = [h for h in ctrl_hits if h[0] == 'key' and h[1] == 'KA' and h[2] == 0 and h[3] == 'pt' and h[4] == 5
            and np.array_equal(M[np.nonzero(L == h[-1])[0][0]], comp)]
    print('control recovered (KA vig pt-period 5, correct perm):', bool(want), flush=True)
    # ---- K4 ------------------------------------------------------------------
    k4_hits = campaign(CT0, 'K4')
    # ---- null: shuffled K4 ---------------------------------------------------
    sh = list(CT0); random.Random(11).shuffle(sh)
    null_hits = campaign(''.join(sh), 'NULL')
    import resource
    print(f'SUMMARY widths={WIDTHS} perms={len(M)} K4_hits={len(k4_hits)} null_hits={len(null_hits)} '
          f'control={"OK" if want else "FAIL"} total {time.time()-t0:.0f}s '
          f'peakRSS(main) {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024:.0f} MB '
          f'peakRSS(worker) {resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss/1024:.0f} MB')
