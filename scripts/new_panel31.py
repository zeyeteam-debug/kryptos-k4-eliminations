"""K4 'ID BY ROWS' test: K4 is carved as a 4x31 block (OBKR top-right of row 0 + 3 full rows).
Hypothesis: the cipher stream was laid into / read out of that block in a non-row order
(columns, serpentine columns, diagonals, spirals, ...), before or after a periodic substitution.
For every route we test M=route and M=route^-1 (and the reversed route), with gen.run
(ST/KA Vig+Beaufort periodic keys, any-alphabet isomorphism check), Pmax=24.
Shapes: width 31 (4 alignments of the 4-letter partial row), widths 30/32/33 likewise,
and the 4x31 block padded with the 27 preceding K3 cells marked '?' (anchor offset 0..27).
Controls: planted KA-Vigenere 'SHADOW' + column-wise fill (97-cell and padded), shuffled-K4 nulls.
"""
import sys, os, random, math, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from stack import CT0, ANCH, KA
import gen

t0 = time.time()
N = 97

# ---------------------------------------------------------------- routes on an R x C rectangle
def spiral(R, C, corner, cw):
    dirs = [(0, 1), (1, 0), (0, -1), (-1, 0)]            # E S W N (clockwise order)
    start = {'TL': (0, 0), 'TR': (0, C - 1), 'BR': (R - 1, C - 1), 'BL': (R - 1, 0)}[corner]
    d = ({'TL': 0, 'TR': 1, 'BR': 2, 'BL': 3} if cw else {'TL': 1, 'TR': 2, 'BR': 3, 'BL': 0})[corner]
    seen = set(); seq = []; r, c = start
    while len(seq) < R * C:
        seq.append((r, c)); seen.add((r, c))
        for _ in range(4):
            nr, nc = r + dirs[d][0], c + dirs[d][1]
            if 0 <= nr < R and 0 <= nc < C and (nr, nc) not in seen:
                break
            d = (d + (1 if cw else -1)) % 4
        else:
            break
        r, c = nr, nc
    assert len(seq) == R * C
    return seq

def rect_routes(R, C):
    rows = range(R); cols = range(C)
    cells = [(r, c) for r in rows for c in cols]
    rt = {}
    rt['rows_LR'] = cells                                            # identity (baseline)
    rt['rows_RL'] = [(r, c) for r in rows for c in reversed(cols)]
    rt['rows_boustro_L'] = [(r, c) for r in rows for c in (cols if r % 2 == 0 else reversed(cols))]
    rt['rows_boustro_R'] = [(r, c) for r in rows for c in (reversed(cols) if r % 2 == 0 else cols)]
    rt['rows_BU_LR'] = [(r, c) for r in reversed(rows) for c in cols]
    rt['rows_BU_RL'] = [(r, c) for r in reversed(rows) for c in reversed(cols)]
    rt['cols_TD_LR'] = [(r, c) for c in cols for r in rows]
    rt['cols_BU_LR'] = [(r, c) for c in cols for r in reversed(rows)]
    rt['cols_TD_RL'] = [(r, c) for c in reversed(cols) for r in rows]
    rt['cols_BU_RL'] = [(r, c) for c in reversed(cols) for r in reversed(rows)]
    for cdir, cl in (('LR', list(cols)), ('RL', list(reversed(cols)))):
        for start in ('TD', 'BU'):
            seq = []
            for k, c in enumerate(cl):
                down = (k % 2 == 0) == (start == 'TD')
                seq += [(r, c) for r in (rows if down else reversed(rows))]
            rt[f'serpcols_{start}_{cdir}'] = seq
    for typ in ('anti', 'main'):
        key = (lambda rc: rc[0] + rc[1]) if typ == 'anti' else (lambda rc: rc[0] - rc[1])
        ks = sorted(set(key(rc) for rc in cells))
        for kord in ('asc', 'desc'):
            kl = ks if kord == 'asc' else list(reversed(ks))
            for rord in ('rasc', 'rdesc'):
                seq = []
                for k in kl:
                    d = [rc for rc in cells if key(rc) == k]
                    d.sort(key=lambda rc: rc[0], reverse=(rord == 'rdesc'))
                    seq += d
                rt[f'diag_{typ}_{kord}_{rord}'] = seq
    for corner in ('TL', 'TR', 'BR', 'BL'):
        for cw in (True, False):
            rt[f'spiral_{corner}_{"cw" if cw else "ccw"}'] = spiral(R, C, corner, cw)
    return rt

# ---------------------------------------------------------------- shapes (cells in carved reading order)
def shapes_for_width(W, n=N):
    R = math.ceil(n / W); part = n - (R - 1) * W
    full = [(r, c) for r in range(R) for c in range(W)]
    out = {}
    out[f'W{W}_first_right'] = [(0, c) for c in range(W - part, W)] + [rc for rc in full if rc[0] > 0]   # sculpture layout for W=31
    out[f'W{W}_first_left'] = [(0, c) for c in range(part)] + [rc for rc in full if rc[0] > 0]
    out[f'W{W}_last_left'] = [rc for rc in full if rc[0] < R - 1] + [(R - 1, c) for c in range(part)]
    out[f'W{W}_last_right'] = [rc for rc in full if rc[0] < R - 1] + [(R - 1, c) for c in range(W - part, W)]
    for v in out.values(): assert len(v) == n
    return R, out

def perms_for_shape(R, C, cells):
    """returns (Ms int16 array, labels). cells = K4 cells in reading order (defines ct index)."""
    idx = {rc: i for i, rc in enumerate(cells)}
    n = len(cells); out = []; labs = []
    for name, route in rect_routes(R, C).items():
        f = [idx[rc] for rc in route if rc in idx]
        assert len(f) == n
        for rev, seq in (('', f), ('_rev', f[::-1])):
            fwd = np.array(seq, dtype=np.int16)
            inv = np.empty(n, dtype=np.int16); inv[fwd] = np.arange(n, dtype=np.int16)
            out.append(fwd); labs.append(f'{name}{rev}/fwd')
            out.append(inv); labs.append(f'{name}{rev}/inv')
    return np.array(out), labs

def dedupe(Ms, labs):
    seen = {}; keep = []
    for i, m in enumerate(Ms):
        k = m.tobytes()
        if k not in seen: seen[k] = i; keep.append(i)
    return Ms[keep], [labs[i] for i in keep]

# ---------------------------------------------------------------- build the 97-cell permutation set
all_M = []; all_L = []
for W in (31, 30, 32, 33):
    R, shp = shapes_for_width(W)
    for sname, cells in shp.items():
        Ms, labs = perms_for_shape(R, W, cells)
        all_M.append(Ms); all_L += [f'{sname}/{l}' for l in labs]
M97 = np.concatenate(all_M); L97 = all_L
n_raw97 = len(L97)
M97, L97 = dedupe(M97, L97)

# padded 4x31 = 124 cells; cells 0..26 (K3 tail + '?') are unknown
RP, CP = 4, 31
cellsP = [(r, c) for r in range(RP) for c in range(CP)]
MP, LP = perms_for_shape(RP, CP, cellsP)
n_rawP = len(LP)
MP, LP = dedupe(MP, LP)
LP = ['PAD4x31/' + l for l in LP]
CTP = '?' * 27 + CT0
assert len(CTP) == 124

print(f'97-cell perms: {len(L97)} unique (of {n_raw97} built; widths 31/30/32/33 x 4 alignments x 30 routes x rev x inv)')
print(f'padded 4x31 perms: {len(LP)} unique (of {n_rawP}); anchor offsets 0..27 -> {len(LP)*28} tests')

PMAX = 24
def run97(ct, tag):
    h = gen.run(ct, ANCH, M97, L97, Pmax=PMAX)
    return h
def runpad(ct, tag):
    hits = []
    for s in range(28):
        anch = {s + i: v for i, v in ANCH.items()}
        for h in gen.run(ct, anch, MP, LP, Pmax=PMAX):
            hits.append(h + (f'off{s}',))
    return hits

def summarize(hits):
    key = [h for h in hits if h[0] == 'key']; anyh = [h for h in hits if h[0] == 'any']
    return key, anyh

# ---------------------------------------------------------------- K4
hk, ha = summarize(run97(CT0, 'K4') + runpad(CTP, 'K4pad'))
print(f'\nK4: key-model hits (red>=6) = {len(hk)}; any-alphabet hits (>=4 strong) = {len(ha)}')
for h in sorted(hk, key=lambda x: -x[5])[:15]: print('   KEY', h)
from collections import Counter
print('   any-alphabet hits by period:', sorted(Counter(h[2] for h in ha).items()))
for h in sorted(ha, key=lambda x: -x[3])[:10]: print('   ANY', h)

# ---------------------------------------------------------------- nulls: shuffled K4 letters
rng = random.Random(97)
nk = []; na = []
for trial in range(5):
    L = list(CT0); rng.shuffle(L); cs = ''.join(L)
    k, a = summarize(run97(cs, 'null') + runpad('?' * 27 + cs, 'nullpad'))
    nk.append(len(k)); na.append(len(a))
print(f'\nNULL (5 shuffles of K4): key-model hits per trial = {nk}; any-alphabet hits per trial = {na}')

# ---------------------------------------------------------------- positive controls
def kavig(pt, key):
    return ''.join(KA[(KA.index(c) + KA.index(key[i % len(key)])) % 26] for i, c in enumerate(pt))
ptr = [rng.choice('ABCDEFGHIJKLMNOPQRSTUVWXYZ') for _ in range(N)]
for i, v in ANCH.items(): ptr[i] = v
pt97 = ''.join(ptr)
c97 = kavig(pt97, 'SHADOW')
# control A: stream written column-wise (top-down, L->R) into the sculpture-shaped 97 cells, carved row-wise
R31, shp31 = shapes_for_width(31)
cells = shp31['W31_first_right']; idx = {rc: i for i, rc in enumerate(cells)}
route = [idx[rc] for rc in rect_routes(4, 31)['cols_TD_LR'] if rc in idx]
carved = [''] * N
for k in range(N): carved[route[k]] = c97[k]
ctrlA = ''.join(carved)
hkA, haA = summarize(run97(ctrlA, 'ctrlA'))
want = [h for h in hkA if h[6] == 'W31_first_right/cols_TD_LR/fwd' and h[1] == 'KA' and h[2] == 0 and h[3] == 'pt' and h[4] == 6]
print(f'\nCONTROL A (KA vig SHADOW, column-wise fill of sculpture block): key hits={len(hkA)}, expected signature found={bool(want)}')
for h in sorted(hkA, key=lambda x: -x[5])[:4]: print('   ', h)
# control B: padded 124-stream, K4 plaintext at stream offset 27, column-wise fill of full 4x31, K3 cells masked
stream = [rng.choice('ABCDEFGHIJKLMNOPQRSTUVWXYZ') for _ in range(124)]
stream[27:] = list(pt97); c124 = kavig(''.join(stream), 'SHADOW')
routeP = [r * 31 + c for (r, c) in rect_routes(4, 31)['cols_TD_LR']]
carvedP = [''] * 124
for k in range(124): carvedP[routeP[k]] = c124[k]
ctrlB = '?' * 27 + ''.join(carvedP[27:])
hkB, haB = summarize(runpad(ctrlB, 'ctrlB'))
wantB = [h for h in hkB if h[6] == 'PAD4x31/cols_TD_LR/fwd' and h[1] == 'KA' and h[3] == 'pt' and h[4] == 6 and h[7] == 'off27']
print(f'CONTROL B (padded 4x31, offset 27, column-wise): key hits={len(hkB)}, expected signature found={bool(wantB)}')
for h in sorted(hkB, key=lambda x: -x[5])[:4]: print('   ', h)
print(f'\nelapsed {time.time()-t0:.1f}s')
