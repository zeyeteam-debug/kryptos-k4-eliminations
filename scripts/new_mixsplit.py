# K4 as two segments of DIFFERENT method families (like K1/K2 = Vigenere vs K3 = double rotation).
# Split s in 34..63 so that A=[0,s) holds EASTNORTHEAST (21-33) and B=[s,97) holds BERLINCLOCK (63-73).
# Case 1: A = pure K3-style double rotation (k3rot.compose over length s, w1,w2 in 2..s-1, both directions,
#         perm and its inverse), NO substitution -> ct[perm[i]]==crib[i] exactly for every crib i in A;
#         B = periodic key alone (gen.run: ST/KA Vig/Beau pt/ct-indexed + any-alphabet, Pmax=24).
# Case 2: the reverse (A periodic key, B pure double rotation over length 97-s).
# Usage: python3 new_mixsplit.py K4 | ctrl | nullN
import numpy as np, random, time, sys
from collections import Counter
from stack import CT0, ANCH
from gen import run
from k3rot import rot
KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'; ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
SPLITS=range(34,64)

_cache={}
def rot_perms(L):
    """All distinct perms m (m[pt]=ct) from compose(L,w1,d1,w2,d2), w in 2..L-1, d in T/F, plus inverses."""
    if L in _cache: return _cache[L]
    R={(w,d):rot(L,w,d) for w in range(2,L) for d in (True,False)}
    seen=set(); out=[]
    for a in R.values():
        for b in R.values():
            idx=a[b]                                   # ct[n]=pt[idx[n]]  (decrypt map: m[ct]=pt)
            inv=np.empty(L,dtype=np.int16); inv[idx]=np.arange(L)   # encrypt map: m[pt]=ct
            for m in (inv, idx.astype(np.int16)):
                k=m.tobytes()
                if k not in seen: seen.add(k); out.append(m)
    P=np.array(out,dtype=np.int16); _cache[L]=P; return P

def trans_hits(seg,crib):
    """# of double-rotation perms of len(seg) with seg[m[i]]==crib[i] for all crib positions i."""
    P=rot_perms(len(seg)); pos=np.array(sorted(crib)); let=np.array([ord(crib[i]) for i in pos])
    c=np.array([ord(x) for x in seg])
    return int(np.all(c[P[:,pos]]==let,axis=1).sum())

def feasible(seg,crib):
    """Necessary condition for ANY pure transposition: seg contains the crib's letter multiset."""
    return not (Counter(crib.values())-Counter(seg))

def key_hits(seg,crib,Pmax=24,RTHR=6):
    I=np.arange(len(seg),dtype=np.int16)[None,:]
    return run(seg,crib,I,['id'],Pmax=Pmax,RTHR=RTHR,ANYSTRONG=4)

def analyse(ct,Pmax=24,RTHR=6):
    """returns dict s -> (case1_transA, case1_keyB, case2_keyA, case2_transB, feasA, feasB)"""
    res={}
    for s in SPLITS:
        A=ct[:s]; B=ct[s:]
        cA={i:c for i,c in ANCH.items() if i<s}; cB={i-s:c for i,c in ANCH.items() if i>=s}
        res[s]=(trans_hits(A,cA), len(key_hits(B,cB,Pmax,RTHR)), len(key_hits(A,cA,Pmax,RTHR)), trans_hits(B,cB),
                feasible(A,cA), feasible(B,cB))
    return res

def summarise(tag,res,t0):
    c1=[s for s in SPLITS if res[s][0] and res[s][1]]; c2=[s for s in SPLITS if res[s][2] and res[s][3]]
    print(f"{tag}: case1 (A rot, B key) both-consistent splits {c1} | case2 (A key, B rot) {c2} | "
          f"halves: transA {sum(1 for s in SPLITS if res[s][0])} keyB {sum(1 for s in SPLITS if res[s][1])} "
          f"keyA {sum(1 for s in SPLITS if res[s][2])} transB {sum(1 for s in SPLITS if res[s][3])} | "
          f"multiset-feasible splits: A {sum(1 for s in SPLITS if res[s][4])} B {sum(1 for s in SPLITS if res[s][5])} | {round(time.time()-t0)}s",flush=True)
    return c1,c2

def plant(case):
    pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
    for i,c in ANCH.items(): pt[i]=c
    if case==1:   # s*=48: A double rotation (7,cw)(5,ccw); B = KA Vigenere key OCHRE (period 5 -> red 6)
        s=48; idx=rot(s,7,True)[rot(s,5,False)]; ctA=''.join(pt[i] for i in idx)
        ctB=''.join(KA[(KA.index(pt[i])+KA.index('OCHRE'[(i-s)%5]))%26] for i in range(s,97))
    else:         # s*=52: A = ST Beaufort key KRYPTOS (period 7 -> red 6); B double rotation (6,ccw)(9,cw)
        s=52; ctA=''.join(ST[(ST.index('KRYPTOS'[i%7])-ST.index(pt[i]))%26] for i in range(s))
        L=97-s; idx=rot(L,6,False)[rot(L,9,True)]; ctB=''.join(pt[s+i] for i in idx)
    return s, ctA+ctB

if __name__=='__main__':
    t0=time.time(); which=sys.argv[1]
    if which=='K4':
        res=analyse(CT0); summarise('K4',res,t0)
        # diagnostics: where does the letter-multiset condition fail?
        needA=Counter('EASTNORTHEAST'); needB=Counter('BERLINCLOCK')
        for s in (34,48,63):
            print(f"  s={s}: A[0:{s}) lacks {dict(needA-Counter(CT0[:s]))}; B[{s}:97) lacks {dict(needB-Counter(CT0[s:]))}")
        print('  E at',[i for i,c in enumerate(CT0) if c=='E'],' L at',[i for i,c in enumerate(CT0) if c=='L'])
        # relaxed key threshold (red>=3) just to show the key half is not the bottleneck
        r3=analyse(CT0,RTHR=3); print('  relaxed RTHR=3: keyB splits',sum(1 for s in SPLITS if r3[s][1]),'keyA splits',sum(1 for s in SPLITS if r3[s][2]))
    elif which=='ctrl':
        for case in (1,2):
            s,ct=plant(case); res=analyse(ct); c1,c2=summarise(f'ctrl case{case} planted s*={s}',res,t0)
            print(f"  at s*={s}: transA {res[s][0]} keyB {res[s][1]} keyA {res[s][2]} transB {res[s][3]}")
    elif which.startswith('null'):
        n=int(which[4:]); tot=np.zeros(4); both=[0,0]
        for k in range(n):
            l=list(CT0); random.Random(1000+k).shuffle(l); res=analyse(''.join(l))
            tot+=[sum(1 for s in SPLITS if res[s][j]) for j in range(4)]
            both[0]+=sum(1 for s in SPLITS if res[s][0] and res[s][1]); both[1]+=sum(1 for s in SPLITS if res[s][2] and res[s][3])
        print(f"null x{n}: mean splits/shuffle transA {tot[0]/n:.3f} keyB {tot[1]/n:.3f} keyA {tot[2]/n:.3f} transB {tot[3]/n:.3f} | "
              f"both case1 {both[0]/n:.3f} case2 {both[1]/n:.3f} | {round(time.time()-t0)}s",flush=True)
