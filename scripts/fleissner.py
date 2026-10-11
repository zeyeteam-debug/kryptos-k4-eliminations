# Turning grille (Fleissner) 10x10, 25 holes, 4 rotations; 97 letters + 3 blank cells (last 3 plaintext slots).
# Model: grille transposition + periodic key (ST/KA, Vigenere/Beaufort, periods 1..12, key indexed by pt or ct position).
# Search: simulated annealing over the 4^25 grilles maximising crib agreement; positive control + shuffled null.
import numpy as np,random,math,sys,time
from stack import CT0,ANCH
KA='KRYPTOSABCDEFGHIJLMNQUVWXZ';ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';n=10
def rot(r,c,d): return (c,n-1-r) if d==1 else (n-1-c,r)
ORB=[]
seen=set()
for r in range(n):
    for c in range(n):
        if (r,c) in seen: continue
        o=[(r,c)]
        for _ in range(3): o.append(rot(*o[-1],1))
        seen.update(o);ORB.append(o)
AN=np.array(sorted(ANCH))
def perm_of(q,d=1,read='row'):
    holes=[ORB[k][q[k]] for k in range(25)];order=[]   # plaintext index -> cell
    H=holes
    for t in range(4):
        cells=sorted(H) if read=='row' else sorted(H,key=lambda x:(x[1],x[0]))
        order+=cells;H=[rot(r,c,d) for r,c in H]
    filled=order[:97];key=(lambda x:x) if read=='row' else (lambda x:(x[1],x[0]))
    rank={cell:i for i,cell in enumerate(sorted(filled,key=key))}
    return np.array([rank[cell] for cell in filled])
PER=list(range(1,13))
def setup(ct):
    out=[]
    for A in (ST,KA):
        lut=np.zeros(128,int);lut[[ord(x) for x in A]]=np.arange(26)
        out.append((lut[np.array([ord(x) for x in ct])],np.array([A.index(ANCH[i]) for i in AN])))
    return out
def score(perm,S):
    best=0
    pos=perm[AN]
    for C,P in S:
        c=C[pos]
        for K in ((c-P)%26,(c+P)%26):
            for idx in (AN,pos):
                for p in PER:
                    cls=idx%p;sc=0
                    # majority per class
                    code=cls*26+K;cnt=np.bincount(code,minlength=p*26).reshape(p,26).max(axis=1).sum()
                    sc=cnt-(24-24)  # matches
                    # penalise trivial: require redundancy
                    if 24-len(np.unique(cls))<6: continue
                    if cnt>best: best=cnt
    return best
def anneal(ct,iters=20000,seed=0,d=1,read='row'):
    rng=random.Random(seed);S=setup(ct);q=[rng.randrange(4) for _ in range(25)]
    cur=score(perm_of(q,d,read),S);best=(cur,q[:]);T=2.0
    for it in range(iters):
        k=rng.randrange(25);old=q[k];q[k]=rng.choice([x for x in range(4) if x!=old])
        s=score(perm_of(q,d,read),S)
        if s>=cur or rng.random()<math.exp((s-cur)/T): cur=s
        else: q[k]=old
        if cur>best[0]: best=(cur,q[:])
        T=max(0.05,2.0*(1-it/iters))
        if best[0]==24: break
    return best
if __name__=='__main__':
    which=sys.argv[1];runs=int(sys.argv[2])
    t=time.time()
    if which=='ctrl':
        rng=random.Random(99);qt=[rng.randrange(4) for _ in range(25)];pm=perm_of(qt)
        pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
        for i,c in ANCH.items(): pt[i]=c
        ct=['?']*97
        for i,c in enumerate(pt): ct[pm[i]]=KA[(KA.index(c)+KA.index('ABSCI'[i%5]))%26]
        ct=''.join(ct)
    elif which=='K4': ct=CT0
    else:
        l=list(CT0);random.Random(int(which[4:])).shuffle(l);ct=''.join(l)
    res=[]
    for s in range(runs):
        for d in (1,-1):
            b=anneal(ct,seed=s,d=d);res.append((b[0],s,d))
    res.sort(reverse=True);print(which,'best scores (of 24):',[r[0] for r in res[:8]],'time',round(time.time()-t),flush=True)
