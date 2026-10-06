import random,math
from stack import solve_mod,CT0,ANCH
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
def blocks(n,o):
    out=[];s=o
    while s+n<=97:
        b=list(range(s,s+n))
        if all(i in ANCH for i in b): out.append(b)
        s+=n
    return out
def test(ct,n,o,A,L=0):
    """c_blk = M p_blk + b[class], class = block index mod L (L=0: pure Hill; L>=1 affine/after periodic key)."""
    bl=blocks(n,o);res=[]
    if not bl: return None
    worst=None
    for row in range(n):
        rows=[];rhs=[]
        for b in bl:
            p=[A.index(ANCH[i]) for i in b];c=A.index(ct[b[row]])
            r=p[:]
            if L: e=[0]*L;e[((b[0]-o)//n)%L]=1;r+=e
            rows.append(r);rhs.append(c)
        ok2,r2=solve_mod(rows,rhs,2);ok13,r13=solve_mod(rows,rhs,13)
        if not(ok2 and ok13): return False,0
        red=len(rows)-max(r2,r13);worst=red if worst is None else min(worst,red)
    return True,worst
def scan(ct):
    hits=[];tested=[]
    for n in (2,3,4,5):
        for o in range(n):
            for an,A in (('ST',ST),('KA',KA)):
                for L in range(0,7):
                    r=test(ct,n,o,A,L)
                    if r is None: continue
                    ok,red=r;tested.append((n,o,an,L,red if ok else None))
                    if ok and red>=2: hits.append((n,o,an,L,red))
    return hits,tested
if __name__=='__main__':
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:98]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
    rng=random.Random(2)
    while True:
        M=[[rng.randrange(26) for _ in range(2)] for _ in range(2)]
        if math.gcd((M[0][0]*M[1][1]-M[0][1]*M[1][0])%26,26)==1: break
    o=1;ct=list(pt[:97])
    for s in range(o,96,2):
        p=[KA.index(pt[s]),KA.index(pt[s+1])]
        for r in range(2): ct[s+r]=KA[(M[r][0]*p[0]+M[r][1]*p[1])%26]
    h,_=scan(''.join(ct));print('control Hill 2x2 KA offset 1:',[x for x in h if x[:3]==(2,1,'KA')][:2])
    h,t=scan(CT0);print('K4 hits (redundancy>=2):',h)
    print('testable configs:',len([x for x in t if True]),' strong ones (redundancy>=2 if consistent):',sorted(set((x[0],x[3]) for x in t if x[4] is None)) [:5],'...')
    import collections;print('per n how many configs consistent at all:',collections.Counter(x[0] for x in t if x[4] is not None))
