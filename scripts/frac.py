import random,sys
from stack import CT0,ANCH
N=97
def blocks(p,o):
    b=[];s=0
    if o: b.append(list(range(0,o)));s=o
    while s<N: b.append(list(range(s,min(N,s+p))));s+=p
    return b
def constraints(ct,p,o,D):
    """D dims. returns list of (ct_letter, [(ptpos,dim) for each output dim]) where all needed pt positions are anchors"""
    out=[]
    for blk in blocks(p,o):
        L=len(blk);seq=[(blk[j],d) for d in range(D) for j in range(L)]
        for k in range(L):
            src=seq[D*k:D*k+D]
            if all(pp in ANCH for pp,_ in src): out.append((ct[blk[k]],[(ANCH[pp],d) for pp,d in src]))
    return out
class UF:
    def __init__(s):s.p={}
    def f(s,x):
        s.p.setdefault(x,x)
        while s.p[x]!=x: s.p[x]=s.p[s.p[x]];x=s.p[x]
        return x
    def u(s,a,b):s.p[s.f(a)]=s.f(b)
def feasible(cons,D,S):
    uf=UF()
    for c,src in cons:
        for outd,(L,d) in enumerate(src): uf.u((c,outd),(L,d))
    letters=sorted({c for c,_ in cons}|{L for _,s in cons for L,_ in s})
    tup={L:tuple(uf.f((L,d)) for d in range(D)) for L in letters}
    # same tuple for two letters -> impossible
    if len(set(tup.values()))<len(letters): return False
    # colouring: each class (dim d) gets value 0..S-1; letters must get distinct tuples; capacity: per (d,value) at most S**(D-1) letters
    cls={};members={}
    for L in letters:
        for d in range(D):
            r=uf.f((L,d));cls.setdefault(d,set()).add(r);members.setdefault(r,set()).add(L)
    # a class may mix dims after union (e.g. row of X == col of Y) -> fine: one value
    allcls=sorted({uf.f((L,d)) for L in letters for d in range(D)},key=lambda r:-len(members[r]))
    cap=S**(D-1);val={}
    def ok_partial():
        # capacity per dim
        for d in range(D):
            cnt={}
            for L in letters:
                r=uf.f((L,d))
                if r in val: cnt[val[r]]=cnt.get(val[r],0)+1
            if any(v>cap for v in cnt.values()): return False
        seen=set()
        for L in letters:
            t=tuple(val.get(uf.f((L,d))) for d in range(D))
            if None not in t:
                if t in seen: return False
                seen.add(t)
        return True
    nodes=[0]
    def bt(i):
        nodes[0]+=1
        if nodes[0]>200000: return True   # give up -> treat as feasible
        if i==len(allcls): return True
        r=allcls[i]
        for v in range(S):
            val[r]=v
            if ok_partial() and bt(i+1): return True
        del val[r];return False
    return bt(0)
def scan(ct,D,S,label):
    res=[]
    for p in range(2,98):
        for o in (range(p) if p<=24 else [0]):
            cons=constraints(ct,p,o,D)
            if len(cons)>=3 and feasible(cons,D,S): res.append((p,o,len(cons)))
    return res
def bifid_enc(pt,sq,S,D,p):
    pos={c:i for i,c in enumerate(sq)}
    coord=lambda c:[(pos[c]//S**(D-1-d))%S for d in range(D)]
    out=[]
    for blk in blocks(p,0):
        seq=[coord(pt[j])[d] for d in range(D) for j in blk]
        for k in range(len(blk)):
            t=seq[D*k:D*k+D];i=sum(v*S**(D-1-d) for d,v in enumerate(t));out.append(sq[i])
    return ''.join(out)
if __name__=='__main__':
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
    for name,D,S,extra in (('Bifid6x6',2,6,'0123456789'),('Trifid3x3x3',3,3,'+')):
        sq=list('ABCDEFGHIJKLMNOPQRSTUVWXYZ'+extra);random.Random(4).shuffle(sq);sq=''.join(sq)
        c=bifid_enc(pt,sq,S,D,7);ctrl=[r for r in scan(c,D,S,'') if r[:2]==(7,0)]
        nl=[];
        for sd in range(5):
            l=list(CT0);random.Random(sd).shuffle(l);nl.append(len(scan(''.join(l),D,S,'')))
        k4=scan(CT0,D,S,'')
        print(name,'| control period 7 found:',ctrl,'| null feasible counts:',nl,'| K4 feasible:',len(k4),k4[:20],flush=True)
