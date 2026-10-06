import numpy as np,itertools,sys,random,math
from stack import CT0,ANCH
from tauto import perms,colmap
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
N=97
def cands(n): return np.array(list(itertools.product(range(26),repeat=n)),dtype=np.int64).T   # n x 26^n
CAND={2:cands(2),3:cands(3)}
def orderB(ct,M,Minv,n,A,dname):
    """Hill on PT order, then transposition. X[i]=ct[minv[i]]"""
    hits=[];ci=np.array([A.index(x) for x in ct])
    for o in range(n):
        bl=[list(range(s,s+n)) for s in range(o,N-n+1,n) if all(i in ANCH for i in range(s,s+n))]
        if len(bl)<=n: continue
        E=np.array([[A.index(ANCH[i]) for i in b] for b in bl])         # k x n
        pred=(E@CAND[n])%26                                             # k x 26^n
        okall=np.ones(len(Minv),bool)
        for r in range(n):
            S=set(map(tuple,pred.T.tolist()))
            pos=np.array([b[r] for b in bl])
            rhs=ci[Minv[:,pos]]                                         # nperm x k
            ok=np.array([tuple(x) in S for x in rhs.tolist()])
            okall&=ok
            if not okall.any(): break
        for x in np.nonzero(okall)[0]: hits.append(('B',n,o,dname,x,len(bl)))
    return hits
def orderA(ct,M,n,A,dname,maxk=None):
    """transposition first (Y[j]=pt[M[j]]), then Hill on Y order."""
    hits=[];ci=np.array([A.index(x) for x in ct]);C=CAND[n]
    known=np.full(N,-1);
    for i,p in ANCH.items(): known[i]=A.index(p)
    Y=known[M]                                                            # nperm x 97 (letter idx or -1)
    for o in range(n):
        starts=np.arange(o,N-n+1,n)
        blk=np.stack([Y[:,starts+r] for r in range(n)],axis=2)           # nperm x nb x n
        full=(blk>=0).all(axis=2)                                        # nperm x nb
        cnt=full.sum(axis=1);cand=np.nonzero(cnt>n)[0]
        for x in cand:
            sel=starts[full[x]];E=blk[x][full[x]]                        # k x n
            pred=(E@C)%26;good=True
            for r in range(n):
                rhs=ci[sel+r]
                if not (pred==rhs[:,None]).all(axis=0).any(): good=False;break
            if good: hits.append(('A',n,o,dname,x,len(sel)))
    return hits
def run(ct,Ws,ns=(2,3),label=''):
    H=[]
    for W in Ws:
        M,Mi=perms(W);M=M.astype(np.int64);Mi=Mi.astype(np.int64)
        for dname,m,mi in (('enc',M,Mi),('dec',Mi,M)):
            for an,A in (('ST',ST),('KA',KA)):
                for n in ns:
                    H+=[h+(an,W) for h in orderB(ct,m,mi,n,A,dname)]
                    if n==2 or W<=8: H+=[h+(an,W) for h in orderA(ct,m,n,A,dname)]
        print(label,'W',W,'hits',len(H),flush=True)
    return H
def hill(text,Mx,A,o):
    n=len(Mx);out=list(text)
    for s in range(o,len(text)-n+1,n):
        p=[A.index(text[s+j]) for j in range(n)]
        for r in range(n): out[s+r]=A[sum(Mx[r][j]*p[j] for j in range(n))%26]
    return ''.join(out)
if __name__=='__main__':
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
    Mx=[[3,5],[7,2]];idx=colmap(5,(2,0,4,1,3))
    x=hill(pt,Mx,KA,1);ctB=''.join(x[i] for i in idx)                  # Hill then transpose
    y=''.join(pt[i] for i in idx);ctA=hill(y,Mx,KA,0)                    # transpose then Hill
    hb=run(ctB,[5],(2,),'ctrlB');print('control B found:',[h for h in hb if h[0]=='B'][:2])
    ha=run(ctA,[5],(2,),'ctrlA');print('control A found:',[h for h in ha if h[0]=='A'][:2])
    Ws=[int(w) for w in sys.argv[1].split(',')]
    l=list(CT0);random.Random(1).shuffle(l);print('null hits',len(run(''.join(l),[w for w in Ws if w<=7],(2,3),'null')))
    h=run(CT0,Ws,(2,3),'K4');print('K4 hits',len(h),h[:10])
