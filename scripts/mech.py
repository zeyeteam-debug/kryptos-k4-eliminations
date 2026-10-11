# Gear / safe-dial / clock-hand mechanisms as K4 key or transposition. Exact linear tests mod 2 and mod 13.
import numpy as np, itertools, sys, time, random
from stack import CT0, ANCH
from gen import colperms
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
AN=np.array(sorted(ANCH));N=97
def nullspace_left(A,m):
    # basis of {y : y A = 0 mod m}, A is r x u ; returns k x r
    M=(A.T%m).astype(np.int64).copy();u,r=M.shape;piv=[];row=0
    for col in range(r):
        p=None
        for i in range(row,u):
            if M[i,col]%m: p=i;break
        if p is None: continue
        M[[row,p]]=M[[p,row]];M[row]=(M[row]*pow(int(M[row,col]),m-2,m))%m
        for i in range(u):
            if i!=row and M[i,col]: M[i]=(M[i]-M[i,col]*M[row])%m
        piv.append(col);row+=1
        if row==u: break
    free=[c for c in range(r) if c not in piv];B=[]
    for f in free:
        y=np.zeros(r,dtype=np.int64);y[f]=1
        for i,pc in enumerate(piv): y[pc]=(-M[i,f])%m
        B.append(y)
    return np.array(B,dtype=np.int64).reshape(len(B),r)
# ---------- models: coefficient matrices over the 24 crib (pt) positions ----------
def gear(periods):
    cols=[]
    for p in periods:
        for r in range(p): cols.append((AN%p==r).astype(int))
    cols.append(np.zeros(len(AN),int))  # placeholder removed below
    return np.array(cols[:-1]).T
def clock(r,s): return np.array([AN, (AN+s)//r, np.ones(len(AN),int)]).T
def clock3(r1,r2,s): return np.array([AN,(AN+s)//r1,(AN+s)//(r1*r2),np.ones(len(AN),int)]).T
YARnums={'YAR_ST0':(24,0,17),'YAR_ST1':(25,1,18),'YAR_KA0':(2,7,1),'YAR_KA1':(3,8,2)}
ERRnums={'QUAE_ST0':(16,20,0,4),'QUAE_ST1':(17,21,1,5),'LOEI_ST0':(11,14,4,8),'LOEI_ST1':(12,15,5,9),
 'QUAE_KA':tuple(KA.index(c) for c in 'QUAE'),'LOEI_KA':tuple(KA.index(c) for c in 'LOEI'),
 'shifts_ST':tuple((ST.index(b)-ST.index(a))%26 for a,b in ('LQ','OU','EA','IE')),
 'shifts_KA':tuple((KA.index(b)-KA.index(a))%26 for a,b in ('LQ','OU','EA','IE'))}
special={**YARnums,**ERRnums}
models=[]
for p1 in range(2,31):
    for p2 in range(p1+1,31): models.append((('gear',p1,p2),gear((p1,p2))))
for nm,t in special.items():
    ps=tuple(sorted(set(x if x>0 else 26 for x in t)-{1}))
    if len(ps)>=1: models.append((('gear',nm)+ps,gear(ps)))
for r in range(2,61):
    for s in range(r): models.append((('clock',r,s),clock(r,s)))
for r1,r2 in ((12,12),(24,12),(12,5),(60,12),(5,12)):
    for s in range(r1): models.append((('clock3',r1,r2,s),clock3(r1,r2,s)))
for nm,t in special.items():  # numbers as tooth ratios for clock hands
    for x in set(t):
        if x>=2 and not any(m[0]==('clock',x,0) for m in models[:0]): pass
ALL=list(models)
def build(sel):
    global models,NM,GM,RED,RED2,info
    models=sel;Ns={2:[],13:[]};grp={2:[],13:[]};info=[]
    for k,(lab,A) in enumerate(models):
        n2=nullspace_left(A,2);n13=nullspace_left(A,13);info.append((lab,len(n13),len(n2)))
        for m,nn in ((2,n2),(13,n13)):
            for y in nn: Ns[m].append(y);grp[m].append(k)
    NM={m:np.array(Ns[m]).T.astype(np.float32) for m in (2,13)};GM={m:np.array(grp[m]) for m in (2,13)}
    RED=np.array([i[1] for i in info]);RED2=np.array([i[2] for i in info])
build(ALL)
def check(ct,Ms,labs,minred=1,variants=None):
    """Ms[n][pt]=ct idx. returns hits and expected false positives."""
    cl=np.array([ord(c) for c in ct]);hits=[];exp=0.0
    w=np.where(RED>=minred,1.0,0.0)*(2.0**-RED2)*(13.0**-RED)
    for an,A in (('ST',ST),('KA',KA)):
        lut=np.zeros(128,int);lut[[ord(c) for c in A]]=np.arange(26)
        P=np.array([A.index(ANCH[i]) for i in AN])
        for rule in (0,1):
            for c0 in range(0,len(Ms),4000):
                M=Ms[c0:c0+4000];C=lut[cl[M[:,AN]]]
                K=((C-P)%26 if rule==0 else (C+P)%26)
                ok=np.ones((len(M),len(models)),bool)
                for m in (2,13):
                    S=np.mod((K%m).astype(np.float32)@NM[m],m)!=0
                    bad=np.zeros((len(M),len(models)),bool)
                    # any nonzero within a model's group
                    cnt=np.add.reduceat(S,np.r_[0,np.nonzero(np.diff(GM[m]))[0]+1],axis=1)
                    ks=GM[m][np.r_[0,np.nonzero(np.diff(GM[m]))[0]+1]]
                    bad[:,ks]=cnt>0;ok&=~bad
                ok&=(RED>=minred)[None,:]
                exp+=len(M)*w.sum()
                for x,k in zip(*np.nonzero(ok)): hits.append((int(RED[k]),models[k][0],an,['vig','beau'][rule],labs[c0+x]))
    return hits,exp
# ---------- permutations ----------
def dial(steps,first):
    L=list(range(N));idx=0;out=[];d=first;j=0
    while L:
        n=steps[j%len(steps)]
        if d=='R': idx=(idx+n-1)%len(L)
        else: idx=(idx-n)%len(L)
        out.append(L.pop(idx))
        if d=='L' and L: idx=idx%len(L)
        if idx>=len(L) and L: idx=0
        d='L' if d=='R' else 'R';j+=1
    return out
def dial_perms(lim):
    P=[];labs=[]
    for k,lo,hi in lim:
        for st in itertools.product(range(lo,hi+1),repeat=k):
            for f in 'LR':
                o=dial(st,f);P.append(o);labs.append(('dial',st,f))
    P=np.array(P,dtype=np.int16);inv=np.argsort(P,axis=1).astype(np.int16)
    return np.concatenate([P,inv]),labs+[l+('inv',) for l in labs]
