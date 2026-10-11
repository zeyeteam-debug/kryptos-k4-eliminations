# PK9 structure with the two Berlin clocks replacing the spiral:
#   Q(any of 370k dictionary words, pt-indexed) -> clock route -> columnar 2..7 (all orders), either order of route/columnar
# Clock routes: Weltzeituhr cylinder (24 faces: write rows of 24, read faces from start face s, east/west, up/down),
#               helical readings (decimation k mod 97, all k), Mengenlehreuhr band orders on 24-lamp blocks (1+4+4+11+4).
import numpy as np,itertools,time
import dictsearch as D
from stack import CT0,ANCH
N=97
def cyl_routes():
    R=[]
    for s in range(24):
        for d in (1,-1):
            for ud in (0,1):
                idx=[]
                for j in range(24):
                    c=(s+d*j)%24;col=list(range(c,N,24))
                    idx+=col[::-1] if ud else col
                R.append((('cyl',s,d,ud),idx))
    return R
def helix_routes():
    return [(('helix',k,st),[(st+k*j)%N for j in range(N)]) for k in range(2,N) for st in ((0,) if k not in (23,24,25) else range(N))]
BANDS=[[0],[1,2,3,4],[5,6,7,8],list(range(9,20)),[20,21,22,23]]
def mengen_routes():
    R=[]
    for perm in itertools.permutations(range(5)):
        for rev in range(32):
            order=[]
            for b,bi in enumerate(perm):
                band=BANDS[bi];order+=band[::-1] if (rev>>b)&1 else band
            idx=[]
            for s in range(0,N,24):
                idx+=[s+o for o in order if s+o<N]
            R.append((('mengen',perm,rev),idx))
    return R
if __name__=='__main__':
    t=time.time()
    routes=[]
    for lab,idx in cyl_routes()+helix_routes()+mengen_routes():
        m=D.perm_from_read(idx,N);routes+=[(lab+('r',),m),(lab+('w',),np.argsort(m))]
    # dedupe
    seen=set();R2=[]
    for lab,m in routes:
        k=m.tobytes()
        if k in seen: continue
        seen.add(k);R2.append((lab,m))
    print('distinct clock routes',len(R2),flush=True)
    # positive control: PK9-like with a clock route: Q(CLEPSYDRA) -> cylinder route -> columnar width 6
    import random
    KA=D.KA;pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
    for i,c in ANCH.items(): pt[i]=c
    x=[KA[(KA.index(c)+KA.index('CLEPSYDRA'[i%9]))%26] for i,c in enumerate(pt)]
    lab,m=[r for r in R2 if r[0][0]=='cyl'][5];y=['?']*N
    for i in range(N): y[m[i]]=x[i]
    M6=D.colperms_all(N,6);c=M6[100];z=['?']*N
    for i in range(N): z[c[i]]=y[i]
    ctl=''.join(z)
    TB=D.table(D.WORDS,ANCH,KA);print('table',round(time.time()-t),flush=True)
    h=D.run(ctl,ANCH,[(lab,m)],[6],KA,TB);print('control',[x for x in h if x[0]=='CLEPSYDRA'][:2],flush=True)
    for an,AL in (('KA',D.KA),('ST',D.ST)):
        TB=D.table(D.WORDS,ANCH,AL)
        h=D.run(CT0,ANCH,R2,[2,3,4,5,6,7],AL,TB);print(an,'K4 hits',len(h),h[:5],round(time.time()-t),flush=True)
