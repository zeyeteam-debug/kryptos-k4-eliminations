import numpy as np,random
from stack import CT0,ANCH
from cyl import run
import cyl
N=97
def from_groups(groups,shifts):
    # groups: lists of positions; each group's letters rotated cyclically by shift (positive=right)
    m=np.empty(N,dtype=np.int16)        # m[pt]=ct
    for g,s in zip(groups,shifts):
        L=len(g)
        for j,p in enumerate(g): m[p]=g[(j+s)%L]
    return m
def family():
    Ms=[];labs=[]
    def add(m,lab):
        Ms.append(m);labs.append(lab+('enc',));Ms.append(np.argsort(m).astype(np.int16));labs.append(lab+('dec',))
    # (1) two pieces split at h, rotated a and b
    for h in range(30,68):
        A=list(range(h));B=list(range(h,N))
        for a in range(len(A)):
            for b in range(len(B)):
                add(from_groups([A,B],[a,b]),('halves',h,a,b))
    # (2) even/odd interleave
    E=list(range(0,N,2));O=list(range(1,N,2))
    for a in range(len(E)):
        for b in range(len(O)): add(from_groups([E,O],[a,b]),('evenodd',a,b))
    # (3) grid rows/cols alternately rotated (gears)
    for W in range(2,25):
        rows=[list(range(r*W,min(N,(r+1)*W))) for r in range((N+W-1)//W)]
        cols=[list(range(c,N,W)) for c in range(W)]
        for kind,grp in (('rows',rows),('cols',cols)):
            Lm=max(len(g) for g in grp)
            for a in range(Lm):
                for b in range(Lm):
                    add(from_groups(grp,[a if i%2==0 else -b for i in range(len(grp))]),(kind,W,a,b))
    return np.array(Ms),labs
if __name__=='__main__':
    Ms,labs=family();print('permutations',len(Ms),flush=True)
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
    KA=cyl.KA;m=from_groups([list(range(48)),list(range(48,97))],[-7,3]);key='CLOCK';ct=['?']*N
    for i in range(N): ct[m[i]]=KA[(KA.index(pt[i])+KA.index(key[i%5]))%26]
    out=[];CH=60000
    def runall(ct):
        h=[]
        for s in range(0,len(Ms),CH): h+=run(ct,Ms[s:s+CH],labs[s:s+CH])
        return h
    h=runall(''.join(ct));print('control (halves 48, 7 left, 3 right, key CLOCK):',[x for x in h if x[4][:4]==('halves',48,41,3)][:2],'total',len(h),flush=True)
    random.seed(4);s=list(CT0);random.shuffle(s);print('null hits',len(runall(''.join(s))),flush=True)
    h=runall(CT0);print('K4 hits',len(h),h[:20])
