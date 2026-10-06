import numpy as np,random
from k3rot import compose
from gen import run
from stack import CT0,ANCH
N=97;KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
Ms=[];labs=[]
for w1 in range(2,97):
    for w2 in range(2,97):
        for d1 in (True,False):
            for d2 in (True,False):
                idx=compose(N,w1,d1,w2,d2)       # ct[n]=pt[idx[n]]
                m=np.empty(N,dtype=np.int16);m[idx]=np.arange(N)   # m[pt]=ct
                Ms.append(m);labs.append((w1,d1,w2,d2,'enc'))
                Ms.append(idx.astype(np.int16));labs.append((w1,d1,w2,d2,'dec'))
Ms=np.array(Ms);print('perms',len(Ms),'distinct',len({m.tobytes() for m in Ms}),flush=True)
def runall(ct):
    h=[]
    for s in range(0,len(Ms),40000): h+=run(ct,ANCH,Ms[s:s+40000],labs[s:s+40000])
    return h
src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
idx=compose(N,8,True,13,True);x=''.join(pt[i] for i in idx)
ct=''.join(KA[(KA.index(c)+KA.index('ABSCISSA'[i%8]))%26] for i,c in enumerate(x))
h=runall(ct);print('control (8,13 + key ABSCISSA after):',[t for t in h if t[-1][:4]==(8,True,13,True)][:2],'total',len(h),flush=True)
l=list(CT0);random.Random(2).shuffle(l);print('null hits',len(runall(''.join(l))),flush=True)
hk=runall(CT0);print('K4 hits',len(hk),hk[:10])
