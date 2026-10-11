import numpy as np,itertools,random
from stack import ANCH
from frac3 import analyse,colperm,N
base=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
for i,c in ANCH.items(): base[i]=c
found=0;ok=0;tries=0
for seed in range(20000):
    rng=random.Random(seed);cells=list(itertools.product(range(3),repeat=3));rng.shuffle(cells)
    cube={chr(65+i):cells[i] for i in range(26)};inv={v:k for k,v in cube.items()}
    for mode in ('read','write'):
        w=rng.randrange(3,8);order=tuple(rng.sample(range(w),w))
        s=[t for c in base for t in cube[c]];idx=colperm(3*N,order);src=idx if mode=='read' else list(np.argsort(idx))
        o=[s[src[k]] for k in range(3*N)];c=''.join(inv.get(tuple(o[3*j:3*j+3]),'?') for j in range(N));tries+=1
        if '?' in c: continue
        found+=1;ok+=not analyse(c,w,order,mode)[1]
    if found>=20: break
print('random cubes/orders tried',tries,'| letter-only outputs',found,'| true order survives in',ok,'of',found)
