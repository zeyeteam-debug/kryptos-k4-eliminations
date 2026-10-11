import numpy as np,itertools,random
from stack import CT0,ANCH
from frac3 import analyse,colperm,N
rng=random.Random(4);cells=list(itertools.product(range(3),repeat=3));rng.shuffle(cells)
cube={chr(65+i):cells[i] for i in range(26)};inv={v:k for k,v in cube.items()}
pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
for i,c in ANCH.items(): pt[i]=c
def enc(pt,order,mode):
    s=[t for c in pt for t in cube[c]];idx=colperm(3*N,order)
    src=idx if mode=='read' else list(np.argsort(idx))
    o=[s[src[k]] for k in range(3*N)];out=''
    for j in range(N):
        tr=tuple(o[3*j:3*j+3]);out+=inv.get(tr,'?')
    return out
ok=0;n=0
for mode in ('read','write'):
    for order in [(2,0,3,1),(4,1,0,3,2),(1,0,2)]:
        c=enc(pt,order,mode)
        if '?' in c: continue
        n+=1;full,contra=analyse(c,len(order),order,mode);ok+=not contra
print('controls (true order) not eliminated:',ok,'of',n)
SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
for mode in ('read','write'):
    e=sum(analyse(SH,w,o,mode)[1] for w in range(2,8) for o in itertools.permutations(range(w)))
    print('null',mode,'eliminated',e,'of 5912')
