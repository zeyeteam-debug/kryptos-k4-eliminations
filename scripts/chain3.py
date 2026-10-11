# "Three turns of the dial": K4 = chain of K1's method (KA Vigenere, PALIMPSEST), K2's (KA Vigenere, ABSCISSA),
# K3's (double rotation, any widths) in any order, each turn either direction, key phase 0 or continued from K1/K2 lengths.
import numpy as np,itertools,random,time
from k3rot import compose
from stack import CT0,ANCH
KA='KRYPTOSABCDEFGHIJLMNQUVWXZ';N=97
AN=np.array(sorted(ANCH));PT=np.array([KA.index(ANCH[i]) for i in AN])
Ms=[];labs=[]
for w1 in range(2,97):
    for w2 in range(2,97):
        for d1 in (True,False):
            for d2 in (True,False):
                idx=compose(N,w1,d1,w2,d2);Ms.append(idx);labs.append((w1,d1,w2,d2,'a'))
                inv=np.empty(N,int);inv[idx]=np.arange(N);Ms.append(inv);labs.append((w1,d1,w2,d2,'b'))
Ms=np.array(Ms,dtype=np.int16)
def keyarr(w,off): return np.array([KA.index(w[(i+off)%len(w)]) for i in range(N)])
KEYS={'P':[('PALIMPSEST',0),('PALIMPSEST',63%10)],'A':[('ABSCISSA',0),('ABSCISSA',369%8)]}
def search(ct,report=8):
    T0=np.array([KA.index(c) for c in ct]);res=[]
    for order in itertools.permutations('PAT'):
        for sP in (1,-1):
            for sA in (1,-1):
                for kp,ka in itertools.product(KEYS['P'],KEYS['A']):
                    X=np.broadcast_to(T0,(len(Ms),N)).copy()
                    for op in order:
                        if op=='T': X=np.take_along_axis(X,Ms.astype(np.int64),axis=1)
                        elif op=='P': X=(X+sP*keyarr(*kp))%26
                        else: X=(X+sA*keyarr(*ka))%26
                    sc=(X[:,AN]==PT).sum(axis=1);b=np.argsort(-sc)[:3]
                    for j in b: res.append((int(sc[j]),''.join(order),sP,sA,kp,ka,labs[j]))
    res.sort(key=lambda x:-x[0]);return res[:report]
if __name__=='__main__':
    t=time.time()
    # positive control: pt -> Vig ABSCISSA -> rotation(13,7) -> Vig PALIMPSEST (forward), decoded by the reverse chain
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK'
    x=np.array([KA.index(c) for c in pt]);x=(x+keyarr('ABSCISSA',0))%26
    idx=compose(N,13,True,7,False);x=x[idx];x=(x+keyarr('PALIMPSEST',3))%26
    ctl=''.join(KA[v] for v in x)
    print('control',search(ctl,2),round(time.time()-t),flush=True)
    print('K4',search(CT0),flush=True)
    l=list(CT0);random.Random(2).shuffle(l);print('null',search(''.join(l),3),round(time.time()-t),flush=True)
