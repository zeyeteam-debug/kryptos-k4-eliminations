# Triple 25+21+17 (Y, U, Q as wheel sizes), pins fixed by letter sets, all offsets, weights solved mod 2 / mod 13
import numpy as np,itertools,random
from hagelin import WL,PINSETS,kvals,idx,ST
from stack import CT0,ANCH,solve_mod
SUB=(25,21,17)
def run(ct,pinsets=PINSETS,neg=True):
    hits=[];t=0
    for nm,letters in pinsets.items():
        for comp in ((False,True) if neg else (False,)):
            pins={p:[int((c in letters)!=comp) for c in WL[p]] for p in SUB}
            for rule in ('beau','vig'):
                K=kvals(ct,rule);k2=list(K%2);k13=list(K%13)
                for offs in itertools.product(*[range(p) for p in SUB]):
                    rows=[[pins[p][(i+o)%p] for p,o in zip(SUB,offs)]+[1] for i in idx]
                    o13,rk=solve_mod(rows,k13,13);t+=1
                    if o13 and solve_mod(rows,k2,2)[0]: hits.append((24-rk,nm,comp,offs,rule))
    return hits,t
pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
for i,c in ANCH.items(): pt[i]=c
pins={p:[int(c in 'QUAE') for c in WL[p]] for p in SUB};o=(3,7,11)
ctl=''.join(ST[(2+5*pins[25][(i+3)%25]+9*pins[21][(i+7)%21]+13*pins[17][(i+11)%17]-ST.index(c))%26] for i,c in enumerate(pt))
h,_=run(ctl,{'QUAE':'QUAE'},neg=False);print('control',[x for x in h if x[3]==o],flush=True)
SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
for nm,ct in (('K4',CT0),('null',SH)):
    h,t=run(ct);print(nm,'tested',t,'consistent',len(h),h[:5],flush=True)
