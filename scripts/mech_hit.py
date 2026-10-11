import numpy as np
from mech import dial,CT0,ANCH
o=dial((86,82,89,92),'L');inv=np.argsort(o)  # label 'inv' => M = argsort(order)
M=inv
for Per in (8,16,24):
    cls={}
    for i in sorted(ANCH):
        c=M[i];cls.setdefault(c%Per,[]).append((i,ANCH[i],CT0[c]))
    same_pt=sum(1 for v in cls.values() for a in range(len(v)) for b in range(a+1,len(v)) if v[a][1]==v[b][1])
    same_ct=sum(1 for v in cls.values() for a in range(len(v)) for b in range(a+1,len(v)) if v[a][2]==v[b][2])
    bad=sum(1 for v in cls.values() for a in range(len(v)) for b in range(a+1,len(v)) if (v[a][1]==v[b][1])!=(v[a][2]==v[b][2]))
    print(Per,'pairs same pt',same_pt,'same ct',same_ct,'violations',bad)
    if Per==8:
        for k in sorted(cls): print(k,[(p,c) for i,p,c in cls[k]])
