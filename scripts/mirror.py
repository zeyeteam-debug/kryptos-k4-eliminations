import numpy as np
from stack import CT0,ANCH
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
hits=[];n=0
for an,A in (('ST',ST),('KA',KA)):
    for src in ('pt','ct'):
        for c in range(40,160):          # reflection centre: partner j=c-i
            pairs=[(i,c-i) for i in ANCH if 0<=c-i<97 and (c-i in ANCH or src=='ct')]
            if len(pairs)<6: continue
            for rule in range(4):
                ok=True
                for i,j in pairs:
                    k=A.index(ANCH[j] if src=='pt' else CT0[j]);p=A.index(ANCH[i]);ct=A.index(CT0[i])
                    e=[(p+k)%26,(k-p)%26,(p-k)%26,(p+k+13)%26][rule]
                    if e!=ct: ok=False;break
                n+=1
                if ok: hits.append((an,src,c,rule,len(pairs)))
print('reflection-autokey models',n,'consistent:',hits)
