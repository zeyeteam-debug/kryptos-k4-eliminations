from stack import CT0,ANCH
import numpy as np
from gen import run
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
Wp=[i for i,c in enumerate(CT0) if c=='W']
starts=[0]+[w+1 for w in Wp]
def segstart(i): return max(s for s in starts if s<=i)
print('segments start at',starts)
loc={i:i-segstart(i) for i in ANCH}
print('EASTNORTHEAST local offsets',[loc[i] for i in range(21,34)]);print('BERLINCLOCK local offsets',[loc[i] for i in range(63,74)])
# 1) ANY key (even non-repeating) that restarts after each W: same local offset -> same key value
pairs=[(a,b) for a in range(21,34) for b in range(63,74) if loc[a]==loc[b]]
print('overlapping local offsets:',[(loc[a],ANCH[a],CT0[a],ANCH[b],CT0[b]) for a,b in pairs])
for an,A in (('ST',ST),('KA',KA)):
    for r,rn in ((0,'vig'),(1,'beau')):
        k=lambda i:((A.index(CT0[i])-A.index(ANCH[i])) if r==0 else (A.index(CT0[i])+A.index(ANCH[i])))%26
        agree=sum(k(a)==k(b) for a,b in pairs);print(f'{an} {rn}: restarting key agrees at {agree}/{len(pairs)} overlaps')
bad=sum((ANCH[a]==ANCH[b])!=(CT0[a]==CT0[b]) for a,b in pairs)
print('any alphabet, any restarting key: contradictions',bad,'(strong checks:',sum(ANCH[a]==ANCH[b] for a,b in pairs),')')
