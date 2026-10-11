# K1-K3 sizes / keys as the dial combination (alternating-direction dial transposition)
import numpy as np,itertools,random
import mech
from mech import *
from gen import run
base={'len':(63,369,336),'len?':(63,372,337),'keys':(10,8),'keys+K3':(10,8,24,8),'keys+K3b':(10,8,8,24),'len+K4':(63,369,336,97)}
combos=set()
for nm,t in base.items():
    for v in (t,tuple(x%97 or 97 for x in t),tuple(x%26 or 26 for x in t)):
        for p in itertools.permutations(v): combos.add(p)
P=[];L=[]
for c in combos:
    for f in 'LR': P.append(dial(c,f));L.append((c,f))
P=np.array(P,dtype=np.int16);P=np.concatenate([P,np.argsort(P,axis=1).astype(np.int16)]);L=L+[l+('inv',) for l in L]
keep_r={12,24,60}|{x for t in special.values() for x in t if x>=2}
sel=[(l,A) for (l,A) in mech.ALL if (l[0]=='clock' and l[1] in keep_r) or l[0] in('clock3','gear')]
mech.build(sel);s2=[m for m,r in zip(sel,mech.RED) if r>=6];mech.build(s2)
SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
for nm,ct in (('K4',CT0),('null',SH)):
    H=run(ct,ANCH,P,L,Pmax=24);h2,e=check(ct,P,L,minred=6)
    print(nm,len(P),'perms | periodic/any-alpha hits',len(H),H[:3],'| clock/gear',len(h2))
