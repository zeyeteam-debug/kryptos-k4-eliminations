# dial combinations taken literally from YAR and from the misspelling letters (incl. values >12, orders permuted)
import numpy as np,itertools,random
import mech
from mech import *
keep_r={12,24,60}|{x for t in special.values() for x in t if x>=2}
sel=[(l,A) for (l,A) in mech.ALL if (l[0]=='clock' and l[1] in keep_r) or l[0] in('clock3','gear')]
mech.build(sel);sel=[m for m,r in zip(sel,mech.RED) if r>=6];mech.build(sel)
combos=set()
for t in special.values():
    t=tuple(x if x>0 else 26 for x in t)
    for p in itertools.permutations(t):
        combos.add(p)
        for k in range(1,4): combos.add(tuple(x+26*k for x in p))  # extra full turns
P=[];L=[]
for st in combos:
    for f in 'LR':
        o=dial(st,f);P.append(o);L.append((st,f))
P=np.array(P,dtype=np.int16);P=np.concatenate([P,np.argsort(P,axis=1).astype(np.int16)]);L=L+[l+('inv',) for l in L]
SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
for nm,ct in (('K4',CT0),('null',SH)):
    h,e=check(ct,P,L,minred=6);print('YAR/error dials',nm,len(P),'perms hits',len(h),'expected',round(e,4),h[:5])
from gen import run
for nm,ct in (('K4',CT0),('null',SH)):
    H=run(ct,ANCH,P,L,Pmax=12);print('YAR/error dials periodic/any-alpha',nm,len(H),H[:4])
