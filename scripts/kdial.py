# KRYPTOS as a safe combination: (a) dial transposition, (b) dial-sweep keystream, (c) null mask from dial stops
import numpy as np,itertools,random
import mech
from mech import *
from gen import run
W='KRYPTOS'
enc={'ST0':[ST.index(c) for c in W],'ST1':[ST.index(c)+1 for c in W],'KA1':[KA.index(c)+1 for c in W],
     'dST':[(ST.index(b)-ST.index(a))%26 for a,b in zip(W,W[1:]+W[0])],
     'dSTabs':[min((ST.index(b)-ST.index(a))%26,(ST.index(a)-ST.index(b))%26) for a,b in zip(W,W[1:]+W[0])],
     'TAB':[ST.index(c)+1 for c in 'KRYPTOSABCDEFGHIJLMNQUVWXZ'][:7]}
combos={}
for nm,v in enc.items():
    v=[x if x>0 else 26 for x in v]
    for rv in (0,1):
        b=v[::-1] if rv else v
        for t in range(4): combos[(nm,rv,'uni',t)]=tuple(x+26*t for x in b)
        combos[(nm,rv,'safe',0)]=tuple(x+26*(len(b)-1-i) for i,x in enumerate(b))   # first number most turns
        combos[(nm,rv,'safe',1)]=tuple(x+26*(len(b)-i) for i,x in enumerate(b))
def dial_abs(targets,first):  # absolute stops: pointer goes to target index (mod remaining) alternating direction (direction irrelevant for absolute; vary origin)
    L=list(range(97));out=[];j=0
    while L:
        idx=targets[j%len(targets)]%len(L);out.append(L.pop(idx));j+=1
    return out
P=[];LB=[]
for k,c in combos.items():
    for f in 'LR': P.append(dial(c,f));LB.append(('step',k,f))
    P.append(dial_abs(c,'L'));LB.append(('abs',k))
P=np.array(P,dtype=np.int16);P=np.concatenate([P,np.argsort(P,axis=1).astype(np.int16)]);LB=LB+[l+('inv',) for l in LB]
SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
print('(a) perms',len(P))
for nm,ct in (('K4',CT0),('null',SH)):
    H=run(ct,ANCH,P,LB,Pmax=12)
    keep_r={12,24,60}|{x for t in special.values() for x in t if x>=2}
    sel=[(l,A) for (l,A) in mech.ALL if (l[0]=='clock' and l[1] in keep_r) or l[0] in('clock3','gear')]
    mech.build(sel);s2=[m for m,r in zip(sel,mech.RED) if r>=6];mech.build(s2)
    h2,e=check(ct,P,LB,minred=6)
    print(' ',nm,'periodic/any-alpha hits',len(H),H[:3],'| clock/gear hits',len(h2),'exp',round(e,4))
# (b) dial-sweep keystream: pointer clicks one letter at a time on a KA or ST ring, going to each letter of KRYPTOS
# alternately R/L with t extra full turns; keystream = letters under pointer at each click (repeat combo cyclically)
def sweep(ring,word,turns,first,start):
    pos=start;ks=[];d=1 if first=='R' else -1;j=0
    while len(ks)<97+30:
        tgt=ring.index(word[j%len(word)]);t=turns[j%len(turns)]
        steps=((tgt-pos)*d)%26+26*t
        if steps==0: steps=26
        for _ in range(steps): pos=(pos+d)%26;ks.append(pos)
        d=-d;j+=1
    return ks
best=[];n=0
for rn,ring in (('ST',ST),('KA',KA)):
    for word in (W,W[::-1]):
        for turns in [(t,) for t in range(4)]+[(4,3,2,1),(3,2,1,0),(6,5,4,3,2,1,0),(1,0)]:
            for first in 'LR':
                for start in range(26):
                    ks=sweep(ring,word,turns,first,start)
                    for off in range(0,30):
                        for an,A in (('ST',ST),('KA',KA)):
                            for rule in range(3):
                                m=0
                                for i,p in ANCH.items():
                                    k=A.index(ring[ks[i+off]]);c=A.index(CT0[i]);pp=A.index(p)
                                    m+=c==[(pp+k)%26,(k-pp)%26,(pp-k)%26][rule]
                                n+=1;best.append((m,rn,word,turns,first,start,off,an,rule))
best.sort(key=lambda x:-x[0]);print('(b) sweep keystreams tested',n,'best',best[:3])
# (c) null mask: dial stops (cumulative, alternating direction on the 97 ring) mark nulls; cribs must stay untouched
crib=set(ANCH);res=[]
for k,c in combos.items():
    for f in 'LR':
        for nn in range(1,15):  # number of nulls
            pos=0;d=1 if f=='R' else -1;nulls=[];j=0;ok=True
            while len(nulls)<nn:
                pos=(pos+d*c[j%len(c)])%97;j+=1;d=-d
                if pos in nulls: continue
                nulls.append(pos)
                if j>200: break
            if len(nulls)<nn or crib&set(nulls): continue
            keep=[i for i in range(97) if i not in nulls];ct2=''.join(CT0[i] for i in keep)
            an2={keep.index(i):ch for i,ch in ANCH.items()}
            H=run(ct2,an2,np.arange(len(ct2),dtype=np.int16)[None,:],[(k,f,nn)],Pmax=24)
            res.append((len(H),k,f,nn,sorted(nulls)))
print('(c) null masks tested (crib-safe)',len(res),'with hits',[r for r in res if r[0]][:5])
