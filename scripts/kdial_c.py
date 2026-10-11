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
SH=None
# (c) null mask: dial stops (cumulative, alternating direction on the 97 ring) mark nulls; cribs must stay untouched
crib=set(ANCH);res=[]
for k,c in combos.items():
    for f in 'LR':
        for nn in range(1,15):  # number of nulls
            pos=0;d=1 if f=='R' else -1;nulls=[];j=0;ok=True
            while len(nulls)<nn:
                pos=(pos+d*c[j%len(c)])%97;j+=1;d=-d
                if j>300: break
                if pos in nulls: continue
                nulls.append(pos)
                if j>200: break
            if len(nulls)<nn or crib&set(nulls): continue
            keep=[i for i in range(97) if i not in nulls];ct2=''.join(CT0[i] for i in keep)
            an2={keep.index(i):ch for i,ch in ANCH.items()}
            H=run(ct2,an2,np.arange(len(ct2),dtype=np.int16)[None,:],[(k,f,nn)],Pmax=24)
            res.append((len(H),k,f,nn,sorted(nulls)))
print('(c) null masks tested (crib-safe)',len(res),'with hits',[r for r in res if r[0]][:5])

# control: insert 7 nulls at crib-safe positions into a KA-vig period-7 cipher, then test
pt='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONX'[:90]
nul=[3,10,40,50,80,88,95];keep=[i for i in range(97) if i not in nul]
ptl=['?']*97
for j,i in enumerate(keep): ptl[i]=pt[j]
for i,ch in ANCH.items(): ptl[i]=ch
ctl=''.join('Q' if ptl[i]=='?' else KA[(KA.index(ptl[i])+KA.index('KRYPTOS'[keep.index(i)%7]))%26] for i in range(97))
ct2=''.join(ctl[i] for i in keep);an2={keep.index(i):ch for i,ch in ANCH.items()}
print('control',[h for h in run(ct2,an2,np.arange(90,dtype=np.int16)[None,:],['c'],Pmax=24) if h[0]=='key'][:2])
