import numpy as np,time,sys
from rk_trans import prep,search,KA,ST
from stack import CT0
from gen import colperms
import cyl,counter
from k3rot import compose
s=open('texts/carter_vol1.txt',errors='ignore').read().upper();t=''.join(c for c in s if 'A'<=c<='Z')
big=t+'|'+t[::-1]
fams=[]
fams.append(('none',np.arange(97,dtype=np.int16)[None,:]))
M=[];
for w1 in range(2,97):
    for w2 in range(2,97):
        for d1 in (True,False):
            for d2 in (True,False):
                idx=compose(97,w1,d1,w2,d2);m=np.empty(97,dtype=np.int16);m[idx]=np.arange(97);M+=[m,idx.astype(np.int16)]
fams.append(('K3rot',np.array(M)))
fams.append(('cyl',cyl.family()[0]));fams.append(('rot',counter.family()[0]))
for W in range(2,10): fams.append((f'col{W}',colperms(97,W)))
for an,A in (('KA',KA),('ST',ST)):
    arr,pos=prep(big,A)
    for name,Ms in fams:
        t0=time.time();h=search(CT0,Ms,[(name,i) for i in range(len(Ms))],arr,pos,A)
        print(an,name,len(Ms),'hits',len(h),h[:5],round(time.time()-t0),'s',flush=True)
