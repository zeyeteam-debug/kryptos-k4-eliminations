# Keyboard layouts as alphabets and as routes
import numpy as np,random,time
import allalpha as AA
from stack import CT0,ANCH
from gen import run
KB={'QWERTY':'QWERTYUIOPASDFGHJKLZXCVBNM','QWERTZ':'QWERTZUIOPASDFGHJKLYXCVBNM','AZERTY':'AZERTYUIOPQSDFGHJKLMWXCVBN','DVORAK':'PYFGCRLAOEUIDHTNSQJKXBMWVZ',
    'QWERTYcol':'QAZWSXEDCRFVTGBYHNUJMIKOLP','QWERTZcol':'QAYWSXEDCRFVTGBZHNUJMIKOLP'}
ALS={}
for k,v in KB.items(): ALS[k]=v;ALS[k+'rev']=v[::-1]
AA.ALS=ALS
t=time.time();SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
M=AA.linear_models()
for nm,ct in (('K4',CT0),('null',SH)):
    labs,K=AA.Kmat(ct)
    print(nm,'linear',len(AA.linear_test(K,M)),'safe',len(AA.safe_test(K)),'masks',len(AA.mask_test(K)[0]),'hill',len(AA.hill_test(ct)),'autokey',len(AA.autokey_test(ct)),round(time.time()-t),flush=True)
T=AA.texts()
for nm,ct in (('K4',CT0),('null',SH)): print(nm,'running key',len(AA.rk_test(ct,T)),flush=True)
import dictsearch as D
sp=[('none',np.arange(97))]
for Wd in range(2,49):
    for corner in ('TL','TR','BL','BR'):
        for cw in (0,1):
            m=D.perm_from_read(D.spiral_order(97,Wd,corner,cw),97);sp+=[((Wd,corner,cw,'r'),m),((Wd,corner,cw,'w'),np.argsort(m))]
tot=0
for nm,A in ALS.items():
    TB=D.table(D.WORDS,ANCH,A);h=D.run(CT0,ANCH,sp,[],A,TB)+D.run(CT0,ANCH,[('none',np.arange(97))],[2,3,4,5,6,7],A,TB);tot+=len(h)
print('PK9-style dictionary pipeline, keyboard alphabets: hits',tot,round(time.time()-t),flush=True)
# keyboard geometry as a route: K4 written in rows of 10,9,7 (the three letter rows), repeated; read by columns (aligned or staggered)
routes=[];labs=[]
for stagger in (0,1):
    for order in ((10,9,7),(7,9,10)):
        rows=[];i=0;r=0
        while i<97:
            L=order[r%3];rows.append(list(range(i,min(97,i+L))));i+=L;r+=1
        off=[0,1,2] if stagger else [0,0,0]
        cells=[]
        for ri,row in enumerate(rows):
            for j,x in enumerate(row): cells.append((j+off[ri%3],ri,x))
        for colfirst in (0,1):
            seq=[x for _,_,x in sorted(cells,key=(lambda z:(z[0],z[1])) if colfirst else (lambda z:(z[0],-z[1])))]
            m=np.array(seq);routes+=[m,np.argsort(m)];labs+=[(stagger,order,colfirst,'a'),(stagger,order,colfirst,'b')]
P=np.array(routes,dtype=np.int16)
for nm,ct in (('K4',CT0),('null',SH)): print(nm,'keyboard-row routes + periodic/any-alpha:',len(run(ct,ANCH,P,labs,Pmax=24)),flush=True)
