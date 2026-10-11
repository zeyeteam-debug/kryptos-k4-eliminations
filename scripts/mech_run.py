import numpy as np,random,time,sys,pickle
import mech
from mech import *
from gen import colperms
t0=time.time()
ID=np.arange(97,dtype=np.int16)[None,:]
SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
def plant(perm,key):  # perm[pt]=ct position ; key list of 97 KA vig shifts
    pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
    for i,c in ANCH.items(): pt[i]=c
    ct=['?']*97
    for i,c in enumerate(pt): ct[perm[i]]=KA[(KA.index(c)+key[i])%26]
    return ''.join(ct)
out={}
# ---- 1) direct, all models
mech.build(mech.ALL)
for nm,ct in (('K4',CT0),('null',SH)):
    h,e=check(ct,ID,['direct']);out[('direct',nm)]=(h,e)
    print('DIRECT',nm,'hits',len(h),'expected',round(e,3),sorted(h,key=lambda x:-x[0])[:8],flush=True)
# control: clock r=12 s=5, a=3,b=7,c=11 ; gear (4,15)
kc=[(3*i+7*((i+5)//12)+11)%26 for i in range(97)];print('ctrl clock',[x for x in check(plant(list(range(97)),kc),ID,['c'])[0] if x[1]==('clock',12,5)])
g1=[random.randrange(26) for _ in range(4)];g2=[random.randrange(26) for _ in range(15)]
kg=[(g1[i%4]+g2[i%15])%26 for i in range(97)];print('ctrl gear',[x for x in check(plant(list(range(97)),kg),ID,['c'])[0] if x[1]==('gear',4,15)])
# ---- 1b) one-error tolerance, direct: drop each crib in turn
from stack import solve_mod
tol=[];RED0=list(mech.RED)
for (lab,A),rd in zip(mech.ALL,RED0):
    if rd<5: continue
    for an,AL in (('ST',ST),('KA',KA)):
        P=[AL.index(ANCH[i]) for i in AN];C=[AL.index(CT0[i]) for i in AN]
        for rule in (0,1):
            K=[((c-p) if rule==0 else (c+p))%26 for c,p in zip(C,P)]
            for drop in range(24):
                rows=[[int(v) for v in A[j]] for j in range(24) if j!=drop];rhs=[K[j] for j in range(24) if j!=drop]
                o2,_=solve_mod(rows,rhs,2);o13,r=solve_mod(rows,rhs,13)
                if o2 and o13 and 23-r>=4: tol.append((23-r,lab,an,rule,int(AN[drop])))
print('ONE-ERROR direct (red>=4):',len(tol),sorted(tol,key=lambda x:-x[0])[:8],flush=True)
# ---- 2) transpositions, strong models only
keep_r={12,24,60}|{x for t in special.values() for x in t if x>=2}
sel=[(l,A) for (l,A) in mech.ALL if (l[0]=='clock' and l[1] in keep_r) or l[0]=='clock3' or l[0]=='gear']
mech.build(sel);sel=[m for m,r in zip(sel,mech.RED) if r>=6];mech.build(sel)
print('transposition models',len(sel),'null cols',mech.NM[13].shape,flush=True)
Pd,Ld=dial_perms([(1,1,97),(2,1,97),(3,1,26),(4,1,12)])
print('dial perms',len(Pd),round(time.time()-t0),flush=True)
Pc=[];Lc=[]
for W in range(2,9):
    M=colperms(97,W);Pc.append(M);Lc+= [('col',W,k) for k in range(len(M))]
Pc=np.concatenate(Pc)
# controls inside dial family
pd=dial((5,9,3),'L');key=[(3*i+7*((i+2)//24)+4)%26 for i in range(97)]
print('ctrl dial+clock',[x for x in check(plant(pd,key),np.array([pd],dtype=np.int16),['ctrl'])[0] if x[1]==('clock',24,2)],flush=True)
pd2=dial((24,1,17),'R');g=[random.randrange(26) for _ in range(5)]
selg=[x for x in sel if x[0]==('gear',2,5)]
for fam,Ps,Ls in (('dial',Pd,Ld),('columnar2-8',Pc,Lc)):
    for nm,ct in (('K4',CT0),('null',SH)):
        h,e=check(ct,Ps,Ls,minred=6);out[(fam,nm)]=(h,e)
        print(fam,nm,'hits',len(h),'expected',round(e,3),sorted(h,key=lambda x:-x[0])[:6],round(time.time()-t0),flush=True)
pickle.dump((out,tol),open('mech_out.pkl','wb'))
