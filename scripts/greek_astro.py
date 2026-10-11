# Greek astronomy / Antikythera numbers
import numpy as np,itertools,random
from stack import CT0,ANCH,solve_mod
from greek import ALS,ST,KA
from gen import run, colperms
AN=sorted(ANCH);SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
# (a) long single periods 25..60 (incl. 38 = 2x19 = an Antikythera gear), pt-indexed, any alphabet + ST/KA/GR
ID=np.arange(97,dtype=np.int16)[None,:]
for nm,ct in (('K4',CT0),('null',SH)):
    H=run(ct,ANCH,ID,['direct'],Pmax=60,RTHR=2,ANYSTRONG=2);H=[h for h in H if (h[0]=='key' and h[4]>=25) or (h[0]=='any' and h[2]>=25)]
    print('(a)',nm,'periods 25-60 consistent:',H[:12])
# also period 38 under alphabets GR (gen.run only does ST/KA): do manually
def per_ok(ct,P,AL,rule):
    d={}
    for i in AN:
        k=((AL.index(ct[i])-AL.index(ANCH[i])) if rule==0 else (AL.index(ct[i])+AL.index(ANCH[i])))%26
        if d.setdefault(i%P,k)!=k: return False
    return True
for an,AL in ALS.items():
    print('(a2)',an,'periods 25-60 consistent (vig/beau):',[(P,r) for P in range(25,61) for r in (0,1) if per_ok(CT0,P,AL,r) and len(AN)-len({i%P for i in AN})>=3])
# (b) calendar-gear keys: k_i = a*i + b*floor((p*i+s)/q) + c  with Greek astronomical ratios
ratios={'Metonic 235/19':(235,19),'19/235':(19,235),'sidereal 254/19':(254,19),'Saros 223/18':(223,18),'draconic 242/223':(242,223),
 'anomalistic 239/223':(239,223),'Callippic 940/76':(940,76),'exeligmos 669/54':(669,54),'Olympiad 4':(1,4),'Games 2':(1,2),
 'gear 38/48':(38,48),'gear 48/24':(48,24),'gear 64/38':(64,38),'gear 127/32':(127,32),'gear 53/96':(53,96),'gear 15/27':(15,27)}
res=[];n=0
for nm,(p,q) in ratios.items():
    for s in range(q):
        f=[(p*i+s)//q for i in AN]
        rows=[[i%26,f_%26,1] for i,f_ in zip(AN,f)]
        for an,AL in ALS.items():
            for rule in (0,1):
                for tn,ct in (('K4',CT0),('null',SH)):
                    K=[((AL.index(ct[i])-AL.index(ANCH[i])) if rule==0 else (AL.index(ct[i])+AL.index(ANCH[i])))%26 for i in AN]
                    o2,_=solve_mod(rows,[k%2 for k in K],2);o13,r=solve_mod(rows,[k%13 for k in K],13);n+=1
                    if o2 and o13 and 24-r>=3: res.append((tn,24-r,nm,s,an,rule))
print('(b) calendar-gear keys tested',n,'consistent K4:',[x for x in res if x[0]=='K4'][:5],'null:',len([x for x in res if x[0]=='null']))
# control for (b)
pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
for i,c in ANCH.items(): pt[i]=c
ctl=''.join(ST[(ST.index(c)+5*i+11*((235*i+7)//19)+3)%26] for i,c in enumerate(pt))
f=[(235*i+7)//19 for i in AN];rows=[[i%26,x%26,1] for i,x in zip(AN,f)];K=[(ST.index(ctl[i])-ST.index(ANCH[i]))%26 for i in AN]
print('(b) control',solve_mod(rows,[k%2 for k in K],2)[0] and solve_mod(rows,[k%13 for k in K],13)[0])
