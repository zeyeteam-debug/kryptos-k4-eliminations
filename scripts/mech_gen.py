import numpy as np,random,time
from mech import dial_perms,dial,CT0,ANCH,KA
from gen import run
Pd,Ld=dial_perms([(1,1,97),(2,1,97),(3,1,26),(4,1,12)])
SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
# control: dial (7,3,11) + KA vig period 7 key KRYPTOS
pd=dial((7,3,11),'R');pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
for i,c in ANCH.items(): pt[i]=c
ct=['?']*97
for i,c in enumerate(pt): ct[pd[i]]=KA[(KA.index(c)+KA.index('ABSCISS'[i%7]))%26]
print('control',[h for h in run(''.join(ct),ANCH,np.array([pd]),['c'],Pmax=12) if h[0]=='key'][:3])
for nm,ct in (('K4',CT0),('null',SH)):
    H=[]
    for c0 in range(0,len(Pd),20000):
        H+=run(ct,ANCH,Pd[c0:c0+20000],Ld[c0:c0+20000],Pmax=12)
    print(nm,'periodic/any-alphabet hits',len(H),H[:6],flush=True)
