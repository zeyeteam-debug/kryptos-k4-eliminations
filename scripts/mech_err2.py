import numpy as np,random
exec(open('mech_err.py').read().split("SH=list")[0])
from gen import run
tot=0;n=0
for s in range(60):
    sh=list(CT0);random.Random(1000+s).shuffle(sh);h=run(''.join(sh),ANCH,P,L,Pmax=12);tot+=len(h);n+=1
print('null shuffles',n,'mean hits per run',tot/n)
