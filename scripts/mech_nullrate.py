import numpy as np,random
from mech import dial_perms,CT0,ANCH
from gen import run
Pd,Ld=dial_perms([(1,1,97),(2,1,97),(3,1,26),(4,1,12)])
tot=0
for s in range(8):
    sh=list(CT0);random.Random(500+s).shuffle(sh);sh=''.join(sh)
    for c0 in range(0,len(Pd),20000): tot+=len([h for h in run(sh,ANCH,Pd[c0:c0+20000],Ld[c0:c0+20000],Pmax=12) if h[0]=='any'])
    print('shuffle',s,'cum any-hits',tot,'perm-tests',(s+1)*len(Pd),flush=True)
