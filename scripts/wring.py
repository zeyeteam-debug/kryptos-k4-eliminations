# (B) hour-ring alignment: key = a*h_i + c, h_i = (s*i + t) mod 24 (local hour on the face that letter i sits on); and city-count keys
import numpy as np
from stack import CT0,ANCH,solve_mod
from wcity import F,keys
AN=sorted(ANCH);K=keys(CT0);n=0;hits=[]
def lin(vals,k):
    rows=[[v%26,1] for v in vals];o2=solve_mod(rows,[x%2 for x in k],2)[0];o13,r=solve_mod(rows,[x%13 for x in k],13);return o2 and o13,24-r
for s in range(1,24):
    for t in range(24):
        for base in (0,1):
            v=[(s*i+t)%24+base for i in AN]
            for kk,k in K.items():
                ok,red=lin(v,list(k));n+=1
                if ok and red>=3: hits.append(('ring',s,t,base,kk,red))
tc=[len(F[z][0].split()) for z in range(-8,13)];bc=[len(F[z][1].split()) for z in range(-8,13)]
for nm,c in (('top',tc),('bottom',bc),('both',[a+b for a,b in zip(tc,bc)])):
    for o in range(len(c)):
        for d in (1,-1):
            v=[c[(o+d*i)%len(c)] for i in AN]
            for kk,k in K.items():
                ok,red=lin(v,list(k));n+=1
                if ok and red>=3: hits.append(('count',nm,o,d,kk,red))
print('(B) ring/count models tested',n,'consistent',hits[:5])
