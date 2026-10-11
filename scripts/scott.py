import numpy as np,random
from stack import CT0,ANCH
from gen import run
KA='KRYPTOSABCDEFGHIJLMNQUVWXZ';ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';N=97
def q3(ct,key,sign):
    return ''.join(KA[(KA.index(c)-sign*KA.index(key[i%len(key)]))%26] for i,c in enumerate(ct))
def aff(s,a,b,direction):
    out=['?']*N
    for x in range(N):
        y=(a+b*x)%N
        if direction==0: out[y]=s[x]
        else: out[x]=s[y]
    return ''.join(out)
def score(p): return sum(p[i]==c for i,c in ANCH.items())
def show(p): return p[21:34]+' | '+p[63:74]
for sign in (1,-1):
    for direction in (0,1):
        for key in ('EMUFPHZLRFA',):
            s1=q3(CT0,key,sign);p=aff(s1,77,38,direction)
            print('sign',sign,'dir',direction,'score',score(p),show(p),'| BERLIN at 63?',p[63:69]=='BERLIN', 'full:',p)
# 3-step variant from the post
for sign in (1,-1):
    for direction in (0,1):
        p=aff(q3(q3(CT0,'EMUFPHZLRFA',sign),'HXDKZEK',sign),77,38,direction);print('3-step',sign,direction,score(p),show(p))
