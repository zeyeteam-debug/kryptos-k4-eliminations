import numpy as np,random
from mech import dial,CT0,ANCH
from stack import solve_mod
def quag(M,ct,Per,mode):
    eqs=[];idx=sorted(ANCH)
    L=sorted(set(ANCH.values())|set(ct[M[i]] for i in idx))
    var={}
    def v(n):
        if n not in var: var[n]=len(var)
        return var[n]
    rows=[]
    for i in idx:
        p=ANCH[i];c=ct[M[i]];r={}
        a=v(('y',c)) if mode=='IV' else v(('x',c));b=v(('x',p));k=v(('k',M[i]%Per))
        r[a]=r.get(a,0)+1;r[b]=r.get(b,0)-1;r[k]=r.get(k,0)-1;rows.append(r)
    n=len(var);R=[[r.get(j,0)%26 for j in range(n)] for r in rows]
    ok2,r2=solve_mod(R,[0]*len(R),2);ok13,r13=solve_mod(R,[0]*len(R),13)
    # homogeneous system: solutions = nullspace. injectivity: is x(a)-x(b) forced to 0 (mod 26) for a!=b ?
    viol=[]
    names=list(var)
    for side in ('x','y'):
        ls=[nm for nm in names if nm[0]==side]
        for i1 in range(len(ls)):
            for i2 in range(i1+1,len(ls)):
                e=[0]*n;e[var[ls[i1]]]=1;e[var[ls[i2]]]=-1
                forced=True
                for m in (2,13):
                    _,ra=solve_mod(R,[0]*len(R),m);_,rb=solve_mod(R+[[x%m for x in e]],[0]*(len(R)+1),m)
                    if rb!=ra: forced=False
                if forced: viol.append((ls[i1],ls[i2]))
    return n,r13,viol
o=dial((86,82,89,92),'L');M=np.argsort(o)
for mode in ('III','IV'):
    print(mode,quag(M,CT0,8,mode))
