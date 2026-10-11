# ADFGVX-style fractionation with a 3x3x3 cube (27 cells >= 26 letters): pt letter -> 3 trits, stream of 291 trits,
# columnar transposition (width w, any order), regroup into trits triples -> ct letters. Unknown cube.
# Constraint: a ct letter whose 3 source trits all come from crib letters gives coordinate equalities.
import numpy as np,itertools,random
from stack import CT0,ANCH
N=97
def colperm(n,order):
    w=len(order);idx=[]
    for c in order: idx+=list(range(c,n,w))
    return idx            # stream position read k-th
def analyse(ct,w,order,mode):
    idx=colperm(3*N,order)
    src=idx if mode=='read' else list(np.argsort(idx))
    # output trit k comes from input trit src[k]; input trit t = (letter t//3, coord t%3)
    parent={}
    def f(x):
        while parent.setdefault(x,x)!=x: x=parent[x]
        return x
    full=0
    for j in range(N):
        s=[src[3*j+q] for q in range(3)]
        if all(t//3 in ANCH for t in s):
            full+=1
            for q in range(3):
                a=('c',ct[j],q);b=('p',ANCH[s[q]//3],s[q]%3)
                # same alphabet letter in pt and ct share a cube cell: identify ('c',L,q)==('p',L,q)
                parent[f(a)]=f(b)
    # identify letter variables across pt/ct roles
    for L in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
        for q in range(3): parent[f(('c',L,q))]=f(('p',L,q))
    letters=sorted({k[1] for k in list(parent)})
    sig={}
    for L in letters:
        t=tuple(f(('p',L,q)) for q in range(3))
        if t in sig and sig[t]!=L: return full,True
        sig[t]=L
    return full,False
if __name__=='__main__':
    for mode in ('read','write'):
        best=0;contra=0;tot=0
        for w in range(2,8):
            for order in itertools.permutations(range(w)):
                full,c=analyse(CT0,w,order,mode);tot+=1;best=max(best,full);contra+=c
        print(mode,'columnar 2-7 all orders:',tot,'| max ct letters fully determined by cribs:',best,'| orders eliminated by a forced collision:',contra)
