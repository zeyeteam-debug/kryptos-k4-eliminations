import numpy as np,random
from gen import run
from stack import CT0,ANCH
N=97
def order_from(digits):
    return sorted(range(len(digits)),key=lambda i:(digits[i],i))
def colmap(W,order):
    L=N//W;r=N%W;lens=[L+(1 if c<r else 0) for c in range(W)];idx=[]
    for c in order: idx+=[row*W+c for row in range(lens[c])]
    m=np.empty(N,dtype=np.int16);m[np.array(idx)]=np.arange(N);return m
keys={'38576577844':'38576577844','3857657784':'3857657784','385765778844':'385765778844','3857065077844':'3857065077844','77844385765':'77844385765','385765':'385765','778844':'778844','675':'675','11091989':'11091989','03111990':'03111990'}
Ms=[];labs=[]
for nm,k in keys.items():
    d=[int(c) for c in k];W=len(d)
    for o,on in ((order_from(d),'asc'),(order_from([-x for x in d]),'desc'),(list(range(W)),'plain')):
        m=colmap(W,o);Ms+= [m,np.argsort(m).astype(np.int16)];labs+=[(nm,on,'enc'),(nm,on,'dec')]
Ms=np.array(Ms)
h=run(CT0,ANCH,Ms,labs,Pmax=12,RTHR=5,ANYSTRONG=3);print('coordinate/date digits as column orders:',len(Ms),'perms; K4 hits',h)
