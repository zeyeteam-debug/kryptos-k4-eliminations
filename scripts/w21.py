import numpy as np,random
from stack import CT0,ANCH
from gen import run
from stackT import scan as stackscan
from biglist import WORDS
N=97;W=21;R=(N+W-1)//W
def cells(): return [(i//W,i%W) for i in range(N)]
def route(order_cells):
    # order_cells: list of indices (row-major positions) in reading order -> transposition idx
    idx=[i for i in order_cells if i<N]
    assert len(idx)==N
    m=np.empty(N,dtype=np.int16);m[np.array(idx)]=np.arange(N);return m
def routes():
    out={}
    grid=[[r*W+c for c in range(W)] for r in range(R)]
    def add(name,seq): out[name]=route(seq)
    for rot in range(W):                       # row-major with cyclic column rotation per row start
        add(f'rowrot{rot}',[r*W+(c+rot)%W for r in range(R) for c in range(W)])
    add('cols',[r*W+c for c in range(W) for r in range(R)])
    add('cols_up',[r*W+c for c in range(W) for r in reversed(range(R))])
    add('cols_rl',[r*W+c for c in reversed(range(W)) for r in range(R)])
    add('serp_rows',[r*W+(c if r%2==0 else W-1-c) for r in range(R) for c in range(W)])
    add('serp_cols',[ (r if c%2==0 else R-1-r)*W+c for c in range(W) for r in range(R)])
    for s in range(-R,W):                      # diagonals
        pass
    d=[];
    for k in range(R+W-1):
        for r in range(R):
            c=k-r
            if 0<=c<W: d.append(r*W+c)
    add('diag',d)
    d=[]
    for k in range(R+W-1):
        for r in range(R):
            c=r+k-(R-1)
            if 0<=c<W: d.append(r*W+c)
    add('antidiag',d)
    # spiral clockwise from top-left
    top,bot,left,right=0,R-1,0,W-1;sp=[]
    while top<=bot and left<=right:
        sp+=[top*W+c for c in range(left,right+1)];top+=1
        sp+=[r*W+right for r in range(top,bot+1)];right-=1
        if top<=bot: sp+=[bot*W+c for c in range(right,left-1,-1)];bot-=1
        if left<=right: sp+=[r*W+left for r in range(bot,top-1,-1)];left+=1
    add('spiral',sp);add('spiral_rev',sp[::-1])
    # keyword columnar width 21 from 21-letter phrases built from word list
    phr=[w for w in WORDS if len(w)>=21]
    for a in WORDS:
        for b in WORDS:
            if len(a)+len(b)==21: phr.append(a+b)
    for p in set(phr):
        k=p[:21];order=sorted(range(W),key=lambda i:(k[i],i))
        add('kw:'+k,[r*W+c for c in order for r in range(R)])
    return out
if __name__=='__main__':
    rt=routes();names=list(rt);Ms=[];labs=[]
    for n in names:
        m=rt[n];Ms+=[m,np.argsort(m).astype(np.int16)];labs+=[(n,'enc'),(n,'dec')]
    Ms=np.array(Ms);print('width-21 transpositions:',len(Ms),flush=True)
    h=run(CT0,ANCH,Ms,labs,Pmax=21,RTHR=5,ANYSTRONG=3);print('periodic key (any alphabet / ST / KA, P<=21):',len(h),h[:8],flush=True)
    l=list(CT0);random.Random(5).shuffle(l);print('null:',len(run(''.join(l),ANCH,Ms,labs,Pmax=21,RTHR=5,ANYSTRONG=3)),flush=True)
    h2=stackscan(CT0,Ms,labs,amax=8,bmax=8);print('two stacked keys around width-21 transposition:',len(h2),h2[:8])
