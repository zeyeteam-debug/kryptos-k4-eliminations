import numpy as np,pickle,random,math,sys
from stack import CT0,ANCH
Q,fl=pickle.load(open('quad.pkl','rb'))
QA=np.full(26**4,fl,dtype=np.float32)
for k,v in Q.items():
    if len(k)==4 and k.isalpha(): QA[((ord(k[0])-65)*26+ord(k[1])-65)*676+(ord(k[2])-65)*26+ord(k[3])-65]=v
def qs(p): return float(QA[((p[:-3]*26+p[1:-2])*26+p[2:-1])*26+p[3:]].sum())
P=19;N=97
AN=np.array(sorted(ANCH));AP=np.array([ord(ANCH[i])-65 for i in AN]);cls=np.arange(N)%P
def make(ct):
    c=np.array([ord(x)-65 for x in ct])
    def decrypt(PA,CAinv,k):        # PA: index->letter ; CAinv: letter->index
        return PA[(CAinv[c]-k[cls])%26]
    def evaluate(PA,CA,kfree):
        CAinv=np.argsort(CA);PAinv=np.argsort(PA)
        # implied k at anchors
        kimp=(CAinv[c[AN]]-PAinv[AP])%26
        k=kfree.copy();viol=0
        for r in range(P):
            sel=kimp[cls[AN]==r]
            if len(sel):
                vals,cnt=np.unique(sel,return_counts=True);k[r]=vals[cnt.argmax()];viol+=len(sel)-cnt.max()
        p=decrypt(PA,CAinv,k)
        return qs(p)-25.0*viol,viol,p,k
    return evaluate
def sa(ct,iters=150000,seed=0,T0=12.0):
    rng=random.Random(seed);ev=make(ct)
    PA=np.array(rng.sample(range(26),26));CA=np.array(rng.sample(range(26),26));kf=np.array([rng.randrange(26) for _ in range(P)])
    s,v,p,k=ev(PA,CA,kf);best=(s,v,p,PA.copy(),CA.copy(),k.copy())
    for it in range(iters):
        T=T0*(1-it/iters)+0.3
        r=rng.random()
        if r<0.45: A=PA
        elif r<0.9: A=CA
        else: A=None
        if A is not None:
            i,j=rng.sample(range(26),2);A[i],A[j]=A[j],A[i]
        else:
            ri=rng.randrange(P);old=kf[ri];kf[ri]=rng.randrange(26)
        t,v2,p2,k2=ev(PA,CA,kf)
        if t>=s or rng.random()<math.exp((t-s)/T): s,v,p,k=t,v2,p2,k2
        else:
            if A is not None: A[i],A[j]=A[j],A[i]
            else: kf[ri]=old
        if s>best[0]: best=(s,v,p.copy(),PA.copy(),CA.copy(),k.copy())
    return best
def txt(p): return ''.join(chr(65+x) for x in p)
if __name__=='__main__':
    mode=sys.argv[1];runs=int(sys.argv[2]);iters=int(sys.argv[3])
    if mode=='ctrl':
        exec(open('board.py').read().split('# K4 cell')[0])
        src=(K1P+K2P+K3P)[150:247];pt=list(src);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
        rng=random.Random(42);PA=rng.sample(range(26),26);CA=rng.sample(range(26),26);key=[rng.randrange(26) for _ in range(P)]
        PAinv=np.argsort(PA);ct=''.join(chr(65+CA[(PAinv[ord(x)-65]+key[i%P])%26]) for i,x in enumerate(pt))
        target=pt
    else:
        ct=CT0;target=None
    res=[]
    for r in range(runs):
        b=sa(ct,iters,seed=r);t=txt(b[2])
        acc=sum(a==c for a,c in zip(t,target))/97 if target else None
        res.append((round(b[0]/94,3),b[1],t,acc));print(r,round(b[0]/94,3),'viol',b[1],t,('acc %.2f'%acc if acc is not None else ''),flush=True)
