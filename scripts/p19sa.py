# Period-19 Quagmire III/IV with UNKNOWN (arbitrary) alphabets: simulated annealing, quadgram English score + exact crib penalty.
import numpy as np,random,math,sys,glob,os,time
from stack import CT0,ANCH
N=97;P=19
def quadgrams():
    f='quad.npy'
    if os.path.exists(f): return np.load(f)
    s=''
    for fn in sorted(glob.glob('texts/nltk_data/packages/corpora/gutenberg/*.txt')):
        s+=''.join(c for c in open(fn,errors='ignore').read().upper() if 'A'<=c<='Z')
    a=(np.frombuffer(s.encode(),np.uint8)-65).astype(np.int64)
    idx=((a[:-3]*26+a[1:-2])*26+a[2:-1])*26+a[3:]
    cnt=np.bincount(idx,minlength=26**4).astype(float)
    q=np.log10((cnt+0.01)/cnt.sum());np.save(f,q);return q
Q=quadgrams()
AN=np.array(sorted(ANCH));PA=np.array([ord(ANCH[i])-65 for i in AN]);CLS=np.arange(N)%P
def score_text(p):
    idx=((p[:-3]*26+p[1:-2])*26+p[2:-1])*26+p[3:];return Q[idx].sum()
class Model:
    def __init__(s,ct,two,beau):
        s.c=np.array([ord(x)-65 for x in ct]);s.two=two;s.beau=beau
    def decrypt(s,pp,pc,k):
        # pp: position of each letter in pt alphabet (26), pc: in ct alphabet. inv arrays for decode
        ip=np.argsort(pp)
        if not s.beau: idx=(pc[s.c]-k[CLS])%26
        else: idx=(k[CLS]-pc[s.c])%26
        return ip[idx]
    def crib_k(s,pp,pc):
        kc=(pc[s.c[AN]]-pp[PA])%26 if not s.beau else (pc[s.c[AN]]+pp[PA])%26
        return kc
def anneal(ct,two,beau,iters,seed,crib_w=6.0):
    rng=np.random.default_rng(seed);M=Model(ct,two,beau)
    pp=rng.permutation(26);pc=pp.copy() if not two else rng.permutation(26)
    def full(pp,pc):
        kc=M.crib_k(pp,pc);k=np.zeros(P,int);conf=0
        cls=CLS[AN]
        for c in range(P):
            v=kc[cls==c]
            if len(v):
                b=np.bincount(v,minlength=26);k[c]=b.argmax();conf+=len(v)-b.max()
            else: k[c]=rng.integers(26)
        return k,conf
    k,conf=full(pp,pc);free=[c for c in range(P) if not np.any(CLS[AN]==c)]
    def sc(pp,pc,k,conf): return score_text(M.decrypt(pp,pc,k))-crib_w*conf
    cur=sc(pp,pc,k,conf);best=(cur,pp.copy(),pc.copy(),k.copy(),conf);T0=12.0
    for it in range(iters):
        T=T0*(1-it/iters)+0.2
        npp,npc,nk=pp.copy(),pc.copy(),k.copy();r=rng.random()
        if r<0.12 and free:
            nk[rng.choice(free)]=rng.integers(26)
            nconf=conf
        else:
            a,b=rng.choice(26,2,replace=False)
            if (not two) or r<0.56:
                npp[[a,b]]=npp[[b,a]]
                if not two: npc=npp
            else: npc[[a,b]]=npc[[b,a]]
            kk,nconf=full(npp,npc)
            for c in range(P):
                if c not in free: nk[c]=kk[c]
        s=sc(npp,npc,nk,nconf)
        if s>=cur or rng.random()<math.exp((s-cur)/T):
            pp,pc,k,conf,cur=npp,npc,nk,nconf,s
            if cur>best[0]: best=(cur,pp.copy(),pc.copy(),k.copy(),conf)
    cur,pp,pc,k,conf=best;txt=''.join(chr(65+x) for x in M.decrypt(pp,pc,k))
    return cur,conf,txt
def plant(two,beau,seed=3):
    rng=np.random.default_rng(seed)
    src=open('texts/carter_vol1.txt',errors='ignore').read().upper();src=''.join(c for c in src if 'A'<=c<='Z')[50000:50097]
    pt=list(src)
    for i,c in ANCH.items(): pt[i]=c
    pp=rng.permutation(26);pc=pp if not two else rng.permutation(26);k=rng.integers(26,size=P);ipc=np.argsort(pc)
    p=np.array([ord(c)-65 for c in pt])
    idx=(pp[p]+k[CLS])%26 if not beau else (k[CLS]-pp[p])%26
    return ''.join(chr(65+x) for x in ipc[idx]),''.join(pt)
if __name__=='__main__':
    which,two,beau,runs,iters=sys.argv[1],int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5])
    t=time.time()
    if which=='ctrl': ct,truth=plant(two,beau);print('truth score',round(score_text(np.array([ord(c)-65 for c in truth])),1),truth)
    elif which=='K4': ct=CT0
    else: l=list(CT0);random.Random(int(which[4:])).shuffle(l);ct=''.join(l)
    res=[anneal(ct,two,beau,iters,s) for s in range(runs)];res.sort(key=lambda r:-r[0])
    for r in res[:3]: print(which,'two' if two else 'one','beau' if beau else 'vig','score %.1f conflicts %d'%(r[0],r[1]),r[2])
    print('time',round(time.time()-t),flush=True)
