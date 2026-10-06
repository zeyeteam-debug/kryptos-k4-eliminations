import numpy as np,json,time
from biglist import WORDS
from layers3 import colidx,search,KA,ST
from layers3b import search_QQT
from stack import CT0,ANCH
def TTQ_fast(ct,anch,words,AL):
    n=len(ct);AN=np.array(sorted(anch));P=np.array([AL.index(anch[i]) for i in AN]);C=np.array([AL.index(c) for c in ct])
    Ws=[w for w in words if len(w)>1]
    KS=np.stack([np.array([AL.index(w[i%len(w)]) for i in range(n)]) for w in Ws])
    inv={w:np.argsort(colidx(n,w)) for w in Ws};hits=[]
    for w1 in Ws:
        a=inv[w1][AN]
        for w2 in Ws:
            q=inv[w2][a];need=(C[q]-P)%26;K=KS[:,q]
            for s in (1,-1):
                ok=np.all(((s*K)%26)==need,axis=1)
                for j in np.nonzero(ok)[0]: hits.append((w1,w2,Ws[j],s))
    return hits
t=time.time()
d=json.load(open('pk_controls.json'))
pk6=open('pk6.txt').read().replace('\n','').strip();p6=d['PK6']['plaintext']
a6={i:p6[i] for i in list(range(21,34))+list(range(63,74))}
print('control PK6 with big list:',[h for h in TTQ_fast(pk6,a6,WORDS+['HANDIWORK','SMITHWORK'],KA) if h[2]=='PORTAL'][:1],round(time.time()-t),'s',flush=True)
for AL,an in ((KA,'KA'),(ST,'ST')):
    h=search(CT0,ANCH,WORDS,AL);print(an,'T+Q+Q:',len(h),h[:5],flush=True)
    h,h3=search_QQT(CT0,ANCH,WORDS,AL);print(an,'Q+Q+T:',len(h),h[:5],'Q+Q+Q:',len(h3),h3[:5],flush=True)
    h=TTQ_fast(CT0,ANCH,WORDS,AL);print(an,'T+T+Q:',len(h),h[:5],round(time.time()-t),'s',flush=True)
