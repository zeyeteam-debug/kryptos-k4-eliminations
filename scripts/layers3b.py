import numpy as np,json
from layers3 import WORDS,colidx,KA,ST
from stack import CT0,ANCH
def search_QQT(ct,anch,words,AL):
    """Q(w1)->Q(w2)->T(w) : key indexed by PT position; and Q^3 with no T"""
    n=len(ct);AN=np.array(sorted(anch));P=np.array([AL.index(anch[i]) for i in AN]);C=np.array([AL.index(c) for c in ct])
    W=[(w,np.array([AL.index(c) for c in w])) for w in words]
    K=np.stack([v[AN%len(v)] for _,v in W])
    hits=[]
    Topts=[('none',np.arange(n))]+[(t+w,(np.argsort(colidx(n,w)) if t=='T:' else colidx(n,w))) for w,_ in W if len(w)>1 for t in ('T:','Tinv:')]
    for tname,pos in Topts:
        need=(C[pos[AN]]-P)%26
        for s1 in (1,-1):
            for s2 in (1,-1):
                d={}
                for j in range(len(W)): d.setdefault(((s2*K[j])%26).tobytes(),[]).append(j)
                lhs=(need[None,:]-s1*K)%26
                for i in range(len(W)):
                    for j in d.get(lhs[i].tobytes(),[]): hits.append((tname,W[i][0],s1,W[j][0],s2))
    # three Q, no transposition
    need=(C[AN]-P)%26;h3=[]
    S2={}
    for i in range(len(W)):
        for j in range(i,len(W)):
            for s in ((1,1),(1,-1),(-1,-1)):
                S2.setdefault(((s[0]*K[i]+s[1]*K[j])%26).tobytes(),[]).append((W[i][0],W[j][0],s))
    for k in range(len(W)):
        for s3 in (1,-1):
            for t in S2.get(((need-s3*K[k])%26).tobytes(),[]): h3.append((W[k][0],s3)+t)
    return hits,h3
for AL,an in ((KA,'KA'),(ST,'ST')):
    h,h3=search_QQT(CT0,ANCH,WORDS,AL);print(an,'Q+Q then T:',len(h),h[:5],'| Q+Q+Q:',len(h3),h3[:5])
