import random,math,pickle,numpy as np
from frac import blocks
from stack import CT0,ANCH
Q,fl=pickle.load(open('quad.pkl','rb'))
QA=np.full(26**4,fl,dtype=np.float32)
for k,v in Q.items():
    if len(k)==4 and k.isalpha(): QA[((ord(k[0])-65)*26+ord(k[1])-65)*676+(ord(k[2])-65)*26+ord(k[3])-65]=v
SYM='ABCDEFGHIJKLMNOPQRSTUVWXYZ+'
def make_dec(ct,p,o):
    B=blocks(p,o);cti=[SYM.index(c) for c in ct]
    def dec(sq):   # sq: list cell->symbol index ; pos[symbol]=cell
        pos=[0]*27
        for cell,s in enumerate(sq): pos[s]=cell
        out=[0]*97
        for blk in B:
            L=len(blk);seq=[]
            for j in blk:
                c=pos[cti[j]];seq+= [c//9,(c//3)%3,c%3]
            for j in range(L):
                out[blk[j]]=sq[seq[j]*9+seq[L+j]*3+seq[2*L+j]]
        return out
    return dec
AN=sorted(ANCH);AV=[SYM.index(ANCH[i]) for i in AN]
def score(o):
    a=np.array(o);ok=a<26
    s=0.0
    aa=np.where(ok,a,0);idx=((aa[:-3]*26+aa[1:-2])*26+aa[2:-1])*26+aa[3:]
    s=float(QA[idx].sum())
    s+=sum(25.0 for i,v in zip(AN,AV) if o[i]==v)
    return s
def sa(ct,p,o,iters=60000,restarts=6,seed=0):
    rng=random.Random(seed);dec=make_dec(ct,p,o);best=(-1e9,None)
    for r in range(restarts):
        sq=list(range(27));rng.shuffle(sq);s=score(dec(sq));T=8.0
        for it in range(iters):
            i,j=rng.sample(range(27),2);sq[i],sq[j]=sq[j],sq[i];t=score(dec(sq))
            if t>=s or rng.random()<math.exp((t-s)/T): s=t
            else: sq[i],sq[j]=sq[j],sq[i]
            T=max(0.2,8.0*(1-it/iters))
        if s>best[0]: best=(s,sq[:])
    o_=dec(best[1]);txt=''.join(SYM[x] for x in o_);anc=sum(o_[i]==v for i,v in zip(AN,AV))
    return round(best[0],1),anc,txt
if __name__=='__main__':
    import sys
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
    sq=list(range(27));random.Random(8).shuffle(sq);pos=[0]*27
    for c,s in enumerate(sq): pos[s]=c
    out=['']*97
    for blk in blocks(13,5):
        L=len(blk);seq=[(pos[SYM.index(pt[j])]//[9,3,1][d])%3 for d in range(3) for j in blk]
        for k in range(L):
            t=seq[3*k:3*k+3];out[blk[k]]=SYM[sq[t[0]*9+t[1]*3+t[2]]]
    ctrl=''.join(out)
    print('control:',sa(ctrl,13,5),flush=True)
    for p,o in ((13,5),(16,0),(16,7),(16,8)):
        print('K4',p,o,sa(CT0,p,o),flush=True)
